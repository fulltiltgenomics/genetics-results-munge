#!/usr/bin/env python3
"""Pinned roll-up tables for the CELLxGENE work list: native cell type -> roll-up CL term,
native tissue -> tissue_general.

    rollup.py [--listing datasets.json] [--download] [--cache-dir DIR] [--out DIR]

Only the publish step applies these tables; the worker keeps native terms, so changing a
roll-up never re-reads an H5AD. Re-running on the same listing and ontology releases writes
byte-identical files (nothing time-dependent is recorded).

Inputs
- the curation API listing (`--listing`, or fetched with `--download` into `--cache-dir`). The
  work list is the cost model's: human, `primary_cell_count > 0`, methylation assays excluded.
  Only the dataset-level `cell_type` and `tissue` term lists are read
- `rollup/cell_type_targets.tsv`, `rollup/pins.tsv` and `rollup/broad_targets.tsv`, all
  hand-curated and reviewed. Every id must exist in the pinned CL release with exactly the label
  written there, so a typo fails the build instead of silently becoming a target nobody rolls up to
- ontologies, pinned by the installed `cellxgene-ontology-guide` for the listing's schema
  version: CL ancestry from that CL release's `cl-basic.obo` (the guide's own ancestor sets
  mix is_a with part_of and expose no edge types), labels and UBERON ancestry from the guide

Cell type rule
- is_a only, from `cl-basic.obo` (which carries the reasoner-inferred is_a closure). part_of
  and develops_from are not followed: part_of puts cells under anatomical structures rather
  than broader cell types, and develops_from would pool a cell with its precursor
- each term maps to its most specific is_a ancestor-or-self in the target list; equally
  specific candidates are broken by target-file order and the losers are kept in
  `other_candidates` for review
- a pin sends its native term, and any descendant with no target at or below the pinned term,
  to the pin's target
- a term with no target ancestor maps to itself (`no_target`), as does a term absent from the
  release; `unknown` stays `unknown`. Mapping to itself rather than to a catch-all "cell" keeps
  such groups distinct instead of pooling unrelated cells

Summary exclusion (`summary_excluded` in `cell_type_rollup.tsv`; publish.py leaves these
native terms out of the cross-dataset summary and keeps them in the per-dataset grain): a native
term that rolls up to a target listed in `broad_targets.tsv` and is itself an is_a
ancestor-or-self of another target. Such a label is coarser than the target list, so its pooled
row would compete with the finer types it contains. A term with no target (`no_target`) stays
in even when it sits above targets ('neural cell'): the exclusion is limited to the reviewed
umbrellas, and leaving those out did not change the marker panel in the rehearsal. The
curated list is needed because being an ancestor of another target is not enough: endothelial
cell, macrophage and neuron are ancestors of targets too, yet are the resolution most datasets
label at. Matching on the native term rather than the roll-up keeps a type that merely has no
finer target (pro-B cell -> hematopoietic cell) in the summary.

Tissue rule: CELLxGENE Census's own (`TissueMapper.get_high_level_tissue` in the census
builder): the first term of the guide's curated `tissue_general` list that is an
ancestor-or-self in the guide's UBERON closure (is_a + part_of), else the term itself. Census
leaves CL (primary cell culture) and Cellosaurus (cell line) tissues as themselves too, and so
does this; it is safe because `tissue_type` is a grouping key, so those never pool with tissue.

Outputs in `--out`: `cell_type_rollup.tsv`, `tissue_rollup.tsv`, `metadata.json` (ontology
releases, input checksums, rule text, per-method counts).
"""
import argparse
import csv
import datetime
import hashlib
import importlib.metadata
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sumstat_utils import fetch

HERE = Path(__file__).resolve().parent / "rollup"
LISTING_URL = "https://api.cellxgene.cziscience.com/curation/v1/datasets"
CL_OBO_URL = "https://github.com/obophenotype/cell-ontology/releases/download/{version}/cl-basic.obo"
HUMAN = "NCBITaxon:9606"
UNKNOWN = "unknown"


def work_list(datasets):
    human = [x for x in datasets if any(o["ontology_term_id"] == HUMAN for o in x["organism"])]
    return [x for x in human if (x.get("primary_cell_count") or 0) > 0
            and not any("snmC" in a["label"] or "methyl" in a["label"].lower() for a in x["assay"])]


def load_obo(path):
    parents, labels, obsolete, version = {}, {}, set(), None
    cur = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("data-version: "):
                version = line.split()[1]
            elif line.startswith("["):
                cur = None
            elif line.startswith("id: ") and cur is None:
                cur = line[4:]
                parents.setdefault(cur, [])
            elif cur and line.startswith("name: "):
                labels[cur] = line[6:]
            elif cur and line.startswith("is_a: "):
                parents[cur].append(line[6:].split()[0])
            elif cur and line == "is_obsolete: true":
                obsolete.add(cur)
    return parents, labels, obsolete, version


def _read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader((ln for ln in f if not ln.startswith("#")), delimiter="\t"))


def read_targets(path, labels):
    rows = _read_tsv(path)
    ids = [r["term_id"] for r in rows]
    dup = {i for i in ids if ids.count(i) > 1}
    bad = [f"{r['term_id']} {r['label']!r} (release: {labels.get(r['term_id'])!r})"
           for r in rows if labels.get(r["term_id"]) != r["label"]]
    if dup or bad:
        raise SystemExit(f"{path}: duplicate {sorted(dup)}; id/label not in the release: {bad}")
    return ids


def read_broad(path, labels, ct):
    rows = _read_tsv(path)
    for r in rows:
        t = r["term_id"]
        problems = [p for p, bad in (
            ("empty reason", not r["reason"].strip()),
            (f"label is {labels.get(t)!r} in the release", labels.get(t) != r["label"]),
            ("not in the target list", t not in ct.rank),
            ("not an is_a ancestor of another target", not ct.finer_targets(t, t)),
            ("listed twice", [x["term_id"] for x in rows].count(t) > 1)) if bad]
        if problems:
            raise SystemExit(f"{path}: {t}: {'; '.join(problems)}")
    return [r["term_id"] for r in rows]


def read_pins(path, labels, targets):
    pins = {}
    for r in _read_tsv(path):
        n, t = r["native_term_id"], r["target_term_id"]
        problems = [p for p, bad in (
            ("empty reason", not r["reason"].strip()),
            (f"native label is {labels.get(n)!r}", labels.get(n) != r["native_label"]),
            (f"target label is {labels.get(t)!r}", labels.get(t) != r["target_label"]),
            ("target is not in the target list", t not in targets),
            ("native term is itself a target", n in targets),
            ("pinned twice", n in pins)) if bad]
        if problems:
            raise SystemExit(f"{path}: pin {n} -> {t}: {'; '.join(problems)}")
        pins[n] = t
    return pins


class CellTypeRollup:
    def __init__(self, parents, labels, obsolete, targets, pins):
        self.parents, self.labels, self.obsolete = parents, labels, obsolete
        self.targets, self.rank = targets, {t: i for i, t in enumerate(targets)}
        self.pins = pins
        self._anc = {}

    def ancestors(self, term):
        """is_a ancestors-or-self."""
        if term not in self._anc:
            out = {term}
            for p in self.parents.get(term, ()):
                out |= self.ancestors(p)
            self._anc[term] = out
        return self._anc[term]

    def finer_targets(self, term, own=None):
        """Targets other than `own` that `term` is an is_a ancestor-or-self of."""
        return [t for t in self.targets if t != own and term in self.ancestors(t)]

    def summary_excluded(self, term, rollup, method, broad):
        return (method in ("self", "ancestor", "pin", "pin_descendant") and rollup in broad
                and bool(self.finer_targets(term, rollup)))

    def _most_specific(self, cands):
        best = [t for t in cands if not any(o != t and t in self.ancestors(o) for o in cands)]
        return sorted(best, key=lambda t: self.rank.get(t, -1))

    def map(self, term):
        """(roll-up term, method, other equally specific candidates)."""
        if term == UNKNOWN:
            return UNKNOWN, "unknown", []
        # obsolete terms are not followed through replaced_by/consider: the work list has none
        # today, and a schema-valid dataset cannot carry one, so a hit here means the listing
        # and the pinned CL release have drifted apart
        if term not in self.parents or term in self.obsolete:
            return term, "not_in_ontology", []
        anc = self.ancestors(term)
        cands = [t for t in self.targets if t in anc]
        pinned = self._most_specific([p for p in self.pins if p in anc])
        if len(pinned) > 1:
            raise SystemExit(f"{term}: falls under unrelated pins {pinned}; add a pin for it")
        if pinned:
            p = pinned[0]
            # only targets at or below the pinned term survive: a descendant that reaches its own
            # closer target through a different parent loses it to the pin
            below = [t for t in cands if p in self.ancestors(t)]
            if not below:
                return self.pins[p], "pin" if term == p else "pin_descendant", []
            cands = below
        if not cands:
            return term, "no_target", []
        best = self._most_specific(cands)
        return best[0], "self" if best[0] == term else "ancestor", best[1:]


def tissue_general(term, parser, general):
    """Census builder semantics: first listed tissue_general that is an ancestor-or-self."""
    if not term.startswith("UBERON:"):
        return term, "cell_culture_self"
    try:
        hits = parser.get_high_level_terms(term, general)
    except (KeyError, ValueError):
        return term, "not_in_ontology"
    return (hits[0], "tissue_general") if hits else (term, "no_tissue_general")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def build(listing, cl_obo, out, targets_path=HERE / "cell_type_targets.tsv", pins_path=HERE / "pins.tsv",
          broad_path=None):
    # beside the targets by default, so a test's synthetic targets get their own list
    broad_path = broad_path or Path(targets_path).parent / "broad_targets.tsv"
    from cellxgene_ontology_guide.curated_ontology_term_lists import (
        get_curated_ontology_term_list,
    )
    from cellxgene_ontology_guide.entities import CuratedOntologyTermList
    from cellxgene_ontology_guide.ontology_parser import OntologyParser

    with open(listing) as f:
        datasets = json.load(f)
    work = work_list(datasets)
    schemas = sorted({x["schema_version"] for x in work})
    if len(schemas) != 1:
        raise SystemExit(f"work list spans schema versions {schemas}; pin one ontology set per build")
    parser = OntologyParser(schemas[0])
    onto = parser.cxg_schema.supported_ontologies

    parents, labels, obsolete, obo_version = load_obo(cl_obo)
    if obo_version != f"cl/releases/{onto['CL']['version'].lstrip('v')}/cl-basic.owl":
        raise SystemExit(f"{cl_obo} is {obo_version}, schema {schemas[0]} pins CL {onto['CL']['version']}")
    targets = read_targets(targets_path, labels)
    pins = read_pins(pins_path, labels, targets)
    ct = CellTypeRollup(parents, labels, obsolete, targets, pins)
    broad = read_broad(broad_path, labels, ct)

    listed_ct = {c["ontology_term_id"]: c["label"] for x in work for c in x["cell_type"]}
    ct_rows = []
    for term in sorted(listed_ct):
        tgt, method, others = ct.map(term)
        ct_rows.append([term, labels.get(term, listed_ct[term]), tgt,
                        labels.get(tgt, listed_ct.get(tgt, tgt)), method, ";".join(others),
                        str(ct.summary_excluded(term, tgt, method, broad)).lower()])

    listed_ts, ts_types = {}, {}
    for x in work:
        for t in x["tissue"]:
            listed_ts[t["ontology_term_id"]] = t["label"]
            ts_types.setdefault(t["ontology_term_id"], set()).add(t["tissue_type"])
    general = get_curated_ontology_term_list(CuratedOntologyTermList.TISSUE_GENERAL)

    def label(term):
        try:
            return parser.get_term_label(term) or listed_ts.get(term, term)
        except (KeyError, ValueError):
            return listed_ts.get(term, term)

    ts_rows = []
    for term in sorted(listed_ts):
        tg, method = tissue_general(term, parser, general)
        ts_rows.append([term, label(term), ",".join(sorted(ts_types[term])), tg, label(tg), method])

    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "cell_type_rollup.tsv", ["cell_type_ontology_term_id", "cell_type", "rollup_term_id",
              "rollup_label", "method", "other_candidates", "summary_excluded"], ct_rows)
    write_tsv(out / "tissue_rollup.tsv", ["tissue_ontology_term_id", "tissue", "tissue_types",
              "tissue_general_ontology_term_id", "tissue_general", "method"], ts_rows)

    def counts(rows, col):
        c = {}
        for r in rows:
            c[r[col]] = c.get(r[col], 0) + 1
        return dict(sorted(c.items()))

    meta = {
        "cellxgene_schema_version": schemas[0],
        "cellxgene_ontology_guide": importlib.metadata.version("cellxgene-ontology-guide"),
        "ontologies": {k: onto[k]["version"] for k in ("CL", "UBERON", "CVCL")},
        "cl_obo": {"url": CL_OBO_URL.format(version=onto["CL"]["version"]), "data_version": obo_version,
                   "sha256": sha256(cl_obo)},
        "listing": {"url": LISTING_URL, "sha256": sha256(listing), "n_datasets": len(datasets),
                    "n_work_list": len(work),
                    "filter": "human, primary_cell_count > 0, no methylation (snmC*/methyl*) assay"},
        "inputs_sha256": {p.name: sha256(p) for p in map(Path, (targets_path, pins_path, broad_path))},
        "n_targets": len(targets),
        "n_pins": len(pins),
        "broad_targets": broad,
        "cell_type_rule": "most specific is_a ancestor-or-self in cell_type_targets.tsv (part_of and "
                          "develops_from not followed); ties by target-file order; pins.tsv overrides; "
                          "no target -> self",
        "summary_excluded_rule": "native term that rolls up to a target in broad_targets.tsv and is an is_a "
                                 "ancestor-or-self of another target; no_target terms are never excluded; "
                                 "left out of the cross-dataset summary only",
        "n_summary_excluded": sum(r[6] == "true" for r in ct_rows),
        "tissue_rule": "CELLxGENE Census TissueMapper: first tissue_general list entry that is an "
                       "ancestor-or-self (guide closure, is_a + part_of); else self; CL and CVCL "
                       "tissues -> self",
        "cell_type_methods": counts(ct_rows, 4),
        "tissue_methods": counts(ts_rows, 5),
    }
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    return work, ct_rows, ts_rows


def coverage(work, ct_rows):
    """Native terms and (estimated) primary cells landing on a target other than themselves.

    The listing gives per-dataset term lists, not per-term counts, so cells are estimated by
    splitting each dataset's primary_cell_count evenly over its listed cell types.
    """
    rolled = {r[0] for r in ct_rows if r[4] in ("ancestor", "pin", "pin_descendant")}
    on_target = rolled | {r[0] for r in ct_rows if r[4] == "self"}
    cells = {}
    for x in work:
        terms = {c["ontology_term_id"] for c in x["cell_type"]}
        for t in terms:
            cells[t] = cells.get(t, 0.0) + x["primary_cell_count"] / len(terms)
    total = sum(cells.values())
    return {"n_native": len(ct_rows),
            "terms_rolled_pct": 100 * len(rolled) / len(ct_rows),
            "terms_on_target_pct": 100 * len(on_target) / len(ct_rows),
            "cells_rolled_pct_est": 100 * sum(cells[t] for t in rolled) / total,
            "cells_on_target_pct_est": 100 * sum(cells[t] for t in on_target) / total}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--listing", type=Path, help="curation API datasets JSON (default: fetch into --cache-dir)")
    ap.add_argument("--cache-dir", type=Path, default=Path(os.environ.get("CELLXGENE_CACHE", "cellxgene-cache")))
    ap.add_argument("--download", action="store_true", help="allow fetching the listing and cl-basic.obo")
    ap.add_argument("--out", type=Path, default=HERE)
    a = ap.parse_args()

    # dated so a cached listing is never mistaken for this week's
    today = datetime.date.today().isoformat()
    listing = a.listing or fetch(LISTING_URL, a.cache_dir / f"curation-datasets-{today}.json", download=a.download)
    from cellxgene_ontology_guide.ontology_parser import OntologyParser
    with open(listing) as f:
        schema = work_list(json.load(f))[0]["schema_version"]
    cl_version = OntologyParser(schema).cxg_schema.supported_ontologies["CL"]["version"]
    cl_obo = fetch(CL_OBO_URL.format(version=cl_version), a.cache_dir / f"cl-basic-{cl_version}.obo",
                   download=a.download)
    work, ct_rows, _ = build(listing, cl_obo, a.out)
    print(json.dumps({k: round(v, 1) for k, v in coverage(work, ct_rows).items()}), file=sys.stderr)


if __name__ == "__main__":
    main()
