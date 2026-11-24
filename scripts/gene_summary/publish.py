#!/usr/bin/env python3
"""Turn checked gene summary records into the release folder results-api serves gene summaries
from (`/gene_summary/{gene}`): one `<SYMBOL>.json` per gene and a `manifest.json`.

A record cites evidence by bundle row id, which means nothing outside the bundle. The served
file carries, for each cited id, a one-line description of the row, so a reader can find it
in the tab the finding names. Records that failed the check are not published.
"""

import argparse
import datetime
import json
import sys
from pathlib import Path

MANIFEST = "manifest.json"


def rows_by_id(node, path: tuple = (), out: dict | None = None) -> dict:
    """Every bundle row keyed by id, with the section it sits in."""
    out = {} if out is None else out
    if isinstance(node, dict):
        if isinstance(node.get("id"), str):
            out[node["id"]] = (path[0], node)
        for k, v in node.items():
            rows_by_id(v, (*path, k), out)
    elif isinstance(node, list):
        for v in node:
            rows_by_id(v, path, out)
    return out


def join(*parts) -> str:
    return " · ".join(str(p) for p in parts if p not in (None, "", []))


def describe(section: str, r: dict) -> str:
    if section == "burden":
        return join(r["dataset"], r["trait"], r["mask"], f"p {r['p']}")
    if section == "exome_variants":
        return join(
            r["variant"], r["annotation"], f"{r['n_significant_traits']} exome-wide significant traits"
        )
    if section == "gwas_credible_sets":
        if "max_pip" in r:
            return join(r["variant"], r["consequence"], f"{r['n_traits']} traits", f"max PIP {r['max_pip']}")
        return join(r["dataset"], r["trait"], r["kind"], f"lead {r['lead_variant']}", f"p {r['p']}")
    if section == "qtl":
        return join(
            r["dataset"], r["tissue_or_cell_type"], r["molecular_trait"],
            f"lead {r['lead_variant']}", f"p {r['p']}",
        )
    if section in ("colocalization_cis", "colocalization_trans"):
        best = r["best"]
        return join(
            r["gwas_dataset"], r["gwas_trait"],
            f"{'cis' if section.endswith('cis') else 'trans'} {best['qtl_type']}",
            best.get("tissue_or_cell_type"), f"PP.H4 {best['pp_h4']}",
        )
    if section == "expression_gtex":
        if "tissue" in r:
            return join("GTEx", r["tissue"], f"median TPM {r['median_tpm']}")
        return join("GTEx", f"{r['n_tissues']} tissues", f"median TPM across tissues {r['median_across_tissues']}")
    if section == "expression_hpa":
        return join("Human Protein Atlas", r["tissue"], r["cell_type"], r["level"])
    if section == "expression_cell_types":
        if "cell_type" in r:
            return join("CELLxGENE", r["cell_type"], r["tissue"], f"{r['percent_cells_expressing']}% of cells")
        return join("CELLxGENE", f"{r['n_groups']} cell type and tissue groups", f"tau {r['tau']}")
    if section == "mendelian":
        return join(r["disease"], r["best_gencc_classification"], ", ".join(r["inheritance"] or []))
    if section == "constraint":
        return join("gnomAD", f"LOEUF {r['loeuf']}", f"pLI {r['pLI']}")
    raise ValueError(f"no description for a row of section {section}")


def served(record: dict, bundle: dict, generated: str) -> dict:
    rows = rows_by_id(bundle)
    s = record["summary"]
    return {
        "gene": record["gene"],
        "generated": generated,
        "model": record["model"],
        "headline": s["headline"],
        "overview": s["overview"],
        "findings": [
            {
                "statement": f["statement"],
                "tab": f["tab"],
                "evidence": [{"id": e, "label": describe(*rows[e])} for e in f["evidence"]],
            }
            for f in s["findings"]
        ],
        "gaps": s["gaps"],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bundles", required=True, help="directory of <SYMBOL>.json bundles")
    ap.add_argument("--records", required=True, help="directory of <SYMBOL>.json summary records")
    ap.add_argument("--out", required=True, help="release directory, to be copied to the bucket folder")
    ap.add_argument(
        "--date", default=datetime.date.today().isoformat(), help="the date the page shows"
    )
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    skipped = 0
    for path in sorted(Path(args.records).glob("*.json")):
        record = json.loads(path.read_text())
        bundle = json.loads((Path(args.bundles) / path.name).read_text())
        # a summary of other data than the bundle beside it would get the wrong row descriptions
        if not record["valid"] or record["bundle_sha256"] != bundle["bundle_sha256"]:
            skipped += 1
            print(f"{path.stem}\tskipped: {'invalid' if not record['valid'] else 'bundle changed'}")
            continue
        (out / path.name).write_text(json.dumps(served(record, bundle, args.date)))
    # results-api refuses to start on a configured release without one; it lists what the
    # folder holds, earlier runs included
    genes = sorted(p.stem for p in out.glob("*.json") if p.name != MANIFEST)
    (out / MANIFEST).write_text(json.dumps({"written": args.date, "n_genes": len(genes), "genes": genes}))
    return 1 if skipped else 0


if __name__ == "__main__":
    sys.exit(main())
