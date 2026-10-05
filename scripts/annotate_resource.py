#!/usr/bin/env python3
"""Stamp gnomAD consequence onto every object of one served credible-set resource.

A served resource is a tree of objects under one prefix: the combined variant-sorted file,
perhaps a gene-indexed QTL copy, one file per phenotype, their tabix indexes and the stats
derived from them. This writes the whole tree again under a NEW prefix with `most_severe`
and `gene_most_severe` taken from a consequence file, keeping every relative name, so that
an API profile can be pointed at the new prefix. `annotate_consequence.py` does the stamping
and the row-by-row verification of each file; this script decides what each object is, runs
it, regenerates what is derived from the annotation, and proves the result before and after
it is uploaded.

Source and destination are each a `gs://bucket/prefix/` or a local directory.

Classification
--------------
Every object under the source prefix gets exactly one class, by the first rule that matches
its name relative to the prefix. An object no rule matches stops the run before anything is
read or written, and so does `--dry-run`.

    placeholder     empty object whose name ends in `/`              skipped
    report          a previous run's own report (see Outputs)        skipped
    combined        the name given as --combined                     stamped, merge mode
    qtl             the name given as --qtl                          stamped, lookup mode
    qtl_unmapped    the --qtl name with `.tsv.gz` -> `.unmapped.tsv` stamped, lookup mode
    stats_aggregate the name given as --stats-file                   regenerated
    stats_json      `*.stats.json`                                   regenerated
    per_phenotype   --per-phenotype-prefix + NAME + --per-phenotype-suffix,
                    NAME holding no `/`                              stamped, merge mode
    index           `.tbi` / `.csi` of a stamped object              regenerated
    skip            matches a --skip glob                            skipped
    copy            `*.log`, or matches a --copy glob                copied byte-identical

The prefix and suffix of the per-phenotype rule are the `prefix` (relative to the source) and
`suffix_95` of the resource's API profile entry, which is how the API builds those names.
The rows the QTL build could not place on a gene are credible-set rows carrying the two
annotation columns, which is why that file is stamped and not copied; lookup mode because
nothing promises their order. A file sent to `copy` whose header names both annotation
columns stops the run: it would carry the old annotation into the new tree.

Stats
-----
`credible_set_stats.tsv` and the per-phenotype `*.stats.json` count coding and LoF credible
sets from `most_severe`, so they are recomputed with `credible_set_stats.py` from the
stamped per-phenotype files, one record per trait of each file. No producer's row order or
file naming is rebuilt from rules. The same computation is first run on the ORIGINAL
per-phenotype files, and every row of the original aggregate and every original stats.json
must be reproduced byte for byte by one of those records; the stamped record of the same
file and trait then takes its place, in the original's row order and under the original's
name. A resource whose stats are not reproduced this way stops the run, as does a record
with a trait that the aggregate does not hold.

Rows of the aggregate that are identical in every field (the same gene with the same counts
in two studies of a QTL resource, whose rows do not name the study) cannot be told apart, so
they are paired with their files in name order.

Outputs
-------
Nothing is written at or under the source, and the run refuses a destination that holds an
object it would not itself write, or one whose content differs from what it is about to
write. Work happens in `--staging`: `original/` (the source tree), `stamped/` (stamped files,
their indexes, regenerated stats), `copied/`, `proof/` (stats regenerated from the original
rows), `state/` (one record per stamped and verified file) and `report/`. Upload starts only
when every object of the resource has been produced and verified locally, so the destination
never holds a partly written or unverified object. A rerun with the same arguments skips
files that have a state record and objects the destination already holds with the same
checksum; a staging directory is bound to its arguments and refuses different ones.

Stamped and regenerated objects of a `gs://` destination carry custom metadata:
`consequence-source`, `consequence-version`, `consequence-file`, `consequence-crc32c` and
`stamped-from`. A local destination cannot hold metadata.

`annotation_report.tsv` (one row per source object) and `annotation_summary.txt` are written
to `report/` and, when every check passed, uploaded beside the stamped objects. They hold no
timings so that a rerun writes the same bytes; timings go to `timing.json` in the staging
directory. The exit status is non-zero when any check failed.

Usage
-----
    scripts/annotate_resource.sh \\
        --source gs://bucket/credible_sets/pgc_scz_finemap/2022/ \\
        --dest gs://bucket/credible_sets/pgc_scz_finemap/2022_vep115/ \\
        --consequence gnomad.v4.1.1.consequence.tsv.bgz --consequence-version 4.1.1-vep115 \\
        --staging /mnt/disks/data/stamp/pgc_scz_finemap \\
        --combined PGC_SCZ_2022_credible_sets.tsv.gz \\
        --per-phenotype-prefix individual/ --per-phenotype-suffix .FINEMAP.munged.tsv
"""

import argparse
import base64
import fnmatch
import gzip
import hashlib
import json
import os
import resource
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict, deque
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
STAMPER = HERE / "annotate_consequence.py"
ANNOTATION = ("most_severe", "gene_most_severe")
INDEX_EXTS = (".tbi", ".csi")
REPORT_TSV, REPORT_SUMMARY = "annotation_report.tsv", "annotation_summary.txt"
STAMP_MODE = {"combined": "merge", "qtl": "lookup", "qtl_unmapped": "lookup", "per_phenotype": "merge"}
REGENERATED = ("index", "stats_aggregate", "stats_json")
SKIPPED = ("placeholder", "report", "skip")
# files stamped one at a time with every thread, since a single one can outweigh the rest
SERIAL = ("combined", "qtl", "qtl_unmapped")
STATS_COLUMNS = {"trait", "trait_original", "dataset", "#dataset", "data_type", "pip", "beta", "aaf",
                 "cs_id", "most_severe"}
# gcloud takes seconds to start, so local files are hashed many to a call
HASH_BATCH = 2000


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def die(msg: str) -> None:
    sys.exit(f"annotate_resource: {msg}")


def run(cmd: list, **kw) -> subprocess.CompletedProcess:
    p = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if p.returncode != 0:
        raise RuntimeError(f"{' '.join(map(str, cmd[:6]))} ... failed:\n{p.stderr[-2000:]}")
    return p


def tree(root: Path) -> list[str]:
    return sorted(str(p.relative_to(root)) for p in root.rglob("*") if p.is_file())


class LocalStore:
    """A directory standing in for a bucket prefix."""

    digest_name = "md5"
    holds_metadata = False

    def __init__(self, path: str):
        self.root = Path(os.path.realpath(path))
        self.url = f"{self.root}/"

    @staticmethod
    def file_digests(paths: list) -> dict:
        out = {}
        for path in paths:
            h = hashlib.md5()
            with open(path, "rb") as f:
                while chunk := f.read(1 << 20):
                    h.update(chunk)
            out[str(path)] = h.hexdigest()
        return out

    metadata_digest = staticmethod(str)

    def list(self) -> dict:
        if not self.root.is_dir():
            return {}
        rels = tree(self.root)
        digests = self.file_digests([self.root / r for r in rels])
        return {r: {"size": (self.root / r).stat().st_size, "digest": digests[str(self.root / r)], "metadata": {}}
                for r in rels}

    def fetch(self, into: Path) -> None:
        for rel in tree(self.root):
            target = into / rel
            if not target.exists() or target.stat().st_size != (self.root / rel).stat().st_size:
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(self.root / rel, target)

    def push(self, source: Path, metadata: dict) -> None:
        for rel in tree(source):
            target = self.root / rel
            if target.exists():
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            # under a scratch name first, so an interrupted copy is never taken for an object
            part = target.with_name(target.name + ".part")
            shutil.copyfile(source / rel, part)
            os.replace(part, target)


class GcsStore:
    digest_name = "crc32c"
    holds_metadata = True

    def __init__(self, url: str):
        self.url = url.rstrip("/") + "/"
        self.bucket, _, self.prefix = self.url[len("gs://"):].partition("/")
        if not self.bucket or not self.prefix:
            die(f"{url}: expected gs://bucket/prefix/")

    @staticmethod
    def file_digests(paths: list) -> dict:
        out = {}
        paths = [str(p) for p in paths]
        for i in range(0, len(paths), HASH_BATCH):
            p = run(["gcloud", "storage", "hash", "--skip-md5", "--format=json", *paths[i:i + HASH_BATCH]])
            out.update({h["url"]: h["crc32c_hash"] for h in json.loads(p.stdout)})
        return out

    @staticmethod
    def metadata_digest(digest: str) -> str:
        # hex, because the padding of base64 is the separator of --custom-metadata
        return base64.b64decode(digest).hex()

    def list(self) -> dict:
        p = subprocess.run(["gcloud", "storage", "objects", "list", f"{self.url}**", "--format=json"],
                           capture_output=True, text=True)
        if p.returncode != 0:
            if "matched no objects" in p.stderr:
                return {}
            die(f"listing {self.url} failed:\n{p.stderr[-2000:]}")
        return {o["name"][len(self.prefix):]: {"size": int(o["size"]), "digest": o["crc32c_hash"],
                                               "metadata": o.get("custom_fields") or {}}
                for o in json.loads(p.stdout)}

    def fetch(self, into: Path) -> None:
        run(["gcloud", "storage", "rsync", "--recursive", self.url, str(into)])

    def push(self, source: Path, metadata: dict) -> None:
        cmd = ["gcloud", "storage", "rsync", "--recursive", "--no-clobber"]
        if metadata:
            cmd.append("--custom-metadata=" + ",".join(f"{k}={v}" for k, v in metadata.items()))
        run([*cmd, str(source), self.url])


def store(location: str):
    return GcsStore(location) if location.startswith("gs://") else LocalStore(location)


def classify(listing: dict, rules) -> dict:
    """{relative name: class}. Stops the run, naming every offender, when an object matches
    no rule or the classes do not add up to a resource this can stamp."""
    unmapped = rules.qtl.removesuffix(".tsv.gz") + ".unmapped.tsv" if rules.qtl else None
    classes, problems = {}, []
    for rel, obj in listing.items():
        middle = rel[len(rules.per_phenotype_prefix or ""):len(rel) - len(rules.per_phenotype_suffix or "")]
        if rel.endswith("/") and obj["size"] == 0:
            cls = "placeholder"
        elif rel in (REPORT_TSV, REPORT_SUMMARY):
            cls = "report"
        elif rel == rules.combined:
            cls = "combined"
        elif rel == rules.qtl:
            cls = "qtl"
        elif rel == unmapped:
            cls = "qtl_unmapped"
        elif rel == rules.stats_file:
            cls = "stats_aggregate"
        elif rel.endswith(".stats.json"):
            cls = "stats_json"
        elif (rules.per_phenotype_suffix and rel.startswith(rules.per_phenotype_prefix)
              and rel.endswith(rules.per_phenotype_suffix) and middle and "/" not in middle):
            cls = "per_phenotype"
        elif rel.endswith(INDEX_EXTS):
            cls = "index"
        elif any(fnmatch.fnmatchcase(rel, g) for g in rules.skip):
            cls = "skip"
        elif any(fnmatch.fnmatchcase(rel, g) for g in ["*.log", *rules.copy]):
            cls = "copy"
        else:
            problems.append(f"{rel}: matches no rule (name it with --copy or --skip if it is neither data nor stats)")
            continue
        classes[rel] = cls

    for flag, name in (("--combined", rules.combined), ("--qtl", rules.qtl)):
        if name and name not in listing:
            problems.append(f"{name}: given as {flag} but not under the source")
    for rel, cls in classes.items():
        if cls == "index" and classes.get(os.path.splitext(rel)[0]) not in STAMP_MODE:
            problems.append(f"{rel}: index of an object that is not stamped")
    found = Counter(classes.values())
    if rules.per_phenotype_suffix and not found["per_phenotype"]:
        problems.append("the per-phenotype prefix and suffix match no object")
    if (found["stats_aggregate"] or found["stats_json"]) and not found["per_phenotype"]:
        problems.append("stats objects but no per-phenotype files to regenerate them from")
    if problems:
        die("cannot classify the source:\n  " + "\n  ".join(sorted(problems)))
    return classes


def action(cls: str) -> str:
    if cls in STAMP_MODE:
        return f"stamp-{STAMP_MODE[cls]}"
    return "regenerate" if cls in REGENERATED else "skip" if cls in SKIPPED else "copy"


def stamper(*argv: str) -> tuple[int, dict, str]:
    """Exit status, the JSON summary on the last stdout line, and stderr."""
    p = subprocess.run([sys.executable, str(STAMPER), *argv], capture_output=True, text=True)
    lines = p.stdout.strip().splitlines()
    try:
        return p.returncode, json.loads(lines[-1]) if lines else {}, p.stderr
    except json.JSONDecodeError:
        return p.returncode or 1, {}, p.stderr


def stats_records(path: Path) -> list[dict]:
    """One credible_set_stats record per trait of a per-phenotype file, in order of first
    appearance. Every column is read as text so that a trait or a credible-set id made of
    digits stays the string the file holds."""
    import polars as pl
    from credible_set_stats import calculate_stats

    with open(path, "rb") as f:
        opener = gzip.open if f.read(2) == b"\x1f\x8b" else open
    with opener(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    df = pl.read_csv(path, separator="\t", null_values=["NA"], infer_schema_length=0,
                     columns=[c for c in header if c in STATS_COLUMNS])
    return [asdict(calculate_stats(part)) for part in df.partition_by("trait", maintain_order=True)]


def stamp_one(job: dict) -> dict:
    """Stamp one file, verify it against its original and leave a state record. Runs in a
    worker process; a failure comes back in the record, not as an exception, so one bad
    file does not hide the state of the others."""
    rel, original, stamped = job["rel"], Path(job["original"]), Path(job["stamped"])
    record = {"rel": rel, "class": job["class"], "source": job["source"]}
    try:
        stamped.parent.mkdir(parents=True, exist_ok=True)
        code, summary, err = stamper("--consequence", job["consequence"], "--mode", STAMP_MODE[job["class"]],
                                     "--input", str(original), "--output", str(stamped),
                                     "--threads", str(job["threads"]), *job["stamper_args"])
        if code != 0:
            raise RuntimeError(f"stamping failed: {err[-1500:]}")
        code, verify, err = stamper("--verify", str(original), str(stamped), "--threads", str(job["threads"]))
        if not verify:
            raise RuntimeError(f"verification did not run: {err[-1500:]}")
        record.update(stamp=summary, verify=verify, size=stamped.stat().st_size)
        record["index"] = next((e for e in INDEX_EXTS if os.path.exists(f"{stamped}{e}")), None)
        # the stamper's own account and the independent comparison must tell the same story
        agree = all(summary[c] == verify[c] for c in ANNOTATION) and summary["rows"] == verify["rows_stamped"]
        if not (verify["ok"] and agree):
            raise RuntimeError("the stamped file does not verify against its original")
        if job["class"] == "per_phenotype" and job["stats"]:
            record["stats"] = {"original": stats_records(original), "stamped": stats_records(stamped)}
        state = Path(job["state"])
        state.parent.mkdir(parents=True, exist_ok=True)
        part = state.with_name(state.name + ".part")
        part.write_text(json.dumps(record))
        os.replace(part, state)
    except Exception as e:
        record["error"] = f"{type(e).__name__}: {e}"
    return record


def resumed(job: dict) -> dict | None:
    """The state record of a file already stamped and verified from the same source object."""
    try:
        record = json.loads(Path(job["state"]).read_text())
    except (OSError, json.JSONDecodeError):
        return None
    stamped = job["stamped"]
    intact = (record.get("source") == job["source"] and os.path.exists(stamped)
              and os.path.getsize(stamped) == record.get("size")
              and (not record.get("index") or os.path.exists(stamped + record["index"]))
              and (job["class"] != "per_phenotype" or not job["stats"] or "stats" in record))
    return record if intact else None


def regenerate_stats(records: dict, classes: dict, dirs: dict, suffix: str) -> dict:
    """Write the stats objects from the stamped records into `stamped/` and from the original
    records into `proof/`, and report whether the proof reproduces the originals."""
    sys.path.insert(0, str(HERE))
    from credible_set_stats import CredibleSetStats, get_tsv_header, stats_to_tsv_row, write_stats_json

    pairs = []
    stems = defaultdict(list)
    for rel in sorted(r for r, c in classes.items() if c == "per_phenotype"):
        stats = records[rel]["stats"]
        if len(stats["original"]) != len(stats["stamped"]):
            die(f"{rel}: stamping changed its number of traits")
        for before, after in zip(stats["original"], stats["stamped"]):
            pair = {"rel": rel, "original": CredibleSetStats(**before), "stamped": CredibleSetStats(**after),
                    "used": False}
            pairs.append(pair)
            stems[os.path.basename(rel).removesuffix(suffix)].append(pair)
    out = {"records": len(pairs), "problems": [], "changed": {}, "reproduced": {}}

    for rel in sorted(r for r, c in classes.items() if c == "stats_aggregate"):
        original = (dirs["original"] / rel).read_text()
        by_row = defaultdict(deque)
        for pair in pairs:
            by_row[stats_to_tsv_row(pair["original"])].append(pair)
        proof, stamped = [get_tsv_header()], [get_tsv_header()]
        for n, row in enumerate(original.splitlines()[1:], 2):
            if not by_row[row]:
                out["problems"].append(f"{rel} line {n} is not reproduced from the per-phenotype files: {row[:200]}")
                continue
            pair = by_row[row].popleft()
            pair["used"] = True
            proof.append(row)
            stamped.append(stats_to_tsv_row(pair["stamped"]))
        for pair in pairs:
            # the producers that drop a null trait from the stats still write its rows
            if not pair["used"] and pair["original"].trait is not None:
                out["problems"].append(f"{pair['rel']}: trait {pair['original'].trait} is missing from {rel}")
        for base, lines in (("proof", proof), ("stamped", stamped)):
            path = dirs[base] / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("".join(line + "\n" for line in lines))
        out["reproduced"][rel] = (dirs["proof"] / rel).read_text() == original
        out["changed"][rel] = sum(a != b for a, b in zip(proof, stamped))

    for rel in sorted(r for r, c in classes.items() if c == "stats_json"):
        original = (dirs["original"] / rel).read_bytes()
        proof, stamped = dirs["proof"] / rel, dirs["stamped"] / rel
        for path in (proof, stamped):
            path.parent.mkdir(parents=True, exist_ok=True)
        # a stats.json is named after its per-phenotype file: that file's name up to the
        # suffix, then a dot. Names hold dots themselves, so every dot is a candidate end
        name = os.path.basename(rel)
        match = None
        for pair in (p for i, c in enumerate(name) if c == "." for p in stems.get(name[:i], [])):
            write_stats_json(pair["original"], str(proof))
            if proof.read_bytes() == original:
                match = pair
                break
        out["reproduced"][rel] = match is not None
        if match is None:
            out["problems"].append(f"{rel} is not reproduced from any per-phenotype file")
            continue
        write_stats_json(match["stamped"], str(stamped))
        out["changed"][rel] = int(stamped.read_bytes() != original)
    return out


def report_rows(classes, listing, records, stats, out_sizes, delivered) -> list[dict]:
    rows = []
    for rel in sorted(classes):
        cls = classes[rel]
        row = {"object": rel, "class": cls, "action": action(cls), "bytes_in": listing[rel]["size"],
               "bytes_out": out_sizes.get(rel, "")}
        ok = cls in SKIPPED or rel in out_sizes
        record = records.get(rel)
        if cls in STAMP_MODE:
            ok = ok and record is not None and "error" not in record
            if ok:
                v = record["verify"]
                row.update(rows_in=v["rows_original"], rows_out=v["rows_stamped"],
                           header_identical=v["header_identical"],
                           rows_differing_outside_annotation=v["rows_differing_outside_annotation"],
                           matched=record["stamp"]["matched"])
                for col in ANNOTATION:
                    row.update({f"{col}_{k}": n for k, n in v[col].items()})
            elif record:
                row["error"] = record["error"].replace("\n", " ")[:300]
        elif cls in ("stats_aggregate", "stats_json"):
            row["original_reproduced"] = stats["reproduced"].get(rel, False)
            row["stats_rows_changed"] = stats["changed"].get(rel, "")
            ok = ok and row["original_reproduced"]
        if delivered is not None and cls not in SKIPPED:
            row["destination_verified"] = rel in delivered
            ok = ok and row["destination_verified"]
        row["ok"] = ok
        rows.append(row)
    return rows


COLUMNS = ["object", "class", "action", "bytes_in", "bytes_out", "rows_in", "rows_out", "header_identical",
           "rows_differing_outside_annotation", "matched",
           *(f"{c}_{k}" for c in ANNOTATION for k in ("na_before", "na_after", "changed_non_na")),
           "original_reproduced", "stats_rows_changed", "destination_verified", "ok", "error"]


def summarise(rows: list[dict], header: dict, checks: dict) -> str:
    lines = [f"{k}: {v}" for k, v in header.items()]
    lines += ["", f"{'class':<16}{'action':<14}{'in':>7}{'out':>7}{'rows in':>13}{'rows out':>13}"
              + "".join(f"{c[:11] + ' NA%':>24}" for c in ANNOTATION) + f"{'changed non-NA':>20}"]
    total = Counter()
    for cls in sorted({r["class"] for r in rows}):
        of = [r for r in rows if r["class"] == cls]
        n = Counter()
        for r in of:
            n.update({k: v for k, v in r.items() if isinstance(v, int) and not isinstance(v, bool)})
        n["in"], n["out"] = len(of), sum(r["bytes_out"] != "" for r in of)
        line = f"{cls:<16}{action(cls):<14}{n['in']:>7}{n['out']:>7}"
        if cls in STAMP_MODE:
            total.update(n)
            line += stamped_columns(n)
        lines.append(line)
    lines.append(f"{'all stamped':<30}{total['in']:>7}{total['out']:>7}{stamped_columns(total)}")
    lines += ["", *(f"{'PASS' if ok else 'FAIL'}  {name}" for name, ok in checks.items()),
              "", "RESULT: " + ("PASS" if all(checks.values()) else "FAIL")]
    return "\n".join(lines) + "\n"


def stamped_columns(n: Counter) -> str:
    def share(k):
        return f"{100 * n[k] / n['rows_in']:.2f}" if n["rows_in"] else "-"

    return (f"{n['rows_in']:>13}{n['rows_out']:>13}"
            + "".join(f"{share(c + '_na_before') + ' -> ' + share(c + '_na_after'):>24}" for c in ANNOTATION)
            + f"{' + '.join(str(n[c + '_changed_non_na']) for c in ANNOTATION):>20}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--source", required=True, help="prefix of the served resource: gs://bucket/prefix/ or a directory")
    p.add_argument("--dest", required=True, help="new prefix to write the stamped resource under")
    p.add_argument("--combined", required=True, help="name of the combined variant-sorted file (the profile's all_cs_file)")
    p.add_argument("--qtl", help="name of the gene-indexed QTL file (the profile's all_cs_qtl_file)")
    p.add_argument("--per-phenotype-prefix", default="",
                   help="the profile's prefix, relative to --source, e.g. individual/")
    p.add_argument("--per-phenotype-suffix", help="the profile's suffix_95, e.g. .SUSIE.munged.tsv")
    p.add_argument("--stats-file", default="credible_set_stats.tsv", help="name of the aggregate stats file")
    p.add_argument("--copy", action="append", default=[], metavar="GLOB", help="copy matching objects as they are")
    p.add_argument("--skip", action="append", default=[], metavar="GLOB", help="leave matching objects out")
    p.add_argument("--consequence", help="bgzip + tabix consequence file")
    p.add_argument("--consequence-version", help="goes into the object metadata, e.g. 4.1.1-vep115")
    p.add_argument("--consequence-source", default="gnomAD")
    p.add_argument("--staging", help="local working directory, on a disk with room for the resource twice over")
    p.add_argument("--workers", type=int, default=3, help="files stamped at a time")
    p.add_argument("--max-rss-gb", type=float, help="passed to annotate_consequence.py")
    p.add_argument("--max-lookup-rows", type=int, help="passed to annotate_consequence.py")
    p.add_argument("--dry-run", action="store_true", help="print the classification and the planned names, touch nothing")
    args = p.parse_args()

    src, dst = store(args.source), store(args.dest)
    if src.url.startswith(dst.url) or dst.url.startswith(src.url):
        die(f"destination {dst.url} and source {src.url} must not contain one another")
    t0 = time.time()
    timing = {}

    def phase(name: str, since: float) -> float:
        timing[name] = round(time.time() - since, 1)
        log(f"{name}: {timing[name]} s")
        return time.time()

    listing = src.list()
    if not listing:
        die(f"nothing under {src.url}")
    classes = classify(listing, args)
    planned = {rel for rel, cls in classes.items() if cls not in SKIPPED}
    if args.dry_run:
        print("class\taction\tsource\tdestination")
        for rel in sorted(classes):
            print(f"{classes[rel]}\t{action(classes[rel])}\t{src.url}{rel}\t{dst.url + rel if rel in planned else '-'}")
        for cls, n in sorted(Counter(classes.values()).items()):
            print(f"# {cls}: {n}", file=sys.stderr)
        return

    for flag in ("consequence", "consequence_version", "staging"):
        if not getattr(args, flag):
            p.error(f"--{flag.replace('_', '-')} is required without --dry-run")
    if not (os.path.exists(args.consequence) and any(os.path.exists(args.consequence + e) for e in INDEX_EXTS)):
        die(f"{args.consequence}: needs the file and its index")
    existing = dst.list()
    foreign = sorted(set(existing) - planned - {REPORT_TSV, REPORT_SUMMARY})
    if foreign:
        die(f"{dst.url} holds {len(foreign)} objects this run would not write, e.g. {foreign[:3]}")
    t = phase("list", t0)

    staging = Path(args.staging).resolve()
    dirs = {name: staging / name for name in ("original", "stamped", "copied", "proof", "state", "report", "tmp")}
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(dirs["tmp"])
    # one thread per worker process, or the stats readers together outnumber the cores
    os.environ["POLARS_MAX_THREADS"] = "1"
    metadata = {
        "consequence-source": args.consequence_source,
        "consequence-version": args.consequence_version,
        "consequence-file": os.path.basename(args.consequence),
        f"consequence-{dst.digest_name}": dst.metadata_digest(dst.file_digests([args.consequence])[args.consequence]),
        "stamped-from": src.url,
    }
    if any(c in str(v) for v in metadata.values() for c in ",="):
        die("metadata values cannot hold a comma or an equals sign")
    binding = {"metadata": metadata, "dest": dst.url,
               "rules": {k: getattr(args, k) for k in ("combined", "qtl", "per_phenotype_prefix",
                                                       "per_phenotype_suffix", "stats_file", "copy", "skip")}}
    bound = staging / "run.json"
    if bound.exists() and json.loads(bound.read_text()) != binding:
        die(f"{staging} belongs to a run with other arguments (see run.json); use a fresh --staging")
    bound.write_text(json.dumps(binding, indent=2))

    src.fetch(dirs["original"])
    short = [r for r in planned if not (dirs["original"] / r).is_file()
             or (dirs["original"] / r).stat().st_size != listing[r]["size"]]
    if short:
        die(f"{len(short)} objects did not arrive whole in {dirs['original']}, e.g. {sorted(short)[:3]}")
    t = phase("download", t)

    stamper_args = [x for flag in ("max_rss_gb", "max_lookup_rows") if getattr(args, flag) is not None
                    for x in (f"--{flag.replace('_', '-')}", str(getattr(args, flag)))]
    want_stats = any(c in ("stats_aggregate", "stats_json") for c in classes.values())
    jobs = [{"rel": rel, "class": cls, "source": listing[rel]["digest"], "consequence": args.consequence,
             "original": str(dirs["original"] / rel), "stamped": str(dirs["stamped"] / rel),
             "state": str(dirs["state"] / f"{rel}.json"), "stamper_args": stamper_args, "stats": want_stats,
             "threads": args.workers if cls in SERIAL else 1}
            for rel, cls in sorted(classes.items()) if cls in STAMP_MODE]
    records = {}
    todo = []
    for job in jobs:
        record = resumed(job)
        if record:
            records[job["rel"]] = record
        else:
            todo.append(job)
    log(f"{len(jobs)} files to stamp, {len(records)} already done")
    for job in (j for j in todo if j["class"] in SERIAL):
        records[job["rel"]] = stamp_one(job)
    t = phase("stamp combined and qtl", t)
    with ProcessPoolExecutor(max_workers=args.workers, initializer=sys.path.insert, initargs=(0, str(HERE))) as pool:
        futures = [pool.submit(stamp_one, j) for j in todo if j["class"] not in SERIAL]
        for n, future in enumerate(as_completed(futures), 1):
            record = future.result()
            records[record["rel"]] = record
            if not n % 500:
                log(f"{n}/{len(futures)} per-phenotype files")
    t = phase("stamp per phenotype", t)
    failed = sorted(r for r, rec in records.items() if "error" in rec)
    for rel in failed[:20]:
        log(f"FAILED {rel}: {records[rel]['error']}")

    for rel in sorted(r for r, c in classes.items() if c == "copy"):
        with open(dirs["original"] / rel, "rb") as f:
            header = f.readline(1 << 16).rstrip(b"\n").split(b"\t")
        if all(c.encode() in header for c in ANNOTATION):
            die(f"{rel}: carries the annotation columns, so a byte-identical copy would keep the old annotation")
        target = dirs["copied"] / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            shutil.copyfile(dirs["original"] / rel, target)

    stats = {"records": 0, "problems": [], "changed": {}, "reproduced": {}}
    if want_stats and not failed:
        stats = regenerate_stats(records, classes, dirs, args.per_phenotype_suffix)
        for problem in stats["problems"][:20]:
            log(f"STATS {problem}")
    t = phase("stats", t)

    def staged(rel: str) -> Path:
        return dirs["copied" if classes[rel] == "copy" else "stamped"] / rel

    out_sizes = {rel: staged(rel).stat().st_size for rel in planned if staged(rel).is_file()}
    counts_in = Counter(classes[r] for r in planned)
    checks = {
        "every file stamped and verified by annotate_consequence.py --verify": not failed,
        "original stats reproduced from the original rows": not stats["problems"] and all(stats["reproduced"].values()),
        "every planned object produced locally": set(out_sizes) == planned,
    }
    delivered = None
    if all(checks.values()):
        digests = dst.file_digests([staged(r) for r in sorted(planned)])
        clash = sorted(r for r in planned if r in existing and existing[r]["digest"] != digests[str(staged(r))])
        if clash:
            die(f"{dst.url} already holds {len(clash)} objects with other content, e.g. {clash[:3]}")
        t = phase("checksum", t)
        dst.push(dirs["stamped"], metadata)
        if any(dirs["copied"].iterdir()):
            dst.push(dirs["copied"], {})
        t = phase("upload", t)
        after = dst.list()
        delivered = {
            r for r in planned
            if r in after and after[r]["digest"] == digests[str(staged(r))]
            and (not dst.holds_metadata or classes[r] == "copy"
                 or all(after[r]["metadata"].get(k) == v for k, v in metadata.items()))
        }
        counts_out = Counter(classes[r] for r in delivered)
        checks["objects in equal objects out for every class"] = counts_in == counts_out
        checks["every destination object has the staged checksum" + (" and metadata" if dst.holds_metadata else "")] = (
            delivered == planned)
        checks["nothing else at the destination"] = not set(after) - planned - {REPORT_TSV, REPORT_SUMMARY}
        t = phase("verify destination", t)

    rows = report_rows(classes, listing, records, stats, out_sizes, delivered)
    checks["every object's own checks"] = all(r["ok"] for r in rows)
    header = {"source": src.url, "destination": dst.url, **metadata,
              "stats records": stats["records"],
              "stats objects changed by the stamp": sum(1 for n in stats["changed"].values() if n)}
    summary = summarise(rows, header, checks)
    (dirs["report"] / REPORT_TSV).write_text(
        "\t".join(COLUMNS) + "\n" + "".join("\t".join(str(r.get(c, "")) for c in COLUMNS) + "\n" for r in rows))
    (dirs["report"] / REPORT_SUMMARY).write_text(summary)
    ok = all(checks.values())
    if ok:
        # only a passing report is published: one beside the objects says they were proven
        dst.push(dirs["report"], metadata)
        published = dst.list()
        local = dst.file_digests([dirs["report"] / n for n in (REPORT_TSV, REPORT_SUMMARY)])
        ok = all(published.get(n, {}).get("digest") == local[str(dirs["report"] / n)]
                 for n in (REPORT_TSV, REPORT_SUMMARY))
        if not ok:
            log(f"the report at {dst.url} is not the one this run wrote")

    wall = time.time() - t0
    timing.update(wall_seconds=round(wall, 1), objects=len(planned), files_stamped=len(todo),
                  objects_per_second=round(len(planned) / wall, 2),
                  peak_child_rss_mb=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss >> 10,
                  peak_stamper_rss_mb=max((r["stamp"]["peak_rss_mb"] for r in records.values() if "stamp" in r), default=0))
    (staging / "timing.json").write_text(json.dumps(timing, indent=2))
    print(summary)
    print(json.dumps(timing))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
