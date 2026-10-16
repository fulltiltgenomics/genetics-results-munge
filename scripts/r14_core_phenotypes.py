#!/usr/bin/env python3
"""Build the FinnGen R14 core GWAS phenotype metadata JSON from the release's pheno_n TSV.

    python3 scripts/r14_core_phenotypes.py --output finngen_r14_pheno.json --stage

The pheweb `pheno-list.json` previously used here lists 97 fewer endpoints than the release
has sumstats for, with no exclusion recorded anywhere, so the release's own
`finngen_R14_pheno_n.tsv` is the source instead: one row per released phenocode, with the
columns every consumer reads (`phenocode`, `name`, `category`, `category_index`,
`num_cases`, `num_controls`). `name` is written as `phenostring`, the key results-api,
results-db and munge_hla.sh read first. A quantitative trait has `num_controls` 0 in the TSV
and keeps it: results-db derives trait_type from a zero control count.

`--previous` names the previously served JSON (local path or gs://). Every phenocode shared
with it must agree on every TSV-derived field, or nothing is written: a regeneration can only
add. The pheweb-computed keys the TSV cannot supply (`gc_lambda`, `num_gw_significant` and the
`*_prev` release-over-release counts) are carried over from it per phenocode, since nothing
else holds them; a phenocode pheweb never saw has none. `assoc_files`, a path inside the
pheweb run's own disk, is not carried.
"""

import argparse
import csv
import io
import json
import subprocess
import sys
from pathlib import Path

from sumstat_utils import upload_to_gcs

PHENO_N_TSV = "gs://finngen-production-library-green/finngen_R14/finngen_R14_analysis_data/finngen_R14_pheno_n.tsv"
STAGE_PATH = "gs://finngen-commons/results_api_data/mapping_files/finngen_r14_pheno.json"
FIELDS = ("phenostring", "category", "category_index", "num_cases", "num_controls")
CARRIED = ("num_gw_significant", "num_cases_prev", "num_controls_prev", "num_gw_significant_prev", "gc_lambda")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pheno-n", default=PHENO_N_TSV, help="release pheno_n TSV, local path or gs://")
    parser.add_argument("--output", required=True, help="path to write the phenotype JSON to")
    parser.add_argument("--previous", metavar="PREVIOUS_JSON",
                        help="previously served JSON: must agree on every shared phenocode, "
                             "and supplies the pheweb-computed keys the TSV lacks")
    parser.add_argument("--stage", action="store_true", help=f"upload the JSON to {STAGE_PATH}")
    return parser.parse_args()


def read_text(path: str) -> str:
    if path.startswith("gs://"):
        return subprocess.run(["gcloud", "storage", "cat", path], check=True, capture_output=True, text=True).stdout
    return Path(path).read_text()


def parse_row(row: dict) -> dict:
    return {
        "phenocode": row["phenocode"],
        "phenostring": row["name"],
        "category": row["category"],
        "category_index": int(row["category_index"]),
        "num_cases": int(row["num_cases"]),
        "num_controls": int(row["num_controls"]),
    }


def build(tsv_text: str) -> list[dict]:
    items = [parse_row(row) for row in csv.DictReader(io.StringIO(tsv_text), delimiter="\t")]
    codes = [item["phenocode"] for item in items]
    dupes = sorted({c for c in codes if codes.count(c) > 1})
    if dupes:
        sys.exit(f"duplicate phenocodes in pheno_n: {dupes}")
    return items


def compare(items: list[dict], previous: list[dict]) -> list[str]:
    """Differences on shared phenocodes, plus the codes each side lacks; empty means additive."""
    new = {item["phenocode"]: item for item in items}
    problems = []
    seen = set()
    for old in previous:
        code = old["phenocode"]
        if code in seen:
            continue
        seen.add(code)
        if code not in new:
            problems.append(f"{code}: in previous JSON, not in pheno_n")
            continue
        for field in FIELDS:
            if str(old.get(field)) != str(new[code][field]):
                problems.append(f"{code}.{field}: previous {old.get(field)!r}, pheno_n {new[code][field]!r}")
    return problems


def carry_over(items: list[dict], previous: list[dict]) -> int:
    """Copy the pheweb-computed keys onto each item pheweb had; returns how many got them."""
    old = {p["phenocode"]: p for p in previous}
    carried = 0
    for item in items:
        source = old.get(item["phenocode"])
        if source is None:
            continue
        item.update({k: source[k] for k in CARRIED if k in source})
        carried += 1
    return carried


def main() -> int:
    args = parse_args()
    items = build(read_text(args.pheno_n))
    if args.previous:
        previous = json.loads(read_text(args.previous))
        problems = compare(items, previous)
        for line in problems:
            print(f"  {line}", file=sys.stderr)
        if problems:
            sys.exit(f"{len(problems)} differences against {args.previous}; nothing written")
        carried = carry_over(items, previous)
        print(f"  {args.previous}: every shared phenocode agrees; {carried} carry pheweb keys, "
              f"{len(items) - carried} added without them")
    Path(args.output).write_text(json.dumps(items, indent=2) + "\n")
    print(f"  wrote {len(items)} phenotypes to {args.output}")
    if args.stage:
        upload_to_gcs(args.output, STAGE_PATH)
    return 0


if __name__ == "__main__":
    sys.exit(main())
