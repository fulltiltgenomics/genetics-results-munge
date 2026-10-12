#!/usr/bin/env python3
"""Build the phenotype metadata JSON for the FinnGen R14 + All of Us + UKBB + MVP
multi-cohort meta-analysis from the release's phenotype mapping TSV.

    python3 scripts/r14_meta_phenotypes.py --output finngen_r14_mvp_ukbb_aou_meta.json --stage

Input: `FinnGen_R14_AoU_UKBB_MVP_phenotype_mapping.tsv` from the green release prefix, one
row per meta phenocode. WIDE endpoints and Labs rows are ordinary phenocodes here. Every
cohort contributes a `<cohort>_n_cases` / `<cohort>_n_controls` pair, `NA` where the cohort
did not take part. Binary rows carry both or neither; quantitative rows (Labs, Quantitative
endpoints) carry `num_controls` NA and only `<cohort>_n_cases`, which is then the cohort's
sample size.

Output is the pheweb-shaped list: `phenocode`, `name`, `category`, `category_index`,
`num_cases`, `num_controls`, `contributing_cohorts`, `cohort_families`, `n_by_cohort`.
`num_controls` (and per-cohort `n_controls`) is omitted rather than written as NA, so a
consumer can tell a quantitative trait by the missing key. The per-phenotype field is
`contributing_cohorts` and not `cohorts` because the results registry already uses
`cohorts` at dataset level for collection tokens (finngen/ukb/mvp/aou).

Cohort presence comes from the sumstats, not the TSV: the header of every
`<phenocode>_meta_out.tsv.gz` is read and a cohort contributes only if its `<cohort>_beta`
column is there. The TSV supplies names, categories and n. A TSV cohort with n but no header
columns is dropped with a WARNING, and `num_cases`/`num_controls` are recomputed from the
cohorts kept (the TSV totals are the per-cohort sums, see docs/r14-meta-serving.md). A header
cohort with no TSV n is an ERROR and nothing is written, since its n_by_cohort would be
incomplete. `--skip-header-scan` gives the TSV-only view for offline use; `--stage` refuses it.
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from sumstat_utils import upload_to_gcs

RELEASE_PREFIX = "gs://finngen-production-library-green/finngen_R14/finngen_R14_analysis_data/meta_analysis"
MAPPING_TSV = f"{RELEASE_PREFIX}/FinnGen_R14_AoU_UKBB_MVP_phenotype_mapping.tsv"
SUMSTATS = RELEASE_PREFIX + "/summary_stats/{phenocode}_meta_out.tsv.gz"
STAGE_PATH = "gs://finngen-commons/results_api_data/mapping_files/finngen_r14_mvp_ukbb_aou_meta.json"

COHORTS = [
    "fg",
    "aou_eur", "aou_afr", "aou_amr", "aou_eas", "aou_sas", "aou_mid",
    "ukbb_eur",
    "mvp_eur", "mvp_afr", "mvp_amr",
]
FAMILIES = ["fg", "aou", "mvp", "ukbb"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mapping", default=MAPPING_TSV, help="phenotype mapping TSV, local path or gs://")
    parser.add_argument("--output", required=True, help="path to write the phenotype JSON to")
    parser.add_argument("--skip-header-scan", action="store_true",
                        help="take cohort presence from the TSV alone (offline use; refused with --stage)")
    parser.add_argument("--workers", type=int, default=16, help="parallel sumstats header reads")
    parser.add_argument("--stage", action="store_true", help=f"upload the JSON to {STAGE_PATH}")
    return parser.parse_args()


def _count(value: str) -> int | None:
    return None if value == "NA" else int(value)


def parse_row(row: dict) -> dict:
    item = {
        "phenocode": row["phenocode"],
        "name": row["name"],
        "category": row["category"],
        "category_index": int(row["category_index"]),
        "num_cases": int(row["num_cases"]),
    }
    num_controls = _count(row["num_controls"])
    if num_controls is not None:
        item["num_controls"] = num_controls

    n_by_cohort = {}
    for cohort in COHORTS:
        n_cases = _count(row[f"{cohort}_n_cases"])
        if n_cases is None:
            continue
        counts = {"n_cases": n_cases}
        n_controls = _count(row[f"{cohort}_n_controls"])
        if n_controls is not None:
            counts["n_controls"] = n_controls
        n_by_cohort[cohort] = counts

    _set_cohorts(item, n_by_cohort)
    return item


def _set_cohorts(item: dict, n_by_cohort: dict) -> None:
    item["contributing_cohorts"] = list(n_by_cohort)
    item["cohort_families"] = [
        f for f in FAMILIES if any(c == f or c.startswith(f + "_") for c in n_by_cohort)
    ]
    item["n_by_cohort"] = n_by_cohort


def restrict_to_header(item: dict, found: set[str]) -> tuple[list[str], list[str]]:
    """Keep only the cohorts whose columns the sumstats header has. Returns (warnings, errors)."""
    tsv = item["n_by_cohort"]
    errors = [f"ERROR {item['phenocode']}: {c} in sumstats header but has no n in the mapping TSV"
              for c in sorted(found - set(tsv))]
    dropped = [c for c in tsv if c not in found]
    warnings = [f"WARNING {item['phenocode']}: dropped {c}, in the mapping TSV with {tsv[c]} "
                f"but not in the sumstats header" for c in dropped]
    if dropped:
        kept = {c: n for c, n in tsv.items() if c in found}
        _set_cohorts(item, kept)
        item["num_cases"] = sum(n["n_cases"] for n in kept.values())
        if "num_controls" in item:
            # KeyError if a kept cohort of a binary row has n_controls NA (none does today)
            item["num_controls"] = sum(n["n_controls"] for n in kept.values())
    return warnings, errors


def read_mapping(path: str) -> list[dict]:
    if path.startswith("gs://"):
        text = subprocess.run(["gcloud", "storage", "cat", path], check=True, capture_output=True, text=True).stdout
    else:
        text = Path(path).read_text()
    rows = csv.DictReader(text.splitlines(), delimiter="\t")
    return [parse_row(row) for row in rows]


def header_cohorts(header: list[str]) -> set[str]:
    """Every `<x>_beta` column that is a single cohort's, known or not, so a cohort the
    TSV never mentions still shows up as a disagreement."""
    return {
        col[: -len("_beta")]
        for col in header
        if col.endswith("_beta") and not col.startswith(("leave_cohort_", "all_inv_var_meta"))
    }


def read_header(phenocode: str, workdir: str) -> list[str]:
    # tabix -H fetches only the first block; reading gs:// needs GCS_OAUTH_TOKEN.
    # tabix also downloads each file's ~2.4 MB .tbi into its cwd, hence a throwaway workdir
    line = subprocess.run(
        ["tabix", "-H", SUMSTATS.format(phenocode=phenocode)],
        check=True, capture_output=True, text=True, cwd=workdir,
    ).stdout.splitlines()[0]
    return line.lstrip("#").split("\t")


def apply_headers(items: list[dict], workers: int) -> int:
    """Restrict every item to its sumstats header's cohorts; returns the number of errors."""
    if "GCS_OAUTH_TOKEN" not in os.environ:
        os.environ["GCS_OAUTH_TOKEN"] = subprocess.run(
            ["gcloud", "auth", "print-access-token"], check=True, capture_output=True, text=True
        ).stdout.strip()

    def apply(item: dict) -> tuple[list[str], list[str]]:
        try:
            found = header_cohorts(read_header(item["phenocode"], workdir))
        except (subprocess.CalledProcessError, IndexError) as e:
            return [], [f"ERROR {item['phenocode']}: header unreadable ({e})"]
        return restrict_to_header(item, found)

    with tempfile.TemporaryDirectory() as workdir, ThreadPoolExecutor(workers) as pool:
        results = list(pool.map(apply, items))
    warnings = [w for ws, _ in results for w in ws]
    errors = [e for _, es in results for e in es]
    for line in warnings + errors:
        print(line, file=sys.stderr)
    changed = sum(bool(ws) for ws, _ in results)
    print(f"header scan: {changed} of {len(items)} phenotypes lost a TSV cohort "
          f"({len(warnings)} cohorts dropped), {len(errors)} errors", file=sys.stderr)
    return len(errors)


def main() -> int:
    args = parse_args()
    if args.stage and args.skip_header_scan:
        sys.exit("--stage needs the header scan: the staged JSON must carry header-derived cohorts")
    items = read_mapping(args.mapping)
    print(f"{len(items)} phenotypes, {sum('num_controls' not in i for i in items)} without num_controls")

    if not args.skip_header_scan and apply_headers(items, args.workers):
        return 1

    Path(args.output).write_text(json.dumps(items, indent=2) + "\n")
    print(f"wrote {args.output}")
    if args.stage:
        upload_to_gcs(args.output, STAGE_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
