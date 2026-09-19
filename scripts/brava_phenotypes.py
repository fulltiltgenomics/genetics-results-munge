#!/usr/bin/env python3
"""Build the phenotype metadata JSON for the BRaVa exome-wide rare variant meta-analysis
(medRxiv 2026.05.21.26353759), one entry per (phenotype, meta-analysis stratum) that has a
gene burden result file.

    python3 scripts/brava_phenotypes.py --output brava_pheno.json
    python3 scripts/brava_phenotypes.py --output brava_pheno.json --check --stage

Output is the pheweb-shaped list the results platform's metadata harmonizers read:
`phenocode`, `phenostring`, `category`, and either `num_cases`/`num_controls` (binary) or
`num_samples` (quantitative). Nothing consumes it yet -- the harmonizer that will is not
written -- so the key names follow the `finngen_r13` and `quantitative_pheweb` harmonizers
in genetics-results-api's `app/services/metadata_harmonizer.py`.

Inputs:
  the preprint's supplementary workbook (`--supplement`, an https URL or a local path)
    Table S1  32 binary phenotypes: Description, Sex, Phenotype ID
    Table S2  11 quantitative phenotypes, same three columns
    Table S4  binary N cases / N controls per (phenotype, ancestry, biobank)
    Table S5  quantitative N per (phenotype, ancestry, biobank)
    Table S6  binary totals per phenotype, 33 rows
    Table S7  quantitative totals per phenotype
  the gene result file listing (`--gene-prefix`, or a saved listing via `--gene-listing`)
    `<code>_gene_meta_analysis_100_cutoff.tsv.gz` is the cross-ancestry meta-analysis and
    `<code>_gene_meta_analysis_100_cutoff.<STRATUM>.tsv.gz` a stratum. Which strata exist
    varies by phenotype -- a phenotype with too few cases in an ancestry has no file for it
    -- so the listing, not a fixed list, is what decides the output.

Facts established from the workbook and from the results themselves rather than taken from
the paper, and re-derived on every run so a revised workbook cannot pass silently:

  - **Table S5's repeated rows are duplicates and are dropped.** 45 of its 207 distinct
    (phenotype, ancestry, biobank) keys appear on more than one row. Summing every row
    reproduces each published Table S7 total exactly -- which is how those totals were
    computed -- but it double-counts: the deduplicated Height EUR sum is 710,271 and the
    largest per-variant NS in `Height_ALL_variant_meta_analysis_100_cutoff.EUR.vcf.gz` is
    710,270, while the raw sum is 920,670. So the deduplicated sums are used, S7 is not, and
    every quantitative `num_samples` here is 0-30% below the published total. `--check`
    prints both. One key, BMI / AMR / all-of-us, carries two different N (57,022 and
    38,834); the smaller is taken -- see the dedupe loop in `per_ancestry` for why.
  - Table S4 has no repeated keys and its sums reproduce Table S6 exactly, for all 33 binary
    phenotypes including HipRep. That reconciliation is a gate in `check_totals`: the binary
    counts are derived the same way as the quantitative ones, and S6 is what proves the
    derivation right.
  - HipRep (hip replacement, a procedure endpoint) is in Table S6 but in neither S1 nor S2,
    so its description comes from S6 and its category defaults to Both.
  - `non_EUR` sums *every* non-EUR ancestry, MID included. MID contributes to the binary
    totals and has no stratum file of its own, so it is only ever reachable through
    `non_EUR`; the quantitative tables have no MID rows at all.

Sex tags in the file codes: six Female-only phenotypes carry `_F` and the rest carry `_ALL`.
The `_F` is KEPT in the phenocode and the `_ALL` is stripped, so a phenocode is either
`BreastCanc_F` or `AFib`, and a stratum appends `|<STRATUM>`.

`--check` compares a sample of the output against the per-variant NS/NC maxima of the
matching VCF. A maximum is a lower bound on the total, not the total -- no variant need be
typed in every sample -- so it prints the comparison and never fails on it. A maximum that
*exceeds* its total is the failure worth looking for, and that is what caught Table S7.
"""

import argparse
import gzip
import json
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from sumstat_utils import fetch, fetch_gs, upload_to_gcs

SUPPLEMENT_URL = (
    "https://www.medrxiv.org/content/medrxiv/early/2026/05/24/2026.05.21.26353759/DC2/"
    "embed/media-2.xlsx?download=true"
)
GENE_PREFIX = "gs://daly-genetics-results/raw/brava/gene/"
VARIANT_PREFIX = "gs://daly-genetics-results/raw/brava/variant/"
STAGE_PATH = "gs://daly-genetics-results/mapping_files/brava_pheno.json"

GENE_FILE = re.compile(
    r"^(?P<code>.+?)_gene_meta_analysis_100_cutoff(?:\.(?P<stratum>[A-Za-z_]+))?\.tsv\.gz$"
)
STRATUM_ORDER = ["AFR", "AMR", "EAS", "EUR", "SAS", "non_EUR"]

CHECK_SAMPLE = ["AFib", "AFib|EUR", "Height"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="path to write the phenotype JSON to")
    parser.add_argument("--supplement", default=SUPPLEMENT_URL, help="supplementary workbook URL or path")
    parser.add_argument("--gene-prefix", default=GENE_PREFIX, help="gs:// prefix holding the gene result files")
    parser.add_argument(
        "--gene-listing",
        help="file of gene result file names (one per line, bare or gs://) to use instead of listing --gene-prefix",
    )
    parser.add_argument("--variant-prefix", default=VARIANT_PREFIX, help="gs:// prefix holding the VCFs --check reads")
    parser.add_argument(
        "--billing-project",
        help="requester-pays project, needed only when a prefix is the source bucket gs://brava-meta-analysis/",
    )
    parser.add_argument("--cache-dir", default=str(Path.home() / "brava_munge" / "cache"))
    parser.add_argument(
        "--check",
        nargs="*",
        metavar="PHENOCODE",
        help=f"compare against the VCFs' per-variant NS/NC maxima; defaults to {' '.join(CHECK_SAMPLE)}",
    )
    parser.add_argument("--stage", action="store_true", help=f"upload the JSON to {STAGE_PATH}")
    return parser.parse_args()


def read_workbook(path: Path) -> dict[str, list[dict]]:
    """Read every sheet this script uses as a list of row dicts keyed by the header row."""
    import openpyxl

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    tables = {}
    for name in ("Table S1", "Table S2", "Table S4", "Table S5", "Table S6", "Table S7"):
        rows = workbook[name].iter_rows(values_only=True)
        header = [str(cell).strip() for cell in next(rows)]
        tables[name] = [
            dict(zip(header, row)) for row in rows if any(cell is not None for cell in row)
        ]
    workbook.close()
    return tables


def list_gene_files(prefix: str, billing_project: str | None) -> list[str]:
    command = ["gcloud", "storage", "ls", prefix.rstrip("/") + "/"]
    if billing_project:
        command.append(f"--billing-project={billing_project}")
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return result.stdout.splitlines()


def parse_listing(lines: list[str]) -> dict[str, set[str | None]]:
    """Map each file code to the strata it has files for; None is the cross-ancestry meta."""
    strata: dict[str, set[str | None]] = defaultdict(set)
    for line in lines:
        name = line.strip().rsplit("/", 1)[-1]
        if not name:
            continue
        match = GENE_FILE.match(name)
        if match is None:
            continue
        strata[match.group("code")].add(match.group("stratum"))
    if not strata:
        raise SystemExit("no gene result files matched; is the prefix empty or still copying?")
    return dict(strata)


SEX_TAG = re.compile(r"_(ALL|F)$")


def phenotype_id(code: str) -> str:
    """The workbook's Phenotype ID for a file code: the trailing sex tag is not part of it."""
    return SEX_TAG.sub("", code)


def phenocode(code: str, stratum: str | None) -> str:
    """`_ALL` is dropped and `_F` kept, so a phenocode carries the sex restriction where there is one."""
    base = code[: -len("_ALL")] if code.endswith("_ALL") else code
    return base if stratum is None else f"{base}|{stratum}"


def source_file_code(base: str) -> str:
    """The file code a phenocode came from: every file code carries a sex tag and `phenocode` above
    stripped the `_ALL` one, so reaching back to a source file puts it on again."""
    return base if base.endswith("_F") else base + "_ALL"


def per_ancestry(tables: dict[str, list[dict]]) -> tuple[dict, dict]:
    """Sum each table by (phenotype, ancestry): binary as (cases, controls), quantitative as N.

    Table S5 is deduplicated first -- see the module docstring for the measurement that
    settles that -- keeping the smaller N where a repeated key disagrees with itself.
    """
    binary: dict[tuple[str, str], list[float]] = defaultdict(lambda: [0.0, 0.0])
    for row in tables["Table S4"]:
        total = binary[(row["Phenotype ID"], row["Ancestry"])]
        total[0] += row["N cases"]
        total[1] += row["N controls"]

    deduplicated: dict[tuple[str, str, str], float] = {}
    for row in tables["Table S5"]:
        key = (row["Phenotype ID"], row["Ancestry"], row["Biobank ID"])
        n = row["N"]
        # BMI/AMR/all-of-us is the one key where the two rows disagree rather than repeat
        # (57,022 vs 38,834). 57,022 equals that biobank's whole AMR sample in Table S8, and
        # every other BMI/AMR biobank sits BELOW its S8 ancestry total (ccpm 9,681 < 9,766;
        # mgbb 1,185 < 3,612; pmbb 556 < 571; uk-biobank 502 < 508) -- so the value equal to
        # the S8 total is the anomalous one here, and 38,834 is the more plausible
        # BMI-measured subset. The VCFs can't settle it: the BMI AMR VCF has no NS-bearing
        # rows, and the meta's max NS (822,584) is below both candidates. Taking the smaller
        # is a judgement, not a derivation, and it affects only BMI|AMR (~18k) and the BMI
        # meta (~2%).
        deduplicated[key] = n if key not in deduplicated else min(deduplicated[key], n)
    quantitative: dict[tuple[str, str], float] = defaultdict(float)
    for (phenotype, ancestry, _), n in deduplicated.items():
        quantitative[(phenotype, ancestry)] += n

    return dict(binary), dict(quantitative)


def check_totals(tables: dict[str, list[dict]], binary: dict, quantitative: dict) -> None:
    """Reconcile the derived per-ancestry sums against the published totals.

    Binary is a gate: Table S4 sums to Table S6 exactly, and that is what proves the
    derivation right, so a disagreement stops the run. Quantitative is not, because the
    deduplication this script applies is a deliberate departure from Table S7; the size of
    that departure is printed instead.
    """
    problems = []
    for row in tables["Table S6"]:
        got = [0.0, 0.0]
        for (phenotype, _), total in binary.items():
            if phenotype == row["Phenotype ID"]:
                got[0] += total[0]
                got[1] += total[1]
        want = [row["N cases"], row["N controls"]]
        if got != want:
            problems.append(f"S4 sums {got} != S6 totals {want} for {row['Phenotype ID']}")
    if problems:
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        raise SystemExit(
            f"{len(problems)} binary phenotypes no longer reconcile with Table S6; the workbook's "
            "shape has changed and the summing rule has to be re-derived"
        )
    print(f"  Table S4 sums reconcile with all {len(tables['Table S6'])} Table S6 totals")

    for row in tables["Table S7"]:
        got = sum(n for (phenotype, _), n in quantitative.items() if phenotype == row["Phenotype ID"])
        print(f"    {row['Phenotype ID']:<8} deduplicated {int(got):>10,}   Table S7 {int(row['N']):>10,}")


def build_items(
    tables: dict[str, list[dict]], strata: dict[str, set[str | None]], binary: dict, quantitative: dict
) -> list[dict]:
    """One item per (phenotype, stratum) with a gene result file.

    Every count is summed from the per-biobank tables, the cross-ancestry meta included, so
    a phenotype's strata always add up to its meta.
    """
    is_quantitative = {row["Phenotype ID"] for row in tables["Table S2"]}

    description = {}
    category = {}
    for sheet in ("Table S6", "Table S7", "Table S1", "Table S2"):
        for row in tables[sheet]:
            description[row["Phenotype ID"]] = row["Description"]
            if "Sex" in row:
                category[row["Phenotype ID"]] = row["Sex"]

    def ancestries(phenotype: str, stratum: str | None) -> list[str]:
        table = quantitative if phenotype in is_quantitative else binary
        present = {ancestry for (name, ancestry) in table if name == phenotype}
        if stratum is None:
            return sorted(present)
        return sorted(present - {"EUR"}) if stratum == "non_EUR" else [stratum]

    items = []
    missing = []
    for code in sorted(strata):
        phenotype = phenotype_id(code)
        if phenotype not in description:
            missing.append(code)
            continue
        for stratum in sorted(strata[code], key=lambda s: (s is not None, STRATUM_ORDER.index(s) if s else -1)):
            item = {
                "phenocode": phenocode(code, stratum),
                "phenostring": description[phenotype] + (f" ({stratum})" if stratum else ""),
                "category": category.get(phenotype, "Both"),
            }
            parts = ancestries(phenotype, stratum)
            if phenotype in is_quantitative:
                item["num_samples"] = int(sum(quantitative[(phenotype, a)] for a in parts))
            else:
                item["num_cases"] = int(sum(binary[(phenotype, a)][0] for a in parts))
                item["num_controls"] = int(sum(binary[(phenotype, a)][1] for a in parts))
            items.append(item)

    if missing:
        raise SystemExit(f"gene result files for phenotypes the workbook does not describe: {missing}")
    return items


def vcf_maxima(path: Path) -> tuple[int | None, int | None]:
    """Largest per-variant NS and NC in a BRaVa meta-analysis VCF; NC is absent on quantitative traits."""
    max_ns = max_nc = None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            keys = fields[8].split(":")
            values = dict(zip(keys, fields[9].split(":")))
            for key, current in (("NS", max_ns), ("NC", max_nc)):
                value = values.get(key, ".")
                if value == ".":
                    continue
                best = int(float(value))
                if current is None or best > current:
                    if key == "NS":
                        max_ns = best
                    else:
                        max_nc = best
    return max_ns, max_nc


def cross_check(items: list[dict], wanted: list[str], args: argparse.Namespace) -> None:
    """Print each sampled item's totals beside its VCF's per-variant maxima, and never fail.

    The two are not the same quantity and are not expected to match; the table exists so a
    total that is off by an order of magnitude, or attached to the wrong phenotype, shows up.
    """
    by_phenocode = {item["phenocode"]: item for item in items}
    cache = Path(args.cache_dir)
    cache.mkdir(parents=True, exist_ok=True)

    print(f"{'phenocode':<22} {'JSON NS':>12} {'VCF max NS':>12} {'JSON NC':>12} {'VCF max NC':>12}")
    for code in wanted:
        item = by_phenocode.get(code)
        if item is None:
            print(f"{code:<22} not in the output")
            continue
        base, _, stratum = code.partition("|")
        name = f"{source_file_code(base)}_variant_meta_analysis_100_cutoff" + (f".{stratum}" if stratum else "") + ".vcf.gz"
        uri = args.variant_prefix.rstrip("/") + "/" + name
        local = fetch_gs(uri, cache, billing_project=args.billing_project)
        max_ns, max_nc = vcf_maxima(local)
        json_ns = item.get("num_samples") or item["num_cases"] + item["num_controls"]
        json_nc = item.get("num_cases")
        print(
            f"{code:<22} {json_ns:>12,} {max_ns or 0:>12,} "
            f"{(f'{json_nc:,}' if json_nc else '-'):>12} {(f'{max_nc:,}' if max_nc else '-'):>12}"
        )


def main() -> None:
    args = parse_args()
    cache = Path(args.cache_dir)

    if args.supplement.startswith(("http://", "https://")):
        workbook = fetch(args.supplement, cache / "brava_supp_tables.xlsx", timeout=300)
    else:
        workbook = Path(args.supplement)
    print(f"Reading {workbook}...")
    tables = read_workbook(workbook)
    print("  " + ", ".join(f"{name}: {len(rows)} rows" for name, rows in tables.items()))
    binary, quantitative = per_ancestry(tables)
    check_totals(tables, binary, quantitative)

    if args.gene_listing:
        lines = Path(args.gene_listing).read_text().splitlines()
        print(f"Reading gene result file names from {args.gene_listing}...")
    else:
        print(f"Listing {args.gene_prefix}...")
        lines = list_gene_files(args.gene_prefix, args.billing_project)
    strata = parse_listing(lines)
    print(f"  {len(strata)} phenotypes, {sum(len(s) for s in strata.values())} gene result files")

    items = build_items(tables, strata, binary, quantitative)
    counts = defaultdict(int)
    for item in items:
        _, _, stratum = item["phenocode"].partition("|")
        counts[stratum or "meta"] += 1
    print(f"  {len(items)} items: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))

    Path(args.output).write_text(json.dumps(items, indent=4) + "\n")
    print(f"wrote {args.output}")

    if args.check is not None:
        cross_check(items, args.check or CHECK_SAMPLE, args)

    if args.stage:
        upload_to_gcs(args.output, STAGE_PATH)


if __name__ == "__main__":
    main()
