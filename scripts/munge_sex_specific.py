#!/usr/bin/env python3
"""Write the phenotype metadata JSON of the FinnGen R12 sex-specific GWAS.

Source (FinnGen green library, R12 analysis data, ``sex_difference/``):
  - ``R12_sex_specific_N_cases_controls_samples.tsv``
        female and male cases/controls/samples for the 2,207 analysed endpoints; cases
        and controls are NA for the three quantitative endpoints
  - the R12 core phenotype metadata JSON (``finngen_r12_pheno.json``), for names

The output is the dataset's phenotype metadata in the pheweb shape results-api reads: the
2,207 analysed endpoints only, so the 262 single-sex endpoints cannot surface as traits
with no data, with ``num_cases`` and ``num_controls`` summed over the sexes (the combined
meta sample) and the quantitative three carrying their sample count as ``num_cases`` with
``num_controls`` 0, exactly as the R12 core file does.

The summary statistics themselves are served in place from the library and the release
goes nowhere near BigQuery, so this is the whole munge.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

NA = "NA"

_N_COLUMNS = [
    "female_N_cases", "female_N_controls", "female_N_samples",
    "male_N_cases", "male_N_controls", "male_N_samples",
]


def read_sex_counts(path: Path) -> dict[str, dict[str, str]]:
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = [c for c in ["PHENO", *_N_COLUMNS] if c not in reader.fieldnames]
        if missing:
            raise SystemExit(f"{path}: missing columns {missing}")
        return {row["PHENO"]: row for row in reader}


def read_core_names(path: Path) -> dict[str, dict]:
    with open(path) as f:
        items = json.load(f)
    return {item["phenocode"]: item for item in items}


def _count(value: str) -> int | None:
    return None if value == NA else int(value)


def write_phenotypes(counts: dict, names: dict, out: Path) -> int:
    items = []
    unnamed = []
    for code, n in counts.items():
        core = names.get(code)
        if core is None:
            unnamed.append(code)
            continue
        f_cases, m_cases = _count(n["female_N_cases"]), _count(n["male_N_cases"])
        if f_cases is None or m_cases is None:
            # quantitative: the whole sample is the analysis, in the core file's shape
            num_cases = int(n["female_N_samples"]) + int(n["male_N_samples"])
            num_controls = 0
        else:
            num_cases = f_cases + m_cases
            num_controls = int(n["female_N_controls"]) + int(n["male_N_controls"])
        item = {
            "phenocode": code,
            "phenostring": core.get("phenostring", ""),
            "num_cases": num_cases,
            "num_controls": num_controls,
            "category": core.get("category", ""),
        }
        if "category_index" in core:
            item["category_index"] = core["category_index"]
        items.append(item)
    if unnamed:
        # an unnamed endpoint would surface as a trait with no name or category; the
        # 2024-05-06 counts file has none, so this is a stop rather than a drop
        raise SystemExit(f"{len(unnamed)} analysed endpoints missing from the core metadata: {unnamed[:10]}")
    items.sort(key=lambda i: i["phenocode"])
    with open(out, "w") as f:
        json.dump(items, f, indent=4)
        f.write("\n")
    return len(items)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sex_counts", type=Path, help="R12_sex_specific_N_cases_controls_samples.tsv")
    parser.add_argument("--pheno-metadata", type=Path, required=True, help="the R12 core finngen_r12_pheno.json")
    parser.add_argument("--output", type=Path, required=True, help="output directory")
    args = parser.parse_args()

    counts = read_sex_counts(args.sex_counts)
    names = read_core_names(args.pheno_metadata)
    args.output.mkdir(parents=True, exist_ok=True)

    pheno = args.output / "finngen_r12_sex_specific_pheno.json"
    n = write_phenotypes(counts, names, pheno)
    print(f"{pheno}: {n} endpoints ({len(counts)} analysed)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
