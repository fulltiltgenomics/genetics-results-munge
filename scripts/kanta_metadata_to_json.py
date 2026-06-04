#!/usr/bin/env python3
"""Convert Kanta metadata TSV to JSON in a format similar to finngen_r13_pheno_202509.json."""

import argparse
import csv
import json
import sys


def parse_row(row):
    analysis_type = row["AnalysisType"]
    n_cases = int(row["N_cases"])
    n_controls = int(row["N_controls"])

    entry = {
        "phenocode": row["OMOPID"],
        "phenostring": row["phenostring"],
        "num_cases": n_cases,
        "num_controls": n_controls,
        "num_total": int(row["N_total"]),
        "category": analysis_type,
        "unit": row["unit"],
    }

    if analysis_type != "Binary":
        entry["mean"] = float(row["Mean"]) if row["Mean"] != "-" else None
        entry["sd"] = float(row["SD"]) if row["SD"] != "-" else None

    if row["Outlier_filter"] != "-":
        entry["outlier_filter"] = row["Outlier_filter"]

    if row["case_def_binary"] not in ("-", "NA", ""):
        entry["case_def_binary"] = row["case_def_binary"]

    if row["Note"] not in ("-", "NA", ""):
        entry["note"] = row["Note"]

    return entry


def main():
    parser = argparse.ArgumentParser(
        description="Convert Kanta metadata TSV to JSON"
    )
    parser.add_argument("input_tsv", help="Input TSV file")
    parser.add_argument("output_json", help="Output JSON file")
    args = parser.parse_args()

    entries = []
    with open(args.input_tsv, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            entries.append(parse_row(row))

    with open(args.output_json, "w") as f:
        json.dump(entries, f)

    print(f"Wrote {len(entries)} entries to {args.output_json}", file=sys.stderr)


if __name__ == "__main__":
    main()
