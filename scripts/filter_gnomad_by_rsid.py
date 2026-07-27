#!/usr/bin/env python3
"""Filter gnomAD to variants matching rsids from a PGC sumstat file.

Reads only the rsid column into memory, then streams gnomAD line by line,
writing matching rows to a local bgzipped file.

Usage:
    python scripts/filter_gnomad_by_rsid.py \
        --input ~/PGC/daner_PGC_SCZ_w3_90_0418b.gz \
        --output ~/PGC/gnomad_filtered.tsv.gz
"""

import argparse
import gzip
import subprocess
import sys

GNOMAD_DEFAULT = "gs://finngen-commons/results_api_data/gnomad/gnomad.genomes.exomes.v4.0.sites.v2.tsv.bgz"


def load_rsids(filepath: str) -> set[str]:
    """Load rsid column from a PGC daner file into a set."""
    rsids = set()
    with gzip.open(filepath, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        snp_idx = header.index("SNP")
        for line in f:
            fields = line.split("\t", snp_idx + 2)
            if len(fields) <= snp_idx:
                continue
            rsid = fields[snp_idx]
            if rsid:
                rsids.add(rsid)
    print(f"Loaded {len(rsids)} unique rsids from {filepath}", flush=True)
    return rsids


def stream_filter(gnomad_path: str, rsids: set[str], output_path: str) -> None:
    """Stream gnomAD, writing rows with matching rsids to output."""
    reader = subprocess.Popen(
        f"gsutil cat {gnomad_path} | zcat",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=4_000_000,
    )

    writer = subprocess.Popen(
        ["bgzip", "-c"],
        stdin=subprocess.PIPE,
        stdout=open(output_path, "wb"),
        text=True,
    )

    header = reader.stdout.readline()
    writer.stdin.write(header)

    header_fields = header.rstrip("\n").split("\t")
    rsid_idx = header_fields.index("rsids")

    n_lines = 0
    n_matches = 0
    for line in reader.stdout:
        n_lines += 1
        if n_lines % 10_000_000 == 0:
            print(f"  {n_lines / 1e6:.0f}M lines scanned, {n_matches} matches", flush=True)

        # rsids column can contain comma-separated values
        fields = line.split("\t", rsid_idx + 2)
        if len(fields) <= rsid_idx:
            continue
        rsid_field = fields[rsid_idx]
        if rsid_field and rsids.intersection(rsid_field.split(",")):
            writer.stdin.write(line)
            n_matches += 1

    writer.stdin.close()
    writer.wait()
    reader.wait()
    print(f"Done. {n_lines / 1e6:.1f}M lines scanned, {n_matches} matches written to {output_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="PGC daner .gz file (rsids read from SNP column)")
    parser.add_argument("--gnomad", default=GNOMAD_DEFAULT, help="GCS path to gnomAD bgz file")
    parser.add_argument("--output", required=True, help="Output bgzipped TSV path")
    args = parser.parse_args()

    rsids = load_rsids(args.input)
    stream_filter(args.gnomad, rsids, args.output)


if __name__ == "__main__":
    main()
