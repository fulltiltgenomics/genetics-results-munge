#!/usr/bin/env python3
"""Rewrite the `trait` / `trait_original` columns of a munged exome result file in place of
re-running its munge.

The served phenotype code of a per-trait burden file is its file name, and the API, the
search index and BigQuery all key on `trait_original`, so a dataset whose munge wrote one
spelling into the rows and whose registry advertises another is unreachable by code until
the two agree. This relabels the rows of an existing `.munged.tsv.gz` (gene burden or
variant layout, told apart by the header) without touching any other column, re-bgzips and
re-indexes it, so the raw inputs need not be on hand.

Every row's `trait_original` must be in the map; an unmapped value stops the run rather
than letting one spelling survive in a file meant to carry the other.

Usage:
  scripts/relabel_traits.py --map crohns_disease=CD="Crohn's disease" IN.tsv.gz OUT.tsv.gz
"""

import argparse
import subprocess
import sys
from pathlib import Path

GENE_TABIX = ["-s5", "-b6", "-e6"]
VARIANT_TABIX = ["-s2", "-b3", "-e3"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="bgzipped TSV with a #-prefixed header")
    parser.add_argument("output", help="where to write the relabelled bgzipped TSV (and its .tbi)")
    parser.add_argument("--map", action="append", required=True, metavar="OLD=CODE=NAME",
                        help="trait_original OLD becomes CODE and trait becomes NAME; repeatable")
    return parser.parse_args()


def parse_map(specs: list[str]) -> dict[str, tuple[str, str]]:
    mapping = {}
    for spec in specs:
        parts = spec.split("=", 2)
        if len(parts) != 3 or not all(parts):
            raise SystemExit(f"--map {spec!r} is not OLD=CODE=NAME")
        mapping[parts[0]] = (parts[1], parts[2])
    return mapping


def tabix_args_for(header: list[str]) -> list[str]:
    if "gene_start_pos" in header:
        return GENE_TABIX
    if "pos" in header:
        return VARIANT_TABIX
    raise SystemExit(f"header carries neither gene_start_pos nor pos: {header}")


def main() -> None:
    args = parse_args()
    mapping = parse_map(args.map)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_rows = 0
    counts: dict[str, int] = {}
    with subprocess.Popen(["gzip", "-cd", args.input], stdout=subprocess.PIPE) as reader, \
         open(out_path, "wb") as out, \
         subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=out) as writer:
        header_line = reader.stdout.readline()
        if not header_line.startswith(b"#"):
            raise SystemExit(f"{args.input} does not start with a #-prefixed header line")
        header = header_line.rstrip(b"\n").lstrip(b"#").decode().split("\t")
        trait_i, orig_i = header.index("trait"), header.index("trait_original")
        tabix_args = tabix_args_for(header)
        writer.stdin.write(header_line)
        for line in reader.stdout:
            fields = line.rstrip(b"\n").split(b"\t")
            old = fields[orig_i].decode()
            try:
                code, name = mapping[old]
            except KeyError:
                raise SystemExit(f"{args.input}: trait_original {old!r} is not in the map") from None
            fields[orig_i] = code.encode()
            fields[trait_i] = name.encode()
            writer.stdin.write(b"\t".join(fields) + b"\n")
            n_rows += 1
            counts[code] = counts.get(code, 0) + 1
        writer.stdin.close()
        if reader.wait() != 0 or writer.wait() != 0:
            raise SystemExit("gzip or bgzip failed")

    subprocess.run(["tabix", "-f", *tabix_args, str(out_path)], check=True)
    print(f"{out_path}: {n_rows} rows, trait_original counts {counts}", file=sys.stderr)


if __name__ == "__main__":
    main()
