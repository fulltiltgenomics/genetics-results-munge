#!/usr/bin/env python3
"""Split the combined Open Targets credible set file into one file per study, with stats.

The combined file is sorted by position, so a study's rows are scattered through it. The rows
are regrouped with an external sort on the accession and read one study at a time: memory is
bounded by the largest study, where reading the whole release into one frame needs more than a
16 GB machine has.

Usage:
    create_open_targets_per_study_files.py COMBINED.tsv.gz OUT_DIR STATS.tsv
"""

import io
import itertools
import os
import subprocess
import sys

import polars as pl

from credible_set_stats import calculate_stats, get_tsv_header, stats_to_tsv_row, write_stats_json

SORT_BUFFER = "2G"


def study_schema(combined: str) -> pl.Schema:
    # one schema for every study, taken as polars infers it for the whole file: left to infer
    # per study, a column that happens to be all NA or all integers in one study would be typed
    # and written differently there than in the next
    return pl.read_csv(combined, separator="\t", null_values=["NA"], n_rows=1000).schema


def rows_by_study(combined: str, key_col: int, tmp_dir: str):
    """yield (accession, raw lines) with each study's rows in their file order"""
    env = {**os.environ, "LC_ALL": "C"}
    unzip = subprocess.Popen(["bgzip", "-dc", combined], stdout=subprocess.PIPE)
    body = subprocess.Popen(["tail", "-n", "+2"], stdin=unzip.stdout, stdout=subprocess.PIPE)
    # stable, so position order survives inside a study
    sort = subprocess.Popen(
        ["sort", "-s", "-t", "\t", f"-k{key_col},{key_col}", "-S", SORT_BUFFER, "-T", tmp_dir],
        stdin=body.stdout,
        stdout=subprocess.PIPE,
        env=env,
    )
    unzip.stdout.close()
    body.stdout.close()
    yield from itertools.groupby(sort.stdout, key=lambda line: line.split(b"\t", key_col)[key_col - 1])
    for proc in (unzip, body, sort):
        if proc.wait() != 0:
            sys.exit(f"{proc.args[0]} exited with {proc.returncode}")


def main(combined: str, out_dir: str, stats_file: str) -> None:
    schema = study_schema(combined)
    header = ("\t".join(schema.names()) + "\n").encode()
    key_col = schema.names().index("trait_original") + 1

    all_stats = []
    for study, lines in rows_by_study(combined, key_col, os.path.dirname(combined) or "."):
        study_data = pl.read_csv(
            io.BytesIO(header + b"".join(lines)), separator="\t", null_values=["NA"], schema=schema
        )
        # files are named by accession: trait carries the free-text study name
        stem = os.path.join(out_dir, f"{study.decode()}.SUSIE.munged")
        study_data.write_csv(f"{stem}.tsv", separator="\t", null_value="NA")
        stats = calculate_stats(study_data)
        write_stats_json(stats, f"{stem}.stats.json")
        all_stats.append(stats)

    with open(stats_file, "w") as f:
        f.write(get_tsv_header() + "\n")
        for s in all_stats:
            f.write(stats_to_tsv_row(s) + "\n")
    print(f"Wrote aggregate stats for {len(all_stats)} studies")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
