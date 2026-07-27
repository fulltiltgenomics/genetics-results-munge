#!/usr/bin/env python3
"""Munge COVID-19 HGI freeze 7 GWAS summary statistics.

Input: COVID19_HGI_*_ALL_leave_23andme_20220403.tsv.gz (GRCh38)
- ALT is the effect allele, no flipping needed
- Streams gnomAD to match by chr:pos:ref:alt for variant consequence and gene name
- Computes -log10(p) from beta/se when p-value underflows
"""

import argparse
import subprocess
import warnings
from pathlib import Path

import numpy as np
import polars as pl
from scipy.special import log_ndtr

from sumstat_utils import write_sumstat_output


GNOMAD_PATH = "/mnt/disks/data/gnomad/gnomad.genomes.exomes.v4.0.sites.v2.tsv.bgz"

GNOMAD_KEEP_COLS = [
    "#chr", "pos", "ref", "alt", "rsids",
    "most_severe", "gene_most_severe", "genome_or_exome",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to COVID HGI .tsv.gz file")
    parser.add_argument("--gnomad", default=GNOMAD_PATH, help="Path to gnomAD file (local bgz or plain)")
    parser.add_argument("--gnomad-source", default="g", choices=["g", "e"],
                        help="Use gnomAD genomes (g) or exomes (e) (default: g)")
    parser.add_argument("--gnomad-filtered", help="Path to previously saved filtered gnomAD TSV (skip streaming)")
    parser.add_argument("--output", help="Output munged TSV path (default: derived from input)")
    return parser.parse_args()


def read_sumstats(filepath: str) -> pl.DataFrame:
    """Read COVID HGI sumstats and compute mlog10p."""
    df = pl.read_csv(
        filepath,
        separator="\t",
        null_values=["NA"],
        truncate_ragged_lines=True,
        schema_overrides={
            "#CHR": pl.Utf8,
            "POS": pl.Int32,
            "all_inv_var_meta_beta": pl.Float64,
            "all_inv_var_meta_sebeta": pl.Float64,
            "all_inv_var_meta_p": pl.Float64,
            "all_inv_var_meta_cases": pl.Int32,
            "all_inv_var_meta_controls": pl.Int32,
            "all_inv_var_meta_effective": pl.Float64,
            "all_meta_AF": pl.Float64,
        },
        ignore_errors=True,
    )

    df = df.with_columns(
        pl.col("#CHR").str.replace("X", "23").cast(pl.Int32, strict=False).alias("chr"),
    )
    df = df.filter(pl.col("chr").is_not_null())

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        df = df.with_columns(
            pl.when(pl.col("all_inv_var_meta_p").eq(0) | pl.col("all_inv_var_meta_p").is_null())
            .then(
                ((-log_ndtr(-(pl.col("all_inv_var_meta_beta") / pl.col("all_inv_var_meta_sebeta")).abs())
                  - np.log(2)) / np.log(10)).round(4)
            )
            .otherwise((-np.log10(pl.col("all_inv_var_meta_p"))).round(4))
            .alias("mlog10p"),
        )

    return df


def _build_variant_keys(df: pl.DataFrame) -> set[str]:
    """Build set of chr:pos:ref:alt keys for gnomAD matching."""
    return set(
        df.select(
            (pl.col("chr").cast(pl.Utf8) + ":" +
             pl.col("POS").cast(pl.Utf8) + ":" +
             pl.col("REF") + ":" +
             pl.col("ALT")).alias("key")
        )["key"].to_list()
    )


def stream_gnomad(gnomad_path: str, variant_keys: set[str], save_path: str) -> pl.DataFrame:
    """Stream gnomAD, keeping rows matching variant_keys (chr:pos:ref:alt).
    Saves filtered rows to save_path for reuse."""
    proc = subprocess.Popen(
        f"zcat '{gnomad_path}'",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=4_000_000,
    )

    header_line = proc.stdout.readline().strip()
    header = header_line.split("\t")
    col_idx = {}
    for col in GNOMAD_KEEP_COLS:
        if col in header:
            col_idx[col] = header.index(col)

    chr_idx = col_idx["#chr"]
    pos_idx = col_idx["pos"]
    ref_idx = col_idx["ref"]
    alt_idx = col_idx["alt"]

    out_header = "\t".join(col_idx.keys())
    outf = open(save_path, "w")
    outf.write(out_header + "\n")

    n_lines = 0
    n_matches = 0
    for line in proc.stdout:
        n_lines += 1
        if n_lines % 10_000_000 == 0:
            print(f"  gnomAD: {n_lines / 1e6:.0f}M lines scanned, {n_matches} matches", flush=True)

        fields = line.split("\t")
        if len(fields) <= alt_idx:
            continue

        key = f"{fields[chr_idx]}:{fields[pos_idx]}:{fields[ref_idx]}:{fields[alt_idx]}"
        if key in variant_keys:
            out_fields = [fields[idx].rstrip("\n") for idx in col_idx.values()]
            outf.write("\t".join(out_fields) + "\n")
            n_matches += 1

    outf.close()
    proc.wait()
    print(f"  gnomAD: done. {n_lines / 1e6:.1f}M lines scanned, {n_matches} matches", flush=True)
    print(f"  saved filtered gnomAD to {save_path}")

    return _read_gnomad_filtered(save_path)


def _read_gnomad_filtered(filepath: str) -> pl.DataFrame:
    """Read a previously saved filtered gnomAD TSV."""
    schema_overrides = {"#chr": pl.Int32, "pos": pl.Int32}
    if filepath.endswith(".gz"):
        proc = subprocess.Popen(["zcat", filepath], stdout=subprocess.PIPE)
        df = pl.read_csv(
            proc.stdout, separator="\t", null_values=["NA"],
            schema_overrides=schema_overrides, ignore_errors=True,
        )
        proc.wait()
    else:
        df = pl.read_csv(
            filepath, separator="\t", null_values=["NA"],
            schema_overrides=schema_overrides, ignore_errors=True,
        )
    return df.rename({"#chr": "chr"})


def join_with_gnomad(df: pl.DataFrame, gnomad: pl.DataFrame, gnomad_source: str) -> pl.DataFrame:
    """Join sumstats with gnomAD on chr:pos:ref:alt to get consequence/gene."""
    gnomad_filt = gnomad.filter(pl.col("genome_or_exome") == gnomad_source)
    source_label = "genomes" if gnomad_source == "g" else "exomes"
    print(f"  using gnomAD {source_label}: {gnomad_filt.height} rows (from {gnomad.height} total)")

    joined = df.join(
        gnomad_filt.select(
            "chr", "pos",
            pl.col("ref").alias("gnomad_ref"),
            pl.col("alt").alias("gnomad_alt"),
            "rsids", "most_severe", "gene_most_severe",
        ),
        left_on=["chr", "POS", "REF", "ALT"],
        right_on=["chr", "pos", "gnomad_ref", "gnomad_alt"],
        how="left",
    )

    n_matched = joined.filter(pl.col("most_severe").is_not_null()).height
    print(f"  gnomAD matched: {n_matched} / {joined.height}")

    return joined


def _prepare_output(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        pl.col("all_inv_var_meta_beta")
            .map_elements(lambda x: f"{x:.3e}" if x is not None else None, return_dtype=pl.Utf8)
            .alias("beta"),
        pl.col("all_inv_var_meta_sebeta")
            .map_elements(lambda x: f"{x:.3e}" if x is not None else None, return_dtype=pl.Utf8)
            .alias("se"),
    ).select(
        pl.col("chr").alias("#chr"),
        pl.col("POS").alias("pos"),
        pl.col("REF").alias("ref"),
        pl.col("ALT").alias("alt"),
        pl.coalesce(pl.col("rsids"), pl.col("rsid")).alias("rsid"),
        "beta", "se", "mlog10p",
        pl.col("all_meta_AF").alias("af"),
        pl.col("all_inv_var_meta_cases").alias("n_cases"),
        pl.col("all_inv_var_meta_controls").alias("n_controls"),
        pl.col("all_inv_var_meta_effective").alias("neff"),
        pl.col("all_meta_N").alias("n_studies"),
        "most_severe", "gene_most_severe",
    ).sort(
        "#chr", "pos", "ref", "alt",
    )


def main():
    args = parse_args()

    input_path = Path(args.input)
    output_dir = input_path.parent
    output_path = args.output or str(output_dir / f"{input_path.stem}.munged.tsv.gz")
    gnomad_save_path = str(output_dir / f"{input_path.stem}.gnomad_filtered.tsv")

    print(f"Reading {args.input}...")
    df = read_sumstats(args.input)
    print(f"  {df.height} variants loaded")

    n_p_zero = df.filter(pl.col("all_inv_var_meta_p").eq(0) | pl.col("all_inv_var_meta_p").is_null()).height
    print(f"  p-value underflows (mlog10p from beta/se): {n_p_zero}")

    if args.gnomad_filtered:
        print(f"Reading filtered gnomAD from {args.gnomad_filtered}...")
        gnomad = _read_gnomad_filtered(args.gnomad_filtered)
    else:
        print("Building variant key set...")
        keys = _build_variant_keys(df)
        print(f"  {len(keys)} unique variant keys")

        print(f"Streaming gnomAD from {args.gnomad}...")
        gnomad = stream_gnomad(args.gnomad, keys, gnomad_save_path)

    print(f"  {gnomad.height} gnomAD rows loaded")

    print("Joining with gnomAD...")
    df = join_with_gnomad(df, gnomad, args.gnomad_source)

    print("Writing output...")
    write_sumstat_output(_prepare_output(df), output_path)

    print("Done.")


if __name__ == "__main__":
    main()
