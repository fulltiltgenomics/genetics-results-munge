#!/usr/bin/env python3
"""Munge GP2 Parkinson's GWAS summary statistics (build 38) and create AF-AF plot against gnomAD."""

import argparse
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import polars as pl
from scipy.special import log_ndtr

from sumstat_utils import (
    GNOMAD_LOCAL_DEFAULT, chrom_code_expr, create_af_af_plot_by_type, gnomad_plot_parser,
    pick_gnomad_rows, read_gnomad_filtered, save_gnomad_by_variant_key, write_sumstat_output,
)


GNOMAD_PATH = GNOMAD_LOCAL_DEFAULT


def parse_args() -> argparse.Namespace:
    return gnomad_plot_parser(
        __doc__, "Path to GP2 .txt.gz file",
        gnomad_default=GNOMAD_PATH, gnomad_help="Path to gnomAD bgz file", af_plot_flag=False,
    ).parse_args()


def read_sumstats(filepath: str) -> pl.DataFrame:
    """Read GP2 sumstats, parse SNP_ID for ref/alt, compute mlog10p."""
    df = pl.read_csv(
        filepath,
        separator="\t",
        null_values=["NA"],
        truncate_ragged_lines=True,
        schema_overrides={
            "chromosome": pl.Utf8,
            "base_pair_position": pl.Int32,
            "effect_allele_frequency": pl.Float64,
            "beta": pl.Float64,
            "standard_error": pl.Float64,
            "p_value": pl.Float64,
            "N_datasets": pl.Int32,
        },
        ignore_errors=True,
    )

    # a contig with no numeric code casts to null and is dropped
    df = df.with_columns(
        chrom_code_expr("chromosome").cast(pl.Int32, strict=False).alias("chromosome"),
    )
    df = df.filter(pl.col("chromosome").is_not_null())

    # parse ref/alt from SNP_ID (format: chr{chr}:{pos}:{ref}:{alt})
    df = df.with_columns(
        pl.col("SNP_ID").str.split(":").list.get(2).alias("snp_id_ref"),
        pl.col("SNP_ID").str.split(":").list.get(3).alias("snp_id_alt"),
    )

    # classify as SNP or indel
    df = df.with_columns(
        ((pl.col("snp_id_ref").str.len_chars() == 1) & (pl.col("snp_id_alt").str.len_chars() == 1)).alias("is_snp"),
    )

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        df = df.with_columns(
            # mlog10p: use beta/se when p underflows to 0
            pl.when(pl.col("p_value").eq(0) | pl.col("p_value").is_null())
            .then(
                ((-log_ndtr(-(pl.col("beta") / pl.col("standard_error")).abs()) - np.log(2)) / np.log(10)).round(4)
            )
            .otherwise((-np.log10(pl.col("p_value"))).round(4))
            .alias("mlog10p"),
        )

    return df


def join_with_gnomad(
    df: pl.DataFrame, gnomad: pl.DataFrame, af_col: str, gnomad_source: str = "any",
) -> pl.DataFrame:
    """Join with gnomAD on chr:pos:alleles, flip SNPs if effect_allele is ref."""
    gnomad_filt = pick_gnomad_rows(gnomad, gnomad_source)

    # join on chr+pos+ref+alt (SNP_ID alleles match gnomAD since both are build 38)
    joined = df.join(
        gnomad_filt.rename({"chr": "gnomad_chr"}),
        left_on=["chromosome", "base_pair_position", "snp_id_ref", "snp_id_alt"],
        right_on=["gnomad_chr", "pos", "ref", "alt"],
        how="inner",
    )

    # determine allele orientation: is effect_allele the ref or alt?
    joined = joined.with_columns(
        pl.when(
            pl.col("effect_allele").str.to_uppercase() == pl.col("snp_id_alt").str.to_uppercase()
        ).then(pl.lit("effect_is_alt"))
        .when(
            pl.col("effect_allele").str.to_uppercase() == pl.col("snp_id_ref").str.to_uppercase()
        ).then(pl.lit("effect_is_ref"))
        .otherwise(pl.lit("mismatch"))
        .alias("allele_match"),
    )

    n_effect_alt = joined.filter(pl.col("allele_match") == "effect_is_alt").height
    n_effect_ref = joined.filter(pl.col("allele_match") == "effect_is_ref").height
    n_mismatch = joined.filter(pl.col("allele_match") == "mismatch").height
    n_snp = joined.filter(pl.col("is_snp")).height
    n_indel = joined.filter(~pl.col("is_snp")).height
    n_unmatched = df.height - joined.height
    print(f"  gnomAD matched: {joined.height}, unmatched: {n_unmatched}")
    print(f"  SNPs: {n_snp}, indels: {n_indel}")
    print(f"  effect=alt (no flip): {n_effect_alt}, effect=ref (flip SNPs): {n_effect_ref}, mismatch: {n_mismatch}")

    # drop mismatches
    result = joined.filter(pl.col("allele_match") != "mismatch")

    # flip beta and AF for SNPs where effect_allele is ref
    # do NOT flip indels
    is_flip = (pl.col("allele_match") == "effect_is_ref") & pl.col("is_snp")
    n_actually_flipped = result.filter(is_flip).height
    print(f"  actually flipped (SNPs with effect=ref): {n_actually_flipped}")

    result = result.with_columns(
        pl.when(is_flip).then(-pl.col("beta")).otherwise(pl.col("beta")).alias("beta"),
        pl.when(is_flip).then((1.0 - pl.col("effect_allele_frequency")).round(6)).otherwise(pl.col("effect_allele_frequency")).alias("af"),
        # after flip, ref/alt columns should reflect gnomAD orientation
        # for flipped SNPs, swap effect/other allele labels
        pl.when(is_flip).then(pl.col("other_allele")).otherwise(pl.col("effect_allele")).alias("final_alt"),
        pl.when(is_flip).then(pl.col("effect_allele")).otherwise(pl.col("other_allele")).alias("final_ref"),
    )

    # for non-flipped rows, af is just effect_allele_frequency
    result = result.with_columns(
        pl.when(~is_flip).then(pl.col("effect_allele_frequency")).otherwise(pl.col("af")).alias("af"),
    )

    return result


def _prepare_output(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        pl.col("beta").map_elements(lambda x: f"{x:.3e}" if x is not None else None, return_dtype=pl.Utf8).alias("beta"),
        pl.col("standard_error").map_elements(lambda x: f"{x:.3e}" if x is not None else None, return_dtype=pl.Utf8).alias("se"),
    ).select(
        pl.col("chromosome").alias("#chr"),
        pl.col("base_pair_position").alias("pos"),
        pl.col("final_ref").alias("ref"),
        pl.col("final_alt").alias("alt"),
        pl.col("rsids").alias("rsid"),
        "beta", "se", "mlog10p",
        "af",
        pl.col("N_datasets").alias("n_datasets"),
    ).sort(
        "#chr", "pos", "ref", "alt",
    )


def main():
    args = parse_args()

    input_stem = Path(args.input).name.replace(".gz", "")
    output_dir = Path(args.input).parent
    output_path = args.output or str(output_dir / f"{input_stem}.munged.tsv.gz")
    plot_path = args.plot or str(output_dir / f"{input_stem}.af_af.png")
    gnomad_save_path = str(output_dir / f"{input_stem}.gnomad_filtered.tsv")

    print(f"Reading {args.input}...")
    df = read_sumstats(args.input)
    print(f"  {df.height} variants loaded")
    print(f"  SNPs: {df.filter(pl.col('is_snp')).height}, indels: {df.filter(~pl.col('is_snp')).height}")

    n_p_zero = df.filter(pl.col("p_value").eq(0) | pl.col("p_value").is_null()).height
    print(f"  p-value underflows (mlog10p from beta/se): {n_p_zero}")

    if args.gnomad_filtered:
        print(f"Reading filtered gnomAD from {args.gnomad_filtered}...")
        gnomad = read_gnomad_filtered(args.gnomad_filtered)
    else:
        # build variant keys for gnomAD matching
        print("Building variant key set...")
        keys = set(
            df.select(
                (pl.col("chromosome").cast(pl.Utf8) + ":" +
                 pl.col("base_pair_position").cast(pl.Utf8) + ":" +
                 pl.col("snp_id_ref") + ":" +
                 pl.col("snp_id_alt")).alias("key")
            )["key"].to_list()
        )
        print(f"  {len(keys)} unique variant keys")

        print(f"Streaming gnomAD from {args.gnomad}...")
        save_gnomad_by_variant_key(args.gnomad, keys, gnomad_save_path)
        gnomad = read_gnomad_filtered(gnomad_save_path)

    print(f"  {gnomad.height} gnomAD rows loaded")

    print("Joining with gnomAD...")
    df = join_with_gnomad(df, gnomad, args.af_col, args.gnomad_source)

    print("Creating AF-AF plot...")
    create_af_af_plot_by_type(df, args.af_col, plot_path, input_stem, "GP2 alt AF (after flip)")

    print("Writing output...")
    write_sumstat_output(_prepare_output(df), output_path)

    print("Done.")


if __name__ == "__main__":
    main()
