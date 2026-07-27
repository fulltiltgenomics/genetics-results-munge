#!/usr/bin/env python3
"""Munge SCHEMA variant results into genebass-compatible format."""

import argparse
import warnings
from pathlib import Path

import numpy as np
import polars as pl
from scipy.special import log_ndtr

from sumstat_utils import write_exome_output


CODING_CONSEQUENCES = {
    "missense_variant_mpc_<2", "missense_variant_mpc_2-3", "missense_variant_mpc_>=3",
    "frameshift_variant", "inframe_deletion", "inframe_insertion",
    "stop_gained", "stop_lost", "stop_retained_variant", "start_lost",
    "splice_donor_variant", "splice_acceptor_variant",
    "coding_sequence_variant", "protein_altering_variant",
}

MLOG10P_THRESHOLD = 4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to SCHEMA variant results .bgz file")
    parser.add_argument("--gencode", default="/mnt/disks/data/gencode.v19.annotation.genes.tsv",
                        help="Path to gencode annotation genes TSV")
    parser.add_argument("--output", help="Output munged TSV path (default: derived from input)")
    return parser.parse_args()


def read_gencode(path: str) -> pl.DataFrame:
    """Read gencode annotation, strip version from gene_id."""
    df = pl.read_csv(path, separator="\t")
    df = df.with_columns(
        pl.col("gene_id").str.split(".").list.first().alias("gene_id_base"),
        pl.col("gene_name").alias("gene"),
    ).select("gene_id_base", "gene")
    return df.unique(subset=["gene_id_base"], keep="first")


def read_schema_variants(path: str) -> pl.DataFrame:
    """Read SCHEMA variant results, filter to meta_super coding variants."""
    print("  scanning file...")
    df = (
        pl.scan_csv(
            path, separator="\t", null_values=["NA"],
            schema_overrides={"p": pl.Float64, "est": pl.Float64, "se": pl.Float64,
                              "ac_case": pl.Int64, "ac_ctrl": pl.Int64,
                              "an_case": pl.Int64, "an_ctrl": pl.Int64},
        )
        .filter(pl.col("group") == "meta_super")
        .filter(pl.col("consequence").is_in(CODING_CONSEQUENCES))
        .collect()
    )
    print(f"  {df.height} coding meta_super rows")

    # parse locus "CHR:POS" → chr, pos
    df = df.with_columns(
        pl.col("locus").str.split(":").list.first()
            .str.replace("X", "23").str.replace("Y", "24").str.replace("MT", "26")
            .cast(pl.Int32).alias("chr"),
        pl.col("locus").str.split(":").list.last().cast(pl.Int32).alias("pos"),
    )

    # parse alleles '["REF","ALT"]' → ref, alt
    df = df.with_columns(
        pl.col("alleles").str.strip_chars("[]").str.replace_all('"', '').str.split(","),
    ).with_columns(
        pl.col("alleles").list.first().alias("ref"),
        pl.col("alleles").list.last().alias("alt"),
    )

    return df


def compute_stats(df: pl.DataFrame) -> pl.DataFrame:
    """Compute mlog10p, allele frequencies, and format columns."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        # allele frequencies
        df = df.with_columns(
            ((pl.col("ac_case") + pl.col("ac_ctrl")) / (pl.col("an_case") + pl.col("an_ctrl"))).alias("af_overall"),
            (pl.col("ac_case") / pl.col("an_case")).alias("af_cases"),
            (pl.col("ac_ctrl") / pl.col("an_ctrl")).alias("af_controls"),
            (pl.col("ac_case") + pl.col("ac_ctrl")).alias("ac"),
            (pl.col("an_case") + pl.col("an_ctrl")).alias("an"),
        )

        # mlog10p from beta/se, fall back to raw p
        df = df.with_columns(
            pl.when(pl.col("est").is_not_null() & pl.col("se").is_not_null() & (pl.col("se") > 0))
                .then(
                    pl.max_horizontal(
                        ((-log_ndtr(-(pl.col("est") / pl.col("se")).abs()) - np.log(2)) / np.log(10)).round(4),
                        0.0,
                    )
                )
                .when(pl.col("p").is_not_null() & (pl.col("p") > 0))
                .then(pl.max_horizontal((-np.log10(pl.col("p"))).round(4), 0.0))
                .when(pl.col("p").is_not_null() & pl.col("p").eq(0))
                .then(324.0)
                .otherwise(None)
                .alias("mlog10p"),
        )

        # format numeric columns as scientific notation
        for col in ["est", "se", "af_overall", "af_cases", "af_controls"]:
            df = df.with_columns(
                pl.col(col).map_elements(lambda x: f"{x:.4e}", return_dtype=pl.Utf8).alias(col),
            )

    return df


def build_output(df: pl.DataFrame, gencode: pl.DataFrame) -> pl.DataFrame:
    """Join with gencode and build final output."""
    n_before = df.height
    df = df.join(gencode, left_on="gene_id", right_on="gene_id_base", how="inner")
    n_unmatched = n_before - df.height
    if n_unmatched > 0:
        print(f"  warning: {n_unmatched} variants not matched to gencode")

    out = df.select(
        pl.lit("SCHEMA").alias("#dataset"),
        "chr", "pos", "ref", "alt",
        "gene",
        pl.col("consequence").alias("annotation"),
        "mlog10p",
        pl.col("est").alias("beta"),
        "se",
        "af_overall", "af_cases", "af_controls",
        "ac", "an",
        pl.lit(None, dtype=pl.Utf8).alias("n_cases"),
        pl.lit(None, dtype=pl.Utf8).alias("n_controls"),
        pl.lit("schizophrenia").alias("trait"),
        pl.lit("schizophrenia").alias("trait_original"),
    ).sort("chr", "pos", "ref", "alt")

    return out


def main():
    args = parse_args()

    output_path = args.output or str(
        Path(args.input).parent / "SCHEMA_variant_results.munged.tsv.gz"
    )

    print(f"Reading gencode from {args.gencode}...")
    gencode = read_gencode(args.gencode)
    print(f"  {gencode.height} genes")

    print(f"Reading SCHEMA variants from {args.input}...")
    df = read_schema_variants(args.input)

    print("Computing stats...")
    df = compute_stats(df)

    print("Building output...")
    out = build_output(df, gencode)

    print("Writing results...")
    write_exome_output(out, output_path, tabix_args=["-s2", "-b3", "-e3"])

    print("Done.")


if __name__ == "__main__":
    main()
