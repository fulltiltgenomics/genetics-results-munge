#!/usr/bin/env python3
"""Munge SCHEMA2 variant results into genebass-compatible format.

SCHEMA2 provides allele counts but no p-values or effect sizes, so we compute
log-odds ratios with Haldane correction and derive p-values from z-scores.
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import polars as pl
from scipy.special import log_ndtr

from sumstat_utils import write_exome_output

CODING_CONSEQUENCES = {
    "missense_variant", "frameshift_variant",
    "inframe_deletion", "inframe_insertion",
    "stop_gained", "stop_lost", "stop_retained_variant", "start_lost",
    "splice_donor_variant", "splice_acceptor_variant",
    "coding_sequence_variant", "protein_altering_variant",
}

MLOG10P_THRESHOLD = 4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to SCHEMA2 variant results .bgz file")
    parser.add_argument("--gencode", default="/mnt/disks/data/gencode.v39.annotation.genes.tsv",
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


def read_schema2_variants(path: str) -> pl.DataFrame:
    """Read SCHEMA2 variant results, filter to coding variants."""
    print("  scanning file...")
    df = (
        pl.scan_csv(
            path, separator="\t", null_values=["NA"],
            schema_overrides={
                "ac_case": pl.Int64, "an_case": pl.Int64,
                "ac_ctrl": pl.Int64, "an_ctrl": pl.Int64,
                "n_de_novo": pl.Int64, "mpc": pl.Float64,
            },
        )
        .filter(pl.col("consequence").is_in(CODING_CONSEQUENCES))
        .collect()
    )
    print(f"  {df.height} coding rows")

    # parse locus "chrN:POS" -> chr, pos
    df = df.with_columns(
        pl.col("locus").str.split(":").list.first()
            .str.replace("^chr", "")
            .str.replace("^X$", "23").str.replace("^Y$", "24").str.replace("^MT?$", "26")
            .cast(pl.Int32).alias("chr"),
        pl.col("locus").str.split(":").list.last().cast(pl.Int64).alias("pos"),
    )

    # parse alleles '["REF","ALT"]' -> ref, alt
    df = df.with_columns(
        pl.col("alleles").str.strip_chars("[]").str.replace_all('"', '').str.split(","),
    ).with_columns(
        pl.col("alleles").list.first().alias("ref"),
        pl.col("alleles").list.last().alias("alt"),
    )

    return df


def compute_stats(df: pl.DataFrame) -> pl.DataFrame:
    """Compute log-OR, SE, mlog10p, and allele frequencies from counts.

    Uses Haldane correction (add 0.5) to handle zero cells in the 2x2 table.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        ac_case = df["ac_case"].to_numpy(zero_copy_only=False).astype(np.float64)
        an_case = df["an_case"].to_numpy(zero_copy_only=False).astype(np.float64)
        ac_ctrl = df["ac_ctrl"].to_numpy(zero_copy_only=False).astype(np.float64)
        an_ctrl = df["an_ctrl"].to_numpy(zero_copy_only=False).astype(np.float64)

        # Haldane-corrected cells
        a = ac_case + 0.5
        b = (an_case - ac_case) + 0.5
        c = ac_ctrl + 0.5
        d = (an_ctrl - ac_ctrl) + 0.5

        beta = np.log(a * d / (b * c))
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

        # mlog10p from z-score: -log10(2 * Phi(-|z|))
        z_abs = np.abs(beta / se)
        # log_ndtr gives log(Phi(x)), so log(p) = log_ndtr(-|z|) + log(2)
        log10p = (log_ndtr(-z_abs) + np.log(2)) / np.log(10)
        mlog10p = np.round(np.maximum(-log10p, 0.0), 4)

        # mask rows where an is 0
        valid = (an_case > 0) & (an_ctrl > 0)
        beta = np.where(valid, beta, np.nan)
        se = np.where(valid, se, np.nan)
        mlog10p = np.where(valid, mlog10p, np.nan)

        # allele frequencies
        af_overall = np.where(
            (an_case + an_ctrl) > 0,
            (ac_case + ac_ctrl) / (an_case + an_ctrl),
            np.nan,
        )
        af_cases = np.where(an_case > 0, ac_case / an_case, np.nan)
        af_controls = np.where(an_ctrl > 0, ac_ctrl / an_ctrl, np.nan)

        df = df.with_columns(
            pl.Series("beta", beta, dtype=pl.Float64),
            pl.Series("se", se, dtype=pl.Float64),
            pl.Series("mlog10p", mlog10p, dtype=pl.Float64),
            pl.Series("af_overall", af_overall, dtype=pl.Float64),
            pl.Series("af_cases", af_cases, dtype=pl.Float64),
            pl.Series("af_controls", af_controls, dtype=pl.Float64),
            (pl.col("ac_case") + pl.col("ac_ctrl")).alias("ac"),
            (pl.col("an_case") + pl.col("an_ctrl")).alias("an"),
        )

        # format numeric columns as scientific notation
        for col in ["beta", "se", "af_overall", "af_cases", "af_controls"]:
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
        pl.lit("SCHEMA2").alias("#dataset"),
        "chr", "pos", "ref", "alt",
        "gene",
        pl.col("consequence").alias("annotation"),
        "mlog10p",
        "beta",
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
        Path(args.input).parent / "SCHEMA2_variant_results.munged.tsv.gz"
    )

    print(f"Reading gencode from {args.gencode}...")
    gencode = read_gencode(args.gencode)
    print(f"  {gencode.height} genes")

    print(f"Reading SCHEMA2 variants from {args.input}...")
    df = read_schema2_variants(args.input)

    print("Computing stats...")
    df = compute_stats(df)

    print("Building output...")
    out = build_output(df, gencode)

    print("Writing results...")
    write_exome_output(out, output_path, tabix_args=["-s2", "-b3", "-e3"])

    print("Done.")


if __name__ == "__main__":
    main()
