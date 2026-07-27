#!/usr/bin/env python3
"""Munge IBD Supplementary Table ST4 (variant results) into genebass-compatible variant format.

Input has per-trait columns (P_meta_CD, BETA_meta_CD, P_meta_UC, BETA_meta_UC, P_meta_IBD, BETA_meta_IBD).
Outputs three files (IBD, CD, UC), each containing all variants with trait-specific stats.
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import norm

from sumstat_utils import write_exome_output

GROUPS = ["IBD", "CD", "UC"]

TRAIT_MAP = {
    "IBD": "inflammatory_bowel_disease",
    "CD": "crohns_disease",
    "UC": "ulcerative_colitis",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True,
                        help="Path to Supplementary_Tables_-_ST4.83_Tier1_8_LD_unresolved.tsv")
    parser.add_argument("--input-extra", nargs="*", default=[],
                        help="Additional supplementary tables with same format (e.g. ST5)")
    parser.add_argument("--output-dir", help="Output directory (default: same as input)")
    parser.add_argument("--n-cases", required=True, nargs=3, type=int, metavar=("IBD", "CD", "UC"),
                        help="Number of cases for IBD, CD, UC")
    parser.add_argument("--n-controls", required=True, type=int, help="Number of controls (shared)")
    return parser.parse_args()


# columns used by downstream processing
KEEP_COLS = [
    "chr", "pos", "ref", "alt", "Type", "Notes", "ID", "PIP",
    "Most_Associated_Trait", "gnomADv4.1_exome_nfe_frq", "gnomADv4.1_genome_nfe_frq",
    "Gene_symbol", "most_severe_consequence", "AAC",
    "P_meta_CD", "BETA_meta_CD", "Direction_CD", "HetP_CD", "P_gwas_tier1_CD",
    "P_meta_UC", "BETA_meta_UC", "Direction_UC", "HetP_UC", "P_gwas_tier1_UC",
    "P_meta_IBD", "BETA_meta_IBD", "Direction_IBD", "HetP_IBD", "P_gwas_tier1_IBD",
    "p_gwas_tier2_eur_cd", "p_gwas_tier2_eur_eas_cd",
    "p_gwas_tier2_eur_uc", "p_gwas_tier2_eur_eas_uc",
    "p_gwas_tier2_eur_ibd", "p_gwas_tier2_eur_eas_ibd",
    "GWAS_index", "FM_coding", "22CD_reported",
]


def read_variants(path: str) -> pl.DataFrame:
    """Read variant file, keeping only shared columns."""
    float_cols = {}
    for g in GROUPS:
        float_cols[f"P_meta_{g}"] = pl.Float64
        float_cols[f"BETA_meta_{g}"] = pl.Float64
        float_cols[f"HetP_{g}"] = pl.Float64
    float_cols["gnomADv4.1_exome_nfe_frq"] = pl.Float64
    float_cols["gnomADv4.1_genome_nfe_frq"] = pl.Float64
    float_cols["PIP"] = pl.Utf8

    df = pl.read_csv(path, separator="\t", null_values=["", "NA", "."],
                     schema_overrides=float_cols)
    # keep only columns that exist and are in KEEP_COLS
    cols = [c for c in KEEP_COLS if c in df.columns]
    return df.select(cols)


def compute_group_stats(df: pl.DataFrame, group: str) -> pl.DataFrame:
    """Compute mlog10p, beta, se for a specific disease group."""
    p_col = f"P_meta_{group}"
    beta_col = f"BETA_meta_{group}"

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        # mlog10p
        df = df.with_columns(
            pl.when(pl.col(p_col).is_not_null() & (pl.col(p_col) > 0))
                .then(pl.max_horizontal((-np.log10(pl.col(p_col))).round(4), 0.0))
                .when(pl.col(p_col).is_not_null() & pl.col(p_col).eq(0))
                .then(324.0)
                .otherwise(None)
                .alias("mlog10p"),
        )

        # se = |beta| / |z| where z = qnorm(p/2)
        pvals = df[p_col].to_numpy()
        betas = df[beta_col].to_numpy()
        z = np.abs(norm.ppf(pvals / 2))
        z[z == 0] = np.nan
        se_arr = np.where(
            np.isfinite(betas) & np.isfinite(z) & (z > 0),
            np.abs(betas) / z, np.nan
        )
        df = df.with_columns(pl.Series("se_raw", se_arr, dtype=pl.Float64))

        # format beta/se as scientific notation
        df = df.with_columns(
            pl.col(beta_col).map_elements(
                lambda x: f"{x:.4e}" if x is not None and np.isfinite(x) else None,
                return_dtype=pl.Utf8,
            ).alias("beta"),
            pl.col("se_raw").map_elements(
                lambda x: f"{x:.4e}" if x is not None and np.isfinite(x) else None,
                return_dtype=pl.Utf8,
            ).alias("se"),
        )

        # gnomAD exome frequency as af_overall
        df = df.with_columns(
            pl.col("gnomADv4.1_exome_nfe_frq").map_elements(
                lambda x: f"{x:.4e}" if x is not None and np.isfinite(x) else None,
                return_dtype=pl.Utf8,
            ).alias("af_overall"),
        )

    return df


def build_output(df: pl.DataFrame, group: str, n_cases: int, n_controls: int) -> pl.DataFrame:
    """Build final output for one disease group."""
    trait = TRAIT_MAP[group]
    out = df.select(
        pl.lit("IBD_exome_2026").alias("#dataset"),
        "chr", "pos", "ref", "alt",
        pl.col("Gene_symbol").alias("gene"),
        pl.col("most_severe_consequence").alias("annotation"),
        "mlog10p",
        "beta",
        "se",
        "af_overall",
        pl.lit(None, dtype=pl.Utf8).alias("af_cases"),
        pl.lit(None, dtype=pl.Utf8).alias("af_controls"),
        pl.lit(None, dtype=pl.Utf8).alias("ac"),
        pl.lit(None, dtype=pl.Utf8).alias("an"),
        pl.lit(n_cases).alias("n_cases"),
        pl.lit(n_controls).alias("n_controls"),
        pl.lit(trait).alias("trait"),
        pl.lit(trait).alias("trait_original"),
    ).sort("chr", "pos", "ref", "alt")

    return out


def main():
    args = parse_args()

    output_dir = args.output_dir or str(Path(args.input).parent)

    print(f"Reading variant data from {args.input}...")
    df = read_variants(args.input)
    print(f"  {df.height} rows from primary input")

    for extra_path in args.input_extra:
        print(f"Reading extra variants from {extra_path}...")
        extra = read_variants(extra_path)
        print(f"  {extra.height} rows")
        df = pl.concat([df, extra], how="align")

    # deduplicate by chr:pos:ref:alt
    n_before = df.height
    df = df.unique(subset=["chr", "pos", "ref", "alt"], keep="first")
    n_dupes = n_before - df.height
    if n_dupes > 0:
        print(f"  removed {n_dupes} duplicate variants")
    print(f"  {df.height} total rows")

    n_cases_map = dict(zip(GROUPS, args.n_cases))

    for group in GROUPS:
        print(f"\nProcessing {group}...")
        group_df = compute_group_stats(df, group)
        out = build_output(group_df, group, n_cases_map[group], args.n_controls)
        print(f"  {out.height} output rows")
        write_exome_output(out, f"{output_dir}/IBD_exome_2026_{group}_variant_results.munged.tsv.gz",
                           tabix_args=["-s2", "-b3", "-e3"])

    print("\nDone.")


if __name__ == "__main__":
    main()
