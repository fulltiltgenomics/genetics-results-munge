#!/usr/bin/env python3
"""Munge IBD Supplementary Table ST7 (Burden test) into genebass-compatible gene burden format.

The input is a spreadsheet export where gene-level info (chr, Gene_symbol, Best_Associated_Trait,
Best_Associated_Mask, LD variants, LD correction) only appears on the first row of each gene block.
Continuation rows have empty values in those first 6 columns.

Outputs three files: IBD, CD, UC.
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import polars as pl

from sumstat_utils import write_exome_output

GROUPS = ["IBD", "CD", "UC"]

TRAIT_MAP = {
    "IBD": "inflammatory_bowel_disease",
    "CD": "crohns_disease",
    "UC": "ulcerative_colitis",
}

P_THRESHOLD = 1e-4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True,
                        help="Path to Supplementary_Tables_-_ST7_Burden_test.tsv")
    parser.add_argument("--gencode", default="/mnt/disks/data/gencode.v45.annotation.genes.tsv",
                        help="Path to gencode annotation genes TSV")
    parser.add_argument("--output-dir", help="Output directory (default: same as input)")
    parser.add_argument("--per-trait-dir", help="If given, also write one unfiltered <trait>.tsv.gz there (local or gs://)")
    parser.add_argument("--n-cases", required=True, nargs=3, type=int, metavar=("IBD", "CD", "UC"),
                        help="Number of cases for IBD, CD, UC")
    parser.add_argument("--n-controls", required=True, type=int, help="Number of controls (shared)")
    return parser.parse_args()


def read_gencode(path: str) -> pl.DataFrame:
    """Read gencode annotation, convert chr to numeric."""
    df = pl.read_csv(path, separator="\t")
    df = df.with_columns(
        pl.col("gene_id").str.split(".").list.first().alias("gene_id_base"),
        pl.col("chrom").cast(pl.Utf8)
            .str.replace("X", "23").str.replace("Y", "24").str.replace("M", "26")
            .cast(pl.Int32).alias("gene_chr"),
        pl.col("gene_start").cast(pl.Int32).alias("gene_start_pos"),
        pl.col("gene_end").cast(pl.Int32).alias("gene_end_pos"),
        pl.col("gene_name").alias("gene"),
    ).select("gene_id_base", "gene", "gene_chr", "gene_start_pos", "gene_end_pos")
    return df.unique(subset=["gene"], keep="first")


def read_burden(path: str) -> pl.DataFrame:
    """Read ST7 burden file, forward-fill gene info from block headers."""
    df = pl.read_csv(path, separator="\t", null_values=["", "NA"],
                     schema_overrides={
                         "P_meta": pl.Float64,
                         "beta_meta": pl.Float64,
                         "se_meta": pl.Float64,
                         "Het_P": pl.Float64,
                         "CAF_UKBB_WES_controls": pl.Float64,
                     })

    # forward-fill the gene block columns
    for col in ["chr", "Gene_symbol", "Best_Associated_Trait", "Best_Associated_Mask",
                "Variants in LD with best associated mask (R2 > 0.001)", "LD correction"]:
        df = df.with_columns(pl.col(col).forward_fill())

    return df


def compute_stats(df: pl.DataFrame) -> pl.DataFrame:
    """Compute mlog10p from p-value, format beta/se."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        df = df.with_columns(
            pl.when(pl.col("P_meta").is_not_null() & (pl.col("P_meta") > 0))
                .then(pl.max_horizontal((-np.log10(pl.col("P_meta"))).round(4), 0.0))
                .when(pl.col("P_meta").is_not_null() & pl.col("P_meta").eq(0))
                .then(324.0)
                .otherwise(None)
                .alias("mlog10p_burden"),
        )

        df = df.with_columns(
            pl.col("beta_meta").map_elements(
                lambda x: f"{x:.4e}" if x is not None else None, return_dtype=pl.Utf8
            ).alias("beta"),
            pl.col("se_meta").map_elements(
                lambda x: f"{x:.4e}" if x is not None else None, return_dtype=pl.Utf8
            ).alias("se"),
        )

    return df


def build_output(df: pl.DataFrame, gencode: pl.DataFrame, group: str,
                 n_cases: int, n_controls: int) -> pl.DataFrame:
    """Filter to group, keep genes where any mask has P_meta < threshold, join with gencode."""
    df = df.filter(pl.col("Trait") == group)

    # keep only genes where at least one mask has P_meta < threshold
    sig_genes = (
        df.filter(pl.col("P_meta") < P_THRESHOLD)
        .select("Gene_symbol").unique()
    )
    df = df.join(sig_genes, on="Gene_symbol", how="inner")

    df = df.with_columns(pl.col("Mask").alias("annotation"))

    # join with gencode on gene name
    n_before = df.height
    df = df.join(gencode, left_on="Gene_symbol", right_on="gene", how="inner")
    n_unmatched = n_before - df.height
    if n_unmatched > 0:
        print(f"  warning: {n_unmatched} rows not found in gencode")

    trait = TRAIT_MAP[group]
    out = df.select(
        pl.lit("IBD_exome_2026").alias("#dataset"),
        pl.lit(trait).alias("trait"),
        pl.col("Gene_symbol").alias("gene"),
        pl.col("gene_id_base").alias("gene_id"),
        "gene_chr",
        "gene_start_pos",
        "gene_end_pos",
        "annotation",
        "mlog10p_burden",

        "beta",
        "se",
        pl.lit(None, dtype=pl.Utf8).alias("total_variants"),
        pl.lit(None, dtype=pl.Utf8).alias("total_variants_pheno"),
        pl.lit(n_cases).alias("n_cases"),
        pl.lit(n_controls).alias("n_controls"),
        pl.lit(trait).alias("trait_original"),
        pl.lit("NA").alias("flags"),
    ).sort("gene_chr", "gene_start_pos", "gene_end_pos", "annotation")

    return out


def main():
    args = parse_args()

    output_dir = args.output_dir or str(Path(args.input).parent)

    print(f"Reading gencode from {args.gencode}...")
    gencode = read_gencode(args.gencode)
    print(f"  {gencode.height} genes")

    print(f"Reading burden data from {args.input}...")
    df = read_burden(args.input)
    print(f"  {df.height} total rows")

    print("Computing stats...")
    df = compute_stats(df)

    n_cases_map = dict(zip(GROUPS, args.n_cases))

    for group in GROUPS:
        print(f"\nProcessing {group}...")
        out = build_output(df, gencode, group, n_cases_map[group], args.n_controls)
        print(f"  {out.height} output rows")
        write_exome_output(out, f"{output_dir}/IBD_exome_2026_{group}_gene_results.munged.tsv.gz",
                           tabix_args=["-s5", "-b6", "-e6"], mlog10p_col="mlog10p_burden",
                           per_trait_dir=args.per_trait_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
