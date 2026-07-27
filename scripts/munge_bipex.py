#!/usr/bin/env python3
"""Munge BipEx2 gene burden results into genebass-compatible long format."""

import argparse
import warnings
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import norm

from sumstat_utils import write_exome_output

ANNOTATION_CLASSES = [
    {"annotation": "missense",        "p_col": "mis_p_value",     "or_col": "mis_odds_ratio",     "case_col": "mis_case_carrier",     "ctrl_col": "mis_ctrl_carrier"},
    {"annotation": "pLoF",            "p_col": "ptv_p_value",     "or_col": "ptv_odds_ratio",     "case_col": "ptv_case_carrier",     "ctrl_col": "ptv_ctrl_carrier"},
    {"annotation": "pLoF|missense",   "p_col": "ptv_mis_p_value", "or_col": "ptv_mis_odds_ratio", "case_col": "ptv_mis_case_carrier", "ctrl_col": "ptv_mis_ctrl_carrier"},
]

GENE_COLS = ["gene_id", "gene", "gene_chr", "gene_start_pos", "gene_end_pos"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to BipEx2 gene results .bgz file")
    parser.add_argument("--gencode", default="/mnt/disks/data/gencode.v39.annotation.genes.tsv",
                        help="Path to gencode annotation genes TSV")
    parser.add_argument("--output-dir", help="Output directory (default: same as input)")
    return parser.parse_args()


def read_gencode(path: str) -> pl.DataFrame:
    """Read gencode annotation, strip version from gene_id, convert chr to numeric."""
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
    df = df.unique(subset=["gene_id_base"], keep="first")
    return df


def read_bipex(path: str) -> pl.DataFrame:
    """Read BipEx2 gene burden results."""
    float_cols = {col for cls in ANNOTATION_CLASSES
                  for col in [cls["p_col"], cls["or_col"]]}.union({"n_cases", "n_controls"})
    carrier_cols = {col for cls in ANNOTATION_CLASSES
                    for col in [cls["case_col"], cls["ctrl_col"]]}
    return pl.read_csv(
        path, separator="\t", null_values=["NA", ""],
        schema_overrides={**{c: pl.Float64 for c in float_cols},
                          **{c: pl.Int64 for c in carrier_cols}},
    )


def melt_to_long(df: pl.DataFrame) -> pl.DataFrame:
    """Melt wide BipEx data to long format with one row per gene per annotation class."""
    frames = []
    for cls in ANNOTATION_CLASSES:
        sub = df.select(
            "gene_id", "flags", "n_cases", "n_controls",
            pl.col(cls["p_col"]).alias("pvalue"),
            pl.col(cls["or_col"]).alias("OR"),
            pl.col(cls["case_col"]).alias("case_carrier"),
            pl.col(cls["ctrl_col"]).alias("ctrl_carrier"),
        ).with_columns(pl.lit(cls["annotation"]).alias("annotation"))
        sub = sub.filter(pl.col("pvalue").is_not_null())
        frames.append(sub)
    return pl.concat(frames)


def compute_stats(df: pl.DataFrame) -> pl.DataFrame:
    """Compute beta, se, and mlog10p from OR and p-value.

    se = |log(OR)| / |z|  where z = qnorm(p/2)
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        df = df.with_columns(
            # beta = log(OR) when OR is finite and > 0
            pl.when(pl.col("OR").is_not_null() & pl.col("OR").is_finite() & (pl.col("OR") > 0))
                .then(np.log(pl.col("OR")))
                .otherwise(None)
                .alias("beta"),
        )

        # se = |beta| / |z| where z = qnorm(p/2)
        # use numpy/scipy on extracted arrays for the z computation
        pvals = df["pvalue"].to_numpy()
        betas = df["beta"].to_numpy()
        z = np.abs(norm.ppf(pvals / 2))
        # avoid division by zero
        z[z == 0] = np.nan
        se = np.where(np.isfinite(betas) & np.isfinite(z), np.abs(betas) / z, np.nan)
        df = df.with_columns(
            pl.Series("se", se, dtype=pl.Float64)
                .map_elements(lambda x: None if (x is None or np.isnan(x)) else x, return_dtype=pl.Float64),
        )

        # mlog10p from p-value
        df = df.with_columns(
            pl.when(pl.col("pvalue").is_not_null() & (pl.col("pvalue") > 0))
                .then(pl.max_horizontal((-np.log10(pl.col("pvalue"))).round(4), 0.0))
                .when(pl.col("pvalue").is_not_null() & pl.col("pvalue").eq(0))
                .then(324.0)
                .otherwise(None)
                .alias("mlog10p_burden"),
        )

        # filter out rows where beta is NA (e.g. OR=Infinity or 0)
        df = df.filter(pl.col("beta").is_not_null())

        # total carriers (cases + controls)
        df = df.with_columns(
            (pl.col("case_carrier") + pl.col("ctrl_carrier")).alias("total_variants"),
            (pl.col("case_carrier") + pl.col("ctrl_carrier")).alias("total_variants_pheno"),
        )

        # format beta/se as scientific notation
        df = df.with_columns(
            pl.col("beta").map_elements(lambda x: f"{x:.4e}", return_dtype=pl.Utf8).alias("beta"),
            pl.col("se").map_elements(lambda x: f"{x:.4e}", return_dtype=pl.Utf8).alias("se"),
        )

    return df


def build_output(df: pl.DataFrame, gencode: pl.DataFrame) -> pl.DataFrame:
    """Join with gencode and build final output DataFrame."""
    df = df.join(gencode, left_on="gene_id", right_on="gene_id_base", how="inner")

    n_unmatched = df.filter(pl.col("gene_chr").is_null()).height
    if n_unmatched > 0:
        print(f"  warning: {n_unmatched} genes not found in gencode")

    out = df.select(
        pl.lit("BipEx2").alias("#dataset"),
        pl.lit("bipolar_disorder").alias("trait"),
        "gene",
        "gene_id",
        "gene_chr",
        "gene_start_pos",
        "gene_end_pos",
        "annotation",
        "mlog10p_burden",

        "beta",
        "se",
        "total_variants",
        "total_variants_pheno",
        pl.col("n_cases").cast(pl.Int64),
        pl.col("n_controls").cast(pl.Int64),
        pl.lit("bipolar_disorder").alias("trait_original"),
        pl.col("flags").fill_null("NA").alias("flags"),
    ).sort("gene_chr", "gene_start_pos", "gene_end_pos", "annotation")

    return out


def main():
    args = parse_args()

    output_dir = args.output_dir or str(Path(args.input).parent)

    print(f"Reading gencode from {args.gencode}...")
    gencode = read_gencode(args.gencode)
    print(f"  {gencode.height} genes")

    print(f"Reading BipEx2 from {args.input}...")
    bipex = read_bipex(args.input)
    print(f"  {bipex.height} rows")

    print("Melting to long format...")
    long = melt_to_long(bipex)
    print(f"  {long.height} rows across {len(ANNOTATION_CLASSES)} annotation classes")

    print("Computing stats...")
    long = compute_stats(long)

    print("Building output...")
    out = build_output(long, gencode)
    print(f"  {out.height} total rows")

    print("Writing results...")
    write_exome_output(out, f"{output_dir}/BipEx2_gene_results.munged.tsv.gz",
                       tabix_args=["-s5", "-b6", "-e6"], mlog10p_col="mlog10p_burden")

    print("Done.")


if __name__ == "__main__":
    main()
