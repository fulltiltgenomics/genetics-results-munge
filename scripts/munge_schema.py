#!/usr/bin/env python3
"""Munge SCHEMA gene burden results into genebass-compatible long format."""

import argparse
import warnings
from pathlib import Path

import numpy as np
import polars as pl

from sumstat_utils import write_exome_output


ANNOTATION_CLASSES = [
    {"annotation": "PTV",      "p_col": "P ca/co (Class 1)", "or_col": "OR (PTV)",     "or_lo": "OR (PTV) lower bound",     "or_hi": "OR (PTV) upper bound"},
    {"annotation": "Class_I",  "p_col": "P ca/co (Class 1)", "or_col": "OR (Class I)",  "or_lo": "OR (Class I) lower bound",  "or_hi": "OR (Class I) upper bound"},
    {"annotation": "Class_II", "p_col": "P ca/co (Class 2)", "or_col": "OR (Class II)", "or_lo": "OR (Class II) lower bound", "or_hi": "OR (Class II) upper bound"},
    {"annotation": "combined", "p_col": "P ca/co (comb)",    "or_col": None,            "or_lo": None,                        "or_hi": None},
    {"annotation": "de_novo",  "p_col": "P de novo",         "or_col": None,            "or_lo": None,                        "or_hi": None},
    {"annotation": "meta",     "p_col": "P meta",            "or_col": None,            "or_lo": None,                        "or_hi": None},
]

GENE_COLS = ["gene_id", "gene", "gene_chr", "gene_start_pos", "gene_end_pos"]

N_CASES = 24248
N_CONTROLS = 97322


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to SCHEMA gene results .bgz file")
    parser.add_argument("--gencode", default="/mnt/disks/data/gencode.v19.annotation.genes.tsv",
                        help="Path to gencode annotation genes TSV")
    parser.add_argument("--output-dir", help="Output directory (default: same as input)")
    parser.add_argument("--per-trait-dir", help="If given, also write one unfiltered <trait>.tsv.gz there (local or gs://)")
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
    # deduplicate on gene_id_base, keep first occurrence
    df = df.unique(subset=["gene_id_base"], keep="first")
    return df


def read_schema(path: str) -> pl.DataFrame:
    """Read SCHEMA gene burden results."""
    # columns with many NA values get inferred as String; force Float64
    float_cols = {col for cls in ANNOTATION_CLASSES
                  for col in [cls["p_col"], cls["or_col"], cls["or_lo"], cls["or_hi"]]
                  if col is not None}
    return pl.read_csv(
        path, separator="\t", null_values=["NA"],
        schema_overrides={c: pl.Float64 for c in float_cols},
    )


def melt_to_long(schema: pl.DataFrame) -> pl.DataFrame:
    """Melt wide SCHEMA data to long format with one row per gene per annotation class."""
    frames = []
    for cls in ANNOTATION_CLASSES:
        cols_to_select = ["gene_id"]
        renames = {cls["p_col"]: "pvalue"}
        cols_to_select.append(cls["p_col"])

        if cls["or_col"]:
            cols_to_select.extend([cls["or_col"], cls["or_lo"], cls["or_hi"]])
            renames[cls["or_col"]] = "OR"
            renames[cls["or_lo"]] = "OR_lo"
            renames[cls["or_hi"]] = "OR_hi"

        sub = schema.select(cols_to_select).rename(renames)

        if not cls["or_col"]:
            sub = sub.with_columns(
                pl.lit(None, dtype=pl.Float64).alias("OR"),
                pl.lit(None, dtype=pl.Float64).alias("OR_lo"),
                pl.lit(None, dtype=pl.Float64).alias("OR_hi"),
            )

        sub = sub.with_columns(pl.lit(cls["annotation"]).alias("annotation"))
        # drop rows where p-value is null
        sub = sub.filter(pl.col("pvalue").is_not_null())
        frames.append(sub)

    return pl.concat(frames)


def compute_stats(df: pl.DataFrame) -> pl.DataFrame:
    """Compute beta, se, and mlog10p from OR/CI and p-values."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        df = df.with_columns(
            # beta = log(OR) when OR is finite and > 0
            pl.when(pl.col("OR").is_not_null() & pl.col("OR").is_finite() & (pl.col("OR") > 0))
                .then(np.log(pl.col("OR")))
                .otherwise(None)
                .alias("beta"),
            # se from CI: (log(upper) - log(lower)) / (2 * 1.96)
            pl.when(
                pl.col("OR_lo").is_not_null() & pl.col("OR_hi").is_not_null()
                & pl.col("OR_lo").is_finite() & pl.col("OR_hi").is_finite()
                & (pl.col("OR_lo") > 0) & (pl.col("OR_hi") > 0)
            )
                .then((np.log(pl.col("OR_hi")) - np.log(pl.col("OR_lo"))) / (2 * 1.96))
                .otherwise(None)
                .alias("se"),
        )

        # mlog10p from raw p-value
        df = df.with_columns(
            pl.when(pl.col("pvalue").is_not_null() & (pl.col("pvalue") > 0))
                .then(pl.max_horizontal((-np.log10(pl.col("pvalue"))).round(4), 0.0))
                .when(pl.col("pvalue").is_not_null() & (pl.col("pvalue").eq(0)))
                .then(324.0)
                .otherwise(None)
                .alias("mlog10p_burden"),
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
        pl.lit("SCHEMA").alias("#dataset"),
        pl.lit("schizophrenia").alias("trait"),
        "gene",
        "gene_id",
        "gene_chr",
        "gene_start_pos",
        "gene_end_pos",
        "annotation",
        "mlog10p_burden",

        "beta",
        "se",
        pl.lit(None, dtype=pl.Utf8).alias("total_variants"),
        pl.lit(None, dtype=pl.Utf8).alias("total_variants_pheno"),
        pl.lit(N_CASES).alias("n_cases"),
        pl.lit(N_CONTROLS).alias("n_controls"),
        pl.lit("schizophrenia").alias("trait_original"),
        pl.lit("NA").alias("flags"),
    ).sort("gene_chr", "gene_start_pos", "gene_end_pos", "annotation")

    return out


def main():
    args = parse_args()

    output_dir = args.output_dir or str(Path(args.input).parent)

    print(f"Reading gencode from {args.gencode}...")
    gencode = read_gencode(args.gencode)
    print(f"  {gencode.height} genes")

    print(f"Reading SCHEMA from {args.input}...")
    schema = read_schema(args.input)
    print(f"  {schema.height} genes")

    print("Melting to long format...")
    long = melt_to_long(schema)
    print(f"  {long.height} rows across {len(ANNOTATION_CLASSES)} annotation classes")

    print("Computing stats...")
    long = compute_stats(long)

    print("Building output...")
    out = build_output(long, gencode)
    print(f"  {out.height} total rows")

    print("Writing results...")
    write_exome_output(out, f"{output_dir}/SCHEMA_gene_results.munged.tsv.gz",
                       tabix_args=["-s5", "-b6", "-e6"], mlog10p_col="mlog10p_burden",
                       per_trait_dir=args.per_trait_dir)

    print("Done.")


if __name__ == "__main__":
    main()
