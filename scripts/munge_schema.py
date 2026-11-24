#!/usr/bin/env python3
"""Munge SCHEMA gene burden results into genebass-compatible long format."""

import warnings

import numpy as np
import polars as pl

from exome_utils import run_gene_burden_munge


ANNOTATION_CLASSES = [
    {"annotation": "PTV",      "p_col": "P ca/co (Class 1)", "or_col": "OR (PTV)",     "or_lo": "OR (PTV) lower bound",     "or_hi": "OR (PTV) upper bound"},
    {"annotation": "Class_I",  "p_col": "P ca/co (Class 1)", "or_col": "OR (Class I)",  "or_lo": "OR (Class I) lower bound",  "or_hi": "OR (Class I) upper bound"},
    {"annotation": "Class_II", "p_col": "P ca/co (Class 2)", "or_col": "OR (Class II)", "or_lo": "OR (Class II) lower bound", "or_hi": "OR (Class II) upper bound"},
    {"annotation": "combined", "p_col": "P ca/co (comb)",    "or_col": None,            "or_lo": None,                        "or_hi": None},
    {"annotation": "de_novo",  "p_col": "P de novo",         "or_col": None,            "or_lo": None,                        "or_hi": None},
    {"annotation": "meta",     "p_col": "P meta",            "or_col": None,            "or_lo": None,                        "or_hi": None},
]

N_CASES = 24248
N_CONTROLS = 97322


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


OUTPUT = dict(
    dataset="SCHEMA",
    trait="Schizophrenia",
    trait_original="SCZ",
    total_variants=pl.lit(None, dtype=pl.Utf8),
    n_cases=pl.lit(N_CASES),
    n_controls=pl.lit(N_CONTROLS),
    flags=pl.lit("NA"),
)


def main():
    run_gene_burden_munge(
        __doc__, "SCHEMA", "/mnt/disks/data/gencode.v19.annotation.genes.tsv",
        read_schema, melt_to_long, len(ANNOTATION_CLASSES), compute_stats, OUTPUT, unit="genes",
    )


if __name__ == "__main__":
    main()
