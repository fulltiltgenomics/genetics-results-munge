#!/usr/bin/env python3
"""Munge SCHEMA2 gene burden results into genebass-compatible long format."""

import polars as pl

from exome_utils import burden_stats_from_or_p, run_gene_burden_munge

ANNOTATION_CLASSES = [
    {"annotation": "PTV",          "p_col": "ptv_p_value",     "or_col": "ptv_odds_ratio"},
    {"annotation": "PTV_missense", "p_col": "ptv_mis_p_value", "or_col": "ptv_mis_odds_ratio"},
]


def read_schema2(path: str) -> pl.DataFrame:
    """Read SCHEMA2 gene burden results."""
    float_cols = {col for cls in ANNOTATION_CLASSES
                  for col in [cls["p_col"], cls["or_col"]]}.union({"n_cases", "n_controls"})
    return pl.read_csv(
        path, separator="\t", null_values=["NA", ""],
        schema_overrides={c: pl.Float64 for c in float_cols},
    )


def melt_to_long(df: pl.DataFrame) -> pl.DataFrame:
    """Melt wide data to long format with one row per gene per annotation class."""
    frames = []
    for cls in ANNOTATION_CLASSES:
        sub = df.select(
            "gene_id", "n_cases", "n_controls",
            pl.col(cls["p_col"]).alias("pvalue"),
            pl.col(cls["or_col"]).alias("OR"),
        ).with_columns(pl.lit(cls["annotation"]).alias("annotation"))
        sub = sub.filter(pl.col("pvalue").is_not_null())
        frames.append(sub)
    return pl.concat(frames)


OUTPUT = dict(
    dataset="SCHEMA2",
    trait="Schizophrenia",
    trait_original="SCZ",
    total_variants=pl.lit(None, dtype=pl.Utf8),
    n_cases=pl.col("n_cases").cast(pl.Int64),
    n_controls=pl.col("n_controls").cast(pl.Int64),
    flags=pl.lit("NA"),
)


def main():
    run_gene_burden_munge(
        __doc__, "SCHEMA2", "/mnt/disks/data/gencode.v39.annotation.genes.tsv",
        read_schema2, melt_to_long, len(ANNOTATION_CLASSES), burden_stats_from_or_p, OUTPUT,
    )


if __name__ == "__main__":
    main()
