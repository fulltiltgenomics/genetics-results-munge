#!/usr/bin/env python3
"""Munge BipEx2 gene burden results into genebass-compatible long format."""

import polars as pl

from exome_utils import burden_stats_from_or_p, run_gene_burden_munge

ANNOTATION_CLASSES = [
    {"annotation": "missense",        "p_col": "mis_p_value",     "or_col": "mis_odds_ratio",     "case_col": "mis_case_carrier",     "ctrl_col": "mis_ctrl_carrier"},
    {"annotation": "pLoF",            "p_col": "ptv_p_value",     "or_col": "ptv_odds_ratio",     "case_col": "ptv_case_carrier",     "ctrl_col": "ptv_ctrl_carrier"},
    {"annotation": "pLoF|missense",   "p_col": "ptv_mis_p_value", "or_col": "ptv_mis_odds_ratio", "case_col": "ptv_mis_case_carrier", "ctrl_col": "ptv_mis_ctrl_carrier"},
]


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


OUTPUT = dict(
    dataset="BipEx2",
    trait="Bipolar disorder",
    trait_original="BIP",
    # total carriers (cases + controls)
    total_variants=pl.col("case_carrier") + pl.col("ctrl_carrier"),
    n_cases=pl.col("n_cases").cast(pl.Int64),
    n_controls=pl.col("n_controls").cast(pl.Int64),
    flags=pl.col("flags").fill_null("NA"),
)


def main():
    run_gene_burden_munge(
        __doc__, "BipEx2", "/mnt/disks/data/gencode.v39.annotation.genes.tsv",
        read_bipex, melt_to_long, len(ANNOTATION_CLASSES), burden_stats_from_or_p, OUTPUT,
    )


if __name__ == "__main__":
    main()
