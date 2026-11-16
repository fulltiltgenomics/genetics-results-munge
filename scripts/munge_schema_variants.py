#!/usr/bin/env python3
"""Munge SCHEMA variant results into genebass-compatible format."""

import warnings

import numpy as np
import polars as pl
from scipy.special import log_ndtr

from exome_utils import run_schema_variant_munge
from sumstat_utils import chrom_code_expr


CODING_CONSEQUENCES = {
    "missense_variant_mpc_<2", "missense_variant_mpc_2-3", "missense_variant_mpc_>=3",
    "frameshift_variant", "inframe_deletion", "inframe_insertion",
    "stop_gained", "stop_lost", "stop_retained_variant", "start_lost",
    "splice_donor_variant", "splice_acceptor_variant",
    "coding_sequence_variant", "protein_altering_variant",
}

MLOG10P_THRESHOLD = 4


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
        chrom_code_expr(pl.col("locus").str.split(":").list.first()).cast(pl.Int32).alias("chr"),
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


def main():
    run_schema_variant_munge(
        __doc__, "SCHEMA", "/mnt/disks/data/gencode.v19.annotation.genes.tsv",
        read_schema_variants, compute_stats, beta_col="est",
    )


if __name__ == "__main__":
    main()
