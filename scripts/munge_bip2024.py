#!/usr/bin/env python3
"""Munge BIP 2024 multi-ancestry GWAS summary statistics and create AF-AF plot against gnomAD."""

import warnings

import matplotlib
matplotlib.use("Agg")
import numpy as np
import polars as pl
from scipy.special import log_ndtr

from sumstat_utils import chrom_code_expr, run_rsid_munge


HRC_FREQ_COLS = ["HRC_AMR_FRQ_A1", "HRC_AFR_FRQ_A1", "HRC_EAS_FRQ_A1", "HRC_EUR_FRQ_A1"]
HRC_RENAME = {
    "HRC_AMR_FRQ_A1": "af_amr",
    "HRC_AFR_FRQ_A1": "af_afr",
    "HRC_EAS_FRQ_A1": "af_eas",
    "HRC_EUR_FRQ_A1": "af_eur",
}


def read_sumstats(filepath: str) -> pl.DataFrame:
    """Read BIP 2024 multi-ancestry file, compute beta/mlog10p, sort."""
    df = pl.read_csv(
        filepath,
        separator=" ",
        null_values=["NA"],
        truncate_ragged_lines=True,
        schema_overrides={
            "CHR": pl.Utf8,
            "BP": pl.Int32,
            "OR": pl.Float64,
            "SE": pl.Float64,
            "P": pl.Float64,
            "INFO": pl.Float64,
            **{c: pl.Float64 for c in HRC_FREQ_COLS},
        },
        ignore_errors=True,
    ).rename({
        "CHR": "chr",
        "BP": "pos",
        "SNP": "rsid",
        "Neff_half": "neff",
        "Nca": "n_cases",
        "Nco": "n_controls",
        **HRC_RENAME,
    })

    # a contig with no numeric code casts to null and is dropped
    df = df.with_columns(
        chrom_code_expr("chr").cast(pl.Int32, strict=False).alias("chr"),
    )
    df = df.filter(pl.col("chr").is_not_null())

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        df = df.with_columns(
            np.log(pl.col("OR")).round(6).alias("beta"),
            pl.col("af_eur").alias("af"),
        ).with_columns(
            pl.when(pl.col("P").eq(0))
            .then(
                ((-log_ndtr(-(pl.col("beta") / pl.col("SE")).abs()) - np.log(2)) / np.log(10)).round(4)
            )
            .otherwise((-np.log10(pl.col("P"))).round(4))
            .alias("mlog10p"),
        )

    return df


def main():
    run_rsid_munge(
        __doc__, "Path to BIP 2024 .gz file", read_sumstats, list(HRC_RENAME.values()), "BIP2024 alt AF (EUR, after flip)",
    )


if __name__ == "__main__":
    main()
