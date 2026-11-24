#!/usr/bin/env python3
"""Munge PGC GWAS summary statistics and create AF-AF plot against gnomAD."""

import warnings

import matplotlib
matplotlib.use("Agg")
import numpy as np
import polars as pl
from scipy.special import log_ndtr

from sumstat_utils import chrom_code_expr, run_rsid_munge

# case and control frequencies of A1: flipped with the effect, written between af and INFO
FREQ_COLS = ["FRQ_A", "FRQ_U"]


def read_pgc_sumstats(filepath: str) -> pl.DataFrame:
    """Read PGC daner file, compute beta/mlog10p, sort.
    Keeps beta as Float64 for later flipping."""
    # detect variable column names from header
    cols = pl.read_csv(filepath, separator="\t", n_rows=0).columns
    frq_a_col = next(c for c in cols if c.startswith("FRQ_A_"))
    frq_u_col = next(c for c in cols if c.startswith("FRQ_U_"))
    neff_col = "Neff_half" if "Neff_half" in cols else "Neff"

    df = pl.read_csv(
        filepath,
        separator="\t",
        null_values=["NA"],
        truncate_ragged_lines=True,
        schema_overrides={
            "CHR": pl.Utf8,
            "BP": pl.Int32,
            "OR": pl.Float64,
            "SE": pl.Float64,
            "P": pl.Float64,
            frq_a_col: pl.Float64,
            frq_u_col: pl.Float64,
            "INFO": pl.Float64,
        },
        ignore_errors=True,
    ).rename({
        frq_a_col: "FRQ_A",
        frq_u_col: "FRQ_U",
        neff_col: "neff",
        "CHR": "chr",
        "BP": "pos",
        "SNP": "rsid",
        "Nca": "n_cases",
        "Nco": "n_controls",
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
        ).with_columns(
            # mlog10p with underflow handling
            pl.when(pl.col("P").eq(0))
            .then(
                ((-log_ndtr(-(pl.col("beta") / pl.col("SE")).abs()) - np.log(2)) / np.log(10)).round(4)
            )
            .otherwise((-np.log10(pl.col("P"))).round(4))
            .alias("mlog10p"),
            # overall af weighted by case/control counts
            ((pl.col("FRQ_A") * pl.col("n_cases") + pl.col("FRQ_U") * pl.col("n_controls"))
             / (pl.col("n_cases") + pl.col("n_controls"))).round(6).alias("af"),
        )

    return df


def main():
    run_rsid_munge(__doc__, "Path to PGC daner .gz file", read_pgc_sumstats, FREQ_COLS, "PGC alt AF (after flip)")


if __name__ == "__main__":
    main()
