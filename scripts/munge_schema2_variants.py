#!/usr/bin/env python3
"""Munge SCHEMA2 variant results into genebass-compatible format.

SCHEMA2 provides allele counts but no p-values or effect sizes, so we compute
log-odds ratios with Haldane correction and derive p-values from z-scores.
"""

import warnings

import numpy as np
import polars as pl
from scipy.special import log_ndtr

from exome_utils import run_schema_variant_munge
from sumstat_utils import chrom_code_expr

CODING_CONSEQUENCES = {
    "missense_variant", "frameshift_variant",
    "inframe_deletion", "inframe_insertion",
    "stop_gained", "stop_lost", "stop_retained_variant", "start_lost",
    "splice_donor_variant", "splice_acceptor_variant",
    "coding_sequence_variant", "protein_altering_variant",
}

MLOG10P_THRESHOLD = 4


def read_schema2_variants(path: str) -> pl.DataFrame:
    """Read SCHEMA2 variant results, filter to coding variants."""
    print("  scanning file...")
    df = (
        pl.scan_csv(
            path, separator="\t", null_values=["NA"],
            schema_overrides={
                "ac_case": pl.Int64, "an_case": pl.Int64,
                "ac_ctrl": pl.Int64, "an_ctrl": pl.Int64,
                "n_de_novo": pl.Int64, "mpc": pl.Float64,
            },
        )
        .filter(pl.col("consequence").is_in(CODING_CONSEQUENCES))
        .collect()
    )
    print(f"  {df.height} coding rows")

    # parse locus "chrN:POS" -> chr, pos
    df = df.with_columns(
        chrom_code_expr(pl.col("locus").str.split(":").list.first()).cast(pl.Int32).alias("chr"),
        pl.col("locus").str.split(":").list.last().cast(pl.Int64).alias("pos"),
    )

    # parse alleles '["REF","ALT"]' -> ref, alt
    df = df.with_columns(
        pl.col("alleles").str.strip_chars("[]").str.replace_all('"', '').str.split(","),
    ).with_columns(
        pl.col("alleles").list.first().alias("ref"),
        pl.col("alleles").list.last().alias("alt"),
    )

    return df


def compute_stats(df: pl.DataFrame) -> pl.DataFrame:
    """Compute log-OR, SE, mlog10p, and allele frequencies from counts.

    Uses Haldane correction (add 0.5) to handle zero cells in the 2x2 table.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        ac_case = df["ac_case"].to_numpy(zero_copy_only=False).astype(np.float64)
        an_case = df["an_case"].to_numpy(zero_copy_only=False).astype(np.float64)
        ac_ctrl = df["ac_ctrl"].to_numpy(zero_copy_only=False).astype(np.float64)
        an_ctrl = df["an_ctrl"].to_numpy(zero_copy_only=False).astype(np.float64)

        # Haldane-corrected cells
        a = ac_case + 0.5
        b = (an_case - ac_case) + 0.5
        c = ac_ctrl + 0.5
        d = (an_ctrl - ac_ctrl) + 0.5

        beta = np.log(a * d / (b * c))
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

        # mlog10p from z-score: -log10(2 * Phi(-|z|))
        z_abs = np.abs(beta / se)
        # log_ndtr gives log(Phi(x)), so log(p) = log_ndtr(-|z|) + log(2)
        log10p = (log_ndtr(-z_abs) + np.log(2)) / np.log(10)
        mlog10p = np.round(np.maximum(-log10p, 0.0), 4)

        # mask rows where an is 0
        valid = (an_case > 0) & (an_ctrl > 0)
        beta = np.where(valid, beta, np.nan)
        se = np.where(valid, se, np.nan)
        mlog10p = np.where(valid, mlog10p, np.nan)

        # allele frequencies
        af_overall = np.where(
            (an_case + an_ctrl) > 0,
            (ac_case + ac_ctrl) / (an_case + an_ctrl),
            np.nan,
        )
        af_cases = np.where(an_case > 0, ac_case / an_case, np.nan)
        af_controls = np.where(an_ctrl > 0, ac_ctrl / an_ctrl, np.nan)

        df = df.with_columns(
            pl.Series("beta", beta, dtype=pl.Float64),
            pl.Series("se", se, dtype=pl.Float64),
            pl.Series("mlog10p", mlog10p, dtype=pl.Float64),
            pl.Series("af_overall", af_overall, dtype=pl.Float64),
            pl.Series("af_cases", af_cases, dtype=pl.Float64),
            pl.Series("af_controls", af_controls, dtype=pl.Float64),
            (pl.col("ac_case") + pl.col("ac_ctrl")).alias("ac"),
            (pl.col("an_case") + pl.col("an_ctrl")).alias("an"),
        )

        # format numeric columns as scientific notation
        for col in ["beta", "se", "af_overall", "af_cases", "af_controls"]:
            df = df.with_columns(
                pl.col(col).map_elements(lambda x: f"{x:.4e}", return_dtype=pl.Utf8).alias(col),
            )

    return df


def main():
    run_schema_variant_munge(
        __doc__, "SCHEMA2", "/mnt/disks/data/gencode.v39.annotation.genes.tsv",
        read_schema2_variants, compute_stats, beta_col="beta",
    )


if __name__ == "__main__":
    main()
