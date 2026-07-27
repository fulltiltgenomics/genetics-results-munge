#!/usr/bin/env python3
"""Munge PGC GWAS summary statistics and create AF-AF plot against gnomAD."""

import argparse
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from scipy.special import log_ndtr

from sumstat_utils import (
    GNOMAD_AF_COLS, GNOMAD_DEFAULT,
    build_rsid_set, read_gnomad_filtered, stream_gnomad_by_rsid,
    write_sumstat_output,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to PGC daner .gz file")
    parser.add_argument("--gnomad", default=GNOMAD_DEFAULT, help="GCS path to gnomAD bgz file")
    parser.add_argument("--af-col", default="AF", choices=GNOMAD_AF_COLS, help="gnomAD AF column for plot (default: AF)")
    parser.add_argument("--output", help="Output munged TSV path (default: derived from input)")
    parser.add_argument("--gnomad-af-plot", action="store_true", help="Stream gnomAD and create AF-AF plot")
    parser.add_argument("--gnomad-filtered", help="Path to previously saved filtered gnomAD TSV (skip streaming)")
    parser.add_argument("--gnomad-source", default="g", choices=["g", "e"], help="Use gnomAD genomes (g) or exomes (e) (default: g)")
    parser.add_argument("--plot", help="Output AF-AF plot PNG path (default: derived from input)")
    return parser.parse_args()


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

    # map X→23 and convert chr to int
    df = df.with_columns(
        pl.col("chr").str.replace("X", "23").cast(pl.Int32, strict=False),
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


def join_with_gnomad(
    pgc: pl.DataFrame, gnomad: pl.DataFrame, af_col: str, gnomad_source: str = "g",
) -> pl.DataFrame:
    """Join PGC with gnomAD on rsid, classify allele orientation, flip to alt-effect."""
    source_label = "genomes" if gnomad_source == "g" else "exomes"
    gnomad_dedup = gnomad.filter(pl.col("genome_or_exome") == gnomad_source)
    print(f"  using gnomAD {source_label}: {gnomad_dedup.height} rows (from {gnomad.height} total)")

    # explode gnomAD rsids (comma-separated) into individual rows
    gnomad_exploded = gnomad_dedup.with_columns(
        pl.col("rsids").str.split(",").alias("rsid_list"),
    ).explode("rsid_list").rename({"rsid_list": "rsid"})

    # join on rsid
    joined = pgc.join(
        gnomad_exploded,
        on="rsid",
        how="inner",
        suffix="_gnomad",
    )

    # classify allele orientation (case-insensitive)
    joined = joined.with_columns(
        pl.when(
            (pl.col("A1").str.to_uppercase() == pl.col("ref").str.to_uppercase())
            & (pl.col("A2").str.to_uppercase() == pl.col("alt").str.to_uppercase())
        ).then(pl.lit("A1_is_ref"))
        .when(
            (pl.col("A1").str.to_uppercase() == pl.col("alt").str.to_uppercase())
            & (pl.col("A2").str.to_uppercase() == pl.col("ref").str.to_uppercase())
        ).then(pl.lit("A1_is_alt"))
        .otherwise(pl.lit("mismatch"))
        .alias("allele_match")
    )

    n_a1_ref = joined.filter(pl.col("allele_match") == "A1_is_ref").height
    n_a1_alt = joined.filter(pl.col("allele_match") == "A1_is_alt").height
    n_mismatch = joined.filter(pl.col("allele_match") == "mismatch").height
    matched_rsids = joined["rsid"].unique().len()
    n_unmatched = pgc["rsid"].unique().len() - matched_rsids
    print(f"  A1=ref (flip): {n_a1_ref}, A1=alt (no flip): {n_a1_alt}, mismatch: {n_mismatch}, no gnomAD match: {n_unmatched}")

    # keep A1_is_ref and A1_is_alt, drop mismatches
    result = joined.filter(pl.col("allele_match") != "mismatch")

    # use gnomAD chr/pos (build 38) for all rows
    result = result.with_columns(
        pl.col("chr_gnomad").alias("chr"),
        pl.col("pos_gnomad").alias("pos"),
    )

    # A1_is_ref: A1 is ref, A2 is alt — flip beta and AFs so effect is on alt
    # A1_is_alt: A1 is alt, A2 is ref — beta/AFs already represent alt effect
    is_flip = pl.col("allele_match") == "A1_is_ref"
    result = result.with_columns(
        pl.when(is_flip).then(-pl.col("beta")).otherwise(pl.col("beta")).alias("beta"),
        pl.when(is_flip).then((1.0 - pl.col("af")).round(6)).otherwise(pl.col("af")).alias("af"),
        pl.when(is_flip).then((1.0 - pl.col("FRQ_A")).round(6)).otherwise(pl.col("FRQ_A")).alias("FRQ_A"),
        pl.when(is_flip).then((1.0 - pl.col("FRQ_U")).round(6)).otherwise(pl.col("FRQ_U")).alias("FRQ_U"),
    )

    return result


def create_af_af_plot(
    df: pl.DataFrame,
    af_col: str,
    output_path: str,
    input_name: str,
) -> None:
    """Create AF-AF scatter plot: sumstat alt AF vs gnomAD alt AF."""
    sample_n = min(200_000, df.height)
    sample = df.sample(n=sample_n, seed=42)

    fig, ax = plt.subplots(figsize=(8, 8))

    ax.scatter(
        sample["af"].to_numpy(),
        sample[af_col].to_numpy(),
        alpha=0.3, s=1, c="blue",
        rasterized=True,
    )

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=0.5)
    ax.set_xlabel("PGC alt AF (after flip)")
    ax.set_ylabel(f"gnomAD {af_col}")
    ax.set_title(f"AF-AF: {input_name} (n={df.height:,})")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  plot saved to {output_path}")


def _prepare_output(df: pl.DataFrame) -> pl.DataFrame:
    """Format and select output columns."""
    return df.with_columns(
        pl.col("beta").map_elements(lambda x: f"{x:.3e}" if x is not None else None, return_dtype=pl.Utf8).alias("beta"),
        pl.col("SE").map_elements(lambda x: f"{x:.3e}" if x is not None else None, return_dtype=pl.Utf8).alias("se"),
    ).select(
        "chr", "pos", "ref", "alt", "rsid",
        "beta", "se", "mlog10p",
        "af", "FRQ_A", "FRQ_U", "INFO",
        "neff", "n_cases", "n_controls",
    ).rename({"chr": "#chr"}).sort(
        "#chr", "pos", "ref", "alt",
    )


def main():
    args = parse_args()

    input_stem = Path(args.input).name.replace(".gz", "")
    output_dir = Path(args.input).parent
    output_path = args.output or str(output_dir / f"{input_stem}.munged.tsv.gz")
    plot_path = args.plot or str(output_dir / f"{input_stem}.af_af.png")
    gnomad_save_path = str(output_dir / f"{input_stem}.gnomad_filtered.tsv")

    print(f"Reading {args.input}...")
    pgc = read_pgc_sumstats(args.input)
    print(f"  {pgc.height} variants loaded")

    if args.gnomad_af_plot:
        if args.gnomad_filtered:
            print(f"Reading filtered gnomAD from {args.gnomad_filtered}...")
            gnomad = read_gnomad_filtered(args.gnomad_filtered)
        else:
            print("Building rsid set...")
            rsids = build_rsid_set(pgc)
            print(f"  {len(rsids)} unique rsids")

            print(f"Streaming gnomAD from {args.gnomad}...")
            gnomad = stream_gnomad_by_rsid(args.gnomad, rsids, gnomad_save_path)

        print(f"  {gnomad.height} gnomAD rows loaded")

        print("Joining with gnomAD...")
        pgc = join_with_gnomad(pgc, gnomad, args.af_col, args.gnomad_source)

        print("Creating AF-AF plot...")
        create_af_af_plot(pgc, args.af_col, plot_path, input_stem)

    print("Writing output...")
    write_sumstat_output(_prepare_output(pgc), output_path)

    print("Done.")


if __name__ == "__main__":
    main()
