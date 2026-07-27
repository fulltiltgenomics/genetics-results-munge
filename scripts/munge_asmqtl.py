#!/usr/bin/env python3
"""Munge ASM-QTL (allele-specific methylation QTL) data and create AF-AF plot against gnomAD.

Handles both Data-S1 (CpG methylation QTLs) and Data-S3 (MDS methylation QTLs),
auto-detecting the type from the header. Each row is a variant-target pair.
"""

import argparse
import subprocess
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from scipy.special import log_ndtr

from sumstat_utils import (
    GNOMAD_AF_COLS, GNOMAD_DEFAULT, GNOMAD_KEEP_COLS,
    build_rsid_set, write_exome_output,
)

# gnomAD columns to keep: standard + annotation
GNOMAD_ANNOT_COLS = ["most_severe", "gene_most_severe"]
GNOMAD_COLS = GNOMAD_KEEP_COLS + GNOMAD_ANNOT_COLS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to Data-S1.tab or Data-S3.tab")
    parser.add_argument("--gnomad", default=GNOMAD_DEFAULT, help="GCS path to gnomAD bgz file")
    parser.add_argument("--af-col", default="AF", choices=GNOMAD_AF_COLS, help="gnomAD AF column for plot (default: AF)")
    parser.add_argument("--output", help="Output munged TSV path (default: derived from input)")
    parser.add_argument("--gnomad-af-plot", action="store_true", help="Stream gnomAD and create AF-AF plot")
    parser.add_argument("--gnomad-filtered", help="Path to previously saved filtered gnomAD TSV (skip streaming)")
    parser.add_argument("--gnomad-source", default="g", choices=["g", "e"], help="Use gnomAD genomes (g) or exomes (e) (default: g)")
    parser.add_argument("--plot", help="Output AF-AF plot PNG path (default: derived from input)")
    parser.add_argument("--norm-mapping", help="Normalization mapping TSV from normalize_asmqtl.py")
    return parser.parse_args()


MIN_GNOMAD_LINES = 800_000_000  # gnomAD genomes+exomes v4.0 has ~855M lines


def stream_gnomad_with_annot(
    gnomad_path: str, rsids: set[str], save_path: str,
) -> pl.DataFrame:
    """Stream gnomAD keeping matching rsids, including annotation columns.

    Supports both GCS (gs://) and local paths. Validates that the full file
    was read by checking total line count.
    """
    if gnomad_path.startswith("gs://"):
        cmd = f"gcloud storage cat {gnomad_path} | zcat"
    else:
        cmd = f"zcat {gnomad_path}"

    proc = subprocess.Popen(
        cmd, shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, bufsize=4_000_000,
    )

    header = proc.stdout.readline().strip().split("\t")
    col_idx = {}
    for col in GNOMAD_COLS:
        if col in header:
            col_idx[col] = header.index(col)

    rsids_idx = col_idx["rsids"]

    outf = open(save_path, "w")
    outf.write("\t".join(col_idx.keys()) + "\n")

    n_lines = 0
    n_matches = 0
    for line in proc.stdout:
        n_lines += 1
        if n_lines % 10_000_000 == 0:
            print(f"  gnomAD: {n_lines / 1e6:.0f}M lines scanned, {n_matches} matches", flush=True)
        fields = line.split("\t")
        if len(fields) <= rsids_idx:
            continue
        rsid_field = fields[rsids_idx]
        if rsid_field and rsids.intersection(rsid_field.split(",")):
            out_fields = [fields[idx].rstrip("\n") for idx in col_idx.values()]
            outf.write("\t".join(out_fields) + "\n")
            n_matches += 1

    outf.close()
    proc.wait()
    print(f"  gnomAD: done. {n_lines / 1e6:.1f}M lines scanned, {n_matches} matches", flush=True)

    if n_lines < MIN_GNOMAD_LINES:
        import os
        os.remove(save_path)
        raise RuntimeError(
            f"gnomAD stream truncated: only {n_lines / 1e6:.1f}M lines read "
            f"(expected >{MIN_GNOMAD_LINES / 1e6:.0f}M). "
            f"Re-run or use a local gnomAD file with --gnomad."
        )

    print(f"  saved filtered gnomAD to {save_path}")
    return read_gnomad_with_annot(save_path)


def read_gnomad_with_annot(filepath: str) -> pl.DataFrame:
    """Read filtered gnomAD TSV including annotation columns."""
    schema_overrides = {
        "#chr": pl.Int32, "pos": pl.Int32,
        **{c: pl.Float64 for c in GNOMAD_AF_COLS},
    }
    if filepath.endswith(".gz"):
        proc = subprocess.Popen(["zcat", filepath], stdout=subprocess.PIPE)
        df = pl.read_csv(
            proc.stdout, separator="\t", null_values=["NA"],
            schema_overrides=schema_overrides, ignore_errors=True,
        )
        proc.wait()
    else:
        df = pl.read_csv(
            filepath, separator="\t", null_values=["NA"],
            schema_overrides=schema_overrides, ignore_errors=True,
        )
    return df.rename({"#chr": "chr"})


def read_asmqtl(filepath: str) -> tuple[pl.DataFrame, str]:
    """Read ASM-QTL file, auto-detect CpG vs MDS type, compute beta/se/mlog10p/af."""
    cols = pl.read_csv(filepath, separator="\t", n_rows=0, infer_schema_length=0).columns
    if "CpG_start" in cols:
        target_type = "cpg"
        target_start_col = "CpG_start"
        target_end_col = "CpG_end"
        ref_methyl_col = "CpG_ref_methylrate"
        alt_methyl_col = "CpG_alt_methylrate"
        ci_col = "SeqVariant_5mCpG_95CI"
    elif "MDS_start" in cols:
        target_type = "mds"
        target_start_col = "MDS_start"
        target_end_col = "MDS_end"
        ref_methyl_col = "MDS_ref_methylrate"
        alt_methyl_col = "MDS_alt_methylrate"
        ci_col = "SeqVariant_MDS_95CI"
    else:
        raise ValueError("Cannot detect file type: no CpG_start or MDS_start column found")

    df = pl.read_csv(
        filepath,
        separator="\t",
        null_values=["NA"],
        schema_overrides={
            "Chrom": pl.Utf8,
            "SeqVariant_start": pl.Float64,
            "SeqVariant_end": pl.Float64,
            "SeqVariant_MAF": pl.Float64,
            "SeqVariant_mCpG_effectsize": pl.Float64,
            "SeqVariant_5mCpG_pvalue": pl.Float64,
            "n_haplotypes": pl.Float64,
            "SeqVariant_LD_count": pl.Float64,
            target_start_col: pl.Float64,
            target_end_col: pl.Float64,
            ref_methyl_col: pl.Float64,
            alt_methyl_col: pl.Float64,
        },
    )

    # chr: strip "chr" prefix, X->23, cast to Int32
    # numeric columns: read as Float64 (handles scientific notation) then cast to Int32
    df = df.with_columns(
        pl.col("Chrom").str.replace("^chr", "").str.replace("^X$", "23").cast(pl.Int32, strict=False).alias("#chr"),
        pl.col("SeqVariant_start").cast(pl.Int32).alias("SeqVariant_start"),
        pl.col("SeqVariant_end").cast(pl.Int32).alias("SeqVariant_end"),
        pl.col("n_haplotypes").cast(pl.Int32).alias("n_haplotypes"),
        pl.col("SeqVariant_LD_count").cast(pl.Int32).alias("SeqVariant_LD_count"),
        pl.col(target_start_col).cast(pl.Int32).alias(target_start_col),
        pl.col(target_end_col).cast(pl.Int32).alias(target_end_col),
    ).filter(pl.col("#chr").is_not_null())

    # parse CI string "lower,upper" into se
    df = df.with_columns(
        pl.col(ci_col).str.split(",").list.get(0).cast(pl.Float64).alias("_ci_lower"),
        pl.col(ci_col).str.split(",").list.get(1).cast(pl.Float64).alias("_ci_upper"),
    ).with_columns(
        ((pl.col("_ci_upper") - pl.col("_ci_lower")) / (2 * 1.96)).alias("_se_raw"),
    )

    # convert minor-allele effect to alt-allele effect
    is_alt_minor = pl.col("SeqVariant_minor_allele") == "ALT"
    df = df.with_columns(
        pl.when(is_alt_minor)
        .then(pl.col("SeqVariant_MAF"))
        .otherwise(1.0 - pl.col("SeqVariant_MAF"))
        .round(6)
        .alias("af"),
        pl.when(is_alt_minor)
        .then(pl.col("SeqVariant_mCpG_effectsize"))
        .otherwise(-pl.col("SeqVariant_mCpG_effectsize"))
        .alias("_beta_raw"),
    )

    # mlog10p with underflow handling
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        df = df.with_columns(
            pl.when(pl.col("SeqVariant_5mCpG_pvalue").eq(0))
            .then(
                ((-log_ndtr(-(pl.col("_beta_raw") / pl.col("_se_raw")).abs()) - np.log(2)) / np.log(10)).round(4)
            )
            .otherwise((-np.log10(pl.col("SeqVariant_5mCpG_pvalue"))).round(4))
            .alias("mlog10p"),
        )

    # rename to output column names
    df = df.rename({
        "SeqVariant_start": "pos",
        "SeqVariant_ref": "ref",
        "SeqVariant_alt": "alt",
        "SeqVariant_rsname": "rsid",
        target_start_col: "target_start",
        target_end_col: "target_end",
        ref_methyl_col: "ref_methylrate",
        alt_methyl_col: "alt_methylrate",
        "SeqVariant_rank": "variant_rank",
        "SeqVariant_LD_count": "ld_count",
        "SeqVariant_vartype": "vartype",
    })

    return df, target_type


def join_gnomad(
    deduped: pl.DataFrame, gnomad: pl.DataFrame, af_col: str, gnomad_source: str = "g",
) -> pl.DataFrame:
    """Join deduplicated variants with gnomAD on rsid.

    Returns joined DataFrame with AF for plot and annotation columns.
    """
    source_label = "genomes" if gnomad_source == "g" else "exomes"
    gnomad_src = gnomad.filter(pl.col("genome_or_exome") == gnomad_source)
    print(f"  using gnomAD {source_label}: {gnomad_src.height} rows")

    gnomad_exploded = gnomad_src.with_columns(
        pl.col("rsids").str.split(",").alias("rsid_list"),
    ).explode("rsid_list").rename({"rsid_list": "rsid"})

    joined = deduped.join(gnomad_exploded, on="rsid", how="inner", suffix="_gnomad")

    # allele orientation: compare ref/alt directly (WGS data, same reference)
    joined = joined.with_columns(
        pl.when(
            (pl.col("ref").str.to_uppercase() == pl.col("ref_gnomad").str.to_uppercase())
            & (pl.col("alt").str.to_uppercase() == pl.col("alt_gnomad").str.to_uppercase())
        ).then(pl.lit("match"))
        .when(
            (pl.col("ref").str.to_uppercase() == pl.col("alt_gnomad").str.to_uppercase())
            & (pl.col("alt").str.to_uppercase() == pl.col("ref_gnomad").str.to_uppercase())
        ).then(pl.lit("swapped"))
        .otherwise(pl.lit("mismatch"))
        .alias("allele_match")
    )

    n_match = joined.filter(pl.col("allele_match") == "match").height
    n_swapped = joined.filter(pl.col("allele_match") == "swapped").height
    n_mismatch = joined.filter(pl.col("allele_match") == "mismatch").height
    n_unmatched = deduped.height - joined.select("rsid").unique().height
    print(f"  match: {n_match}, swapped (flip): {n_swapped}, mismatch: {n_mismatch}, no gnomAD match: {n_unmatched}")

    result = joined.filter(pl.col("allele_match") != "mismatch")
    is_swap = pl.col("allele_match") == "swapped"
    result = result.with_columns(
        pl.when(is_swap).then(1.0 - pl.col("af")).otherwise(pl.col("af")).round(6).alias("af"),
    )

    return result


def create_af_af_plot(
    df: pl.DataFrame, af_col: str, output_path: str, input_name: str,
) -> None:
    """Create AF-AF scatter plot: study alt AF vs gnomAD alt AF."""
    sample_n = min(200_000, df.height)
    sample = df.sample(n=sample_n, seed=42)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(
        sample["af"].to_numpy(),
        sample[af_col].to_numpy(),
        alpha=0.3, s=1, c="blue", rasterized=True,
    )
    ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=0.5)
    ax.set_xlabel("ASM-QTL alt AF (after flip)")
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
    # ensure annotation columns exist (NA if no gnomAD match)
    for col in GNOMAD_ANNOT_COLS:
        if col not in df.columns:
            df = df.with_columns(pl.lit(None).cast(pl.Utf8).alias(col))

    return df.with_columns(
        pl.col("_beta_raw").map_elements(lambda x: f"{x:.3e}" if x is not None else None, return_dtype=pl.Utf8).alias("beta"),
        pl.col("_se_raw").map_elements(lambda x: f"{x:.3e}" if x is not None else None, return_dtype=pl.Utf8).alias("se"),
    ).select(
        "#chr", "pos", "ref", "alt", "rsid",
        "beta", "se", "mlog10p", "af",
        "most_severe", "gene_most_severe",
        "target_start", "target_end", "ref_methylrate", "alt_methylrate",
        "n_haplotypes", "variant_rank", "ld_count", "vartype",
    ).sort("#chr", "pos", "ref", "alt", "target_start")


def main():
    args = parse_args()

    input_stem = Path(args.input).name.replace(".tab", "")
    output_dir = Path(args.input).parent
    output_path = args.output or str(output_dir / f"{input_stem}.munged.tsv.gz")
    plot_path = args.plot or str(output_dir / f"{input_stem}.af_af.png")
    gnomad_save_path = str(output_dir / f"{input_stem}.gnomad_filtered.tsv")

    print(f"Reading {args.input}...")
    df, target_type = read_asmqtl(args.input)
    n_unique = df.unique(subset=["#chr", "pos", "ref", "alt"]).height
    print(f"  {df.height} variant-{target_type} pairs loaded ({n_unique} unique variants)")

    # apply normalization mapping if provided
    if args.norm_mapping:
        print(f"Applying normalization from {args.norm_mapping}...")
        norm = pl.read_csv(args.norm_mapping, separator="\t",
                           schema_overrides={"pos": pl.Int32, "norm_pos": pl.Int32})
        # chrom in mapping has "chr" prefix, #chr in df is Int32
        norm = norm.with_columns(
            pl.col("chrom").str.replace("^chr", "").str.replace("^X$", "23").cast(pl.Int32).alias("#chr"),
        )
        print(f"  {norm.height} indels to normalize")
        n_before = df.height
        df = df.join(norm, on=["#chr", "pos", "ref", "alt"], how="left")
        assert df.height == n_before, f"row count changed: {n_before} -> {df.height}"
        n_normed = df.filter(pl.col("norm_pos").is_not_null()).height
        df = df.with_columns(
            pl.when(pl.col("norm_pos").is_not_null()).then(pl.col("norm_pos")).otherwise(pl.col("pos")).alias("pos"),
            pl.when(pl.col("norm_ref").is_not_null()).then(pl.col("norm_ref")).otherwise(pl.col("ref")).alias("ref"),
            pl.when(pl.col("norm_alt").is_not_null()).then(pl.col("norm_alt")).otherwise(pl.col("alt")).alias("alt"),
        ).drop("chrom", "norm_pos", "norm_ref", "norm_alt")
        n_unique_after = df.unique(subset=["#chr", "pos", "ref", "alt"]).height
        print(f"  {n_normed} rows normalized ({n_unique_after} unique variants after)")

    if args.gnomad_af_plot:
        deduped = df.filter(pl.col("rsid").is_not_null()).unique(subset=["#chr", "pos", "ref", "alt"])
        print(f"  {deduped.height} unique variants with rsids for gnomAD matching")

        if args.gnomad_filtered:
            print(f"Reading filtered gnomAD from {args.gnomad_filtered}...")
            gnomad = read_gnomad_with_annot(args.gnomad_filtered)
        else:
            rsids = build_rsid_set(deduped)
            print(f"  {len(rsids)} unique rsids")
            print(f"Streaming gnomAD from {args.gnomad}...")
            gnomad = stream_gnomad_with_annot(args.gnomad, rsids, gnomad_save_path)

        print(f"  {gnomad.height} gnomAD rows loaded")

        print("Joining with gnomAD...")
        joined = join_gnomad(deduped, gnomad, args.af_col, args.gnomad_source)

        print("Creating AF-AF plot...")
        create_af_af_plot(joined, args.af_col, plot_path, input_stem)

        # extract variant-level annotations from gnomAD match
        annot = joined.select(
            "#chr", "pos", "ref", "alt", "most_severe", "gene_most_severe",
        ).unique(subset=["#chr", "pos", "ref", "alt"], keep="first")
        n_before = df.height
        df = df.join(annot, on=["#chr", "pos", "ref", "alt"], how="left")
        assert df.height == n_before, f"row count changed after annotation join: {n_before} -> {df.height}"
        n_annotated = df.filter(pl.col("most_severe").is_not_null()).height
        n_variants_annotated = df.filter(pl.col("most_severe").is_not_null()).unique(subset=["#chr", "pos", "ref", "alt"]).height
        print(f"  {n_variants_annotated} unique variants annotated with gnomAD consequence")

    print("Writing output...")
    write_exome_output(_prepare_output(df), output_path,
                       tabix_args=["-s1", "-b2", "-e2"])

    print("Done.")


if __name__ == "__main__":
    main()
