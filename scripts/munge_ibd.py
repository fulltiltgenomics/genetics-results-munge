#!/usr/bin/env python3
"""Munge IBD/CD/UC meta-analysis GWAS summary statistics and create AF-AF plot against gnomAD."""

import argparse
import resource
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import random

import numpy as np
import polars as pl
from scipy.special import log_ndtr

from sumstat_utils import (
    GNOMAD_AF_COLS, GNOMAD_DEFAULT, GNOMAD_KEEP_COLS,
    read_gnomad_filtered, upload_to_gcs,
)


def _log_memory():
    """Print current RSS memory usage."""
    rss_bytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    print(f"  memory: {rss_bytes / 1024**3:.1f} GB RSS")


CHROMOSOMES = [str(c) for c in range(1, 23)] + ["X"]
CHR_MAP = {"X": 23}
FILENAME_SUFFIX_TEMPLATE = "_{phenotype}_meta_eur_eas_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz"

# columns to rename from input to output
RENAME_MAP = {
    "A1": "ref",
    "A2": "alt",
    "avgA2FREQ": "af",
    "avgA2FREQ_CASES": "af_cases",
    "avgA2FREQ_CONTROLS": "af_controls",
    "N_CASES": "n_cases",
    "N_CONTROLS": "n_controls",
}

# columns that are redundant to derived columns and should be dropped
DROP_COLS = {"MarkerName", "Position_b38", "P-value", "BETA", "SE", "avgA2FREQ", "avgA2FREQ_CASES", "avgA2FREQ_CONTROLS"}

# output column order: derived columns first, then kept original columns
DERIVED_COLS = ["#chr", "pos", "ref", "alt", "beta", "se", "mlog10p", "af", "af_cases", "af_controls"]
RENAMED_COLS = ["INFO", "n_cases", "n_controls"]
KEPT_COLS = ["neff", "Direction_ed", "HetISq", "HetChiSq", "HetDf", "HetPVal", "TotalSampleSize", "N", "rate_total_sample", "rate_Neff"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phenotype", default="ibd", help="Phenotype name used in filenames and labels (default: ibd)")
    parser.add_argument("--input-dir", required=True, help="Directory containing per-chromosome files")
    parser.add_argument("--output", help="Output munged TSV path (default: derived from input-dir)")
    parser.add_argument("--skip-munge", action="store_true", help="Skip munging, use existing munged output")
    parser.add_argument("--gnomad", default=GNOMAD_DEFAULT, help="GCS path to gnomAD bgz file")
    parser.add_argument("--af-col", default="AF", choices=GNOMAD_AF_COLS, help="gnomAD AF column for plot (default: AF)")
    parser.add_argument("--gnomad-af-plot", action="store_true", help="Stream gnomAD and create AF-AF plot")
    parser.add_argument("--gnomad-filter-only", action="store_true", help="Stream gnomAD and save filtered TSV, without creating plot")
    parser.add_argument("--gnomad-filtered", help="Path to previously saved filtered gnomAD TSV (skip streaming)")
    parser.add_argument("--gnomad-source", default="g", choices=["g", "e"], help="Use gnomAD genomes (g) or exomes (e) (default: g)")
    parser.add_argument("--plot", help="Output AF-AF plot PNG path (default: derived from output)")
    parser.add_argument("--max-memory-gb", type=int, default=24, help="Max virtual memory in GB (default: 24)")
    return parser.parse_args()


def read_and_transform(filepath: str, chr_int: int) -> pl.DataFrame:
    """Read one chromosome file, transform columns, sort by position."""
    df = pl.read_csv(
        filepath,
        separator="\t",
        null_values=["NA"],
        truncate_ragged_lines=True,
        schema_overrides={
            "Position_b38": pl.Utf8,
            "BETA": pl.Float64,
            "SE": pl.Float64,
            "P-value": pl.Float64,
            "avgA2FREQ": pl.Float64,
            "avgA2FREQ_CASES": pl.Float64,
            "avgA2FREQ_CONTROLS": pl.Float64,
            "INFO": pl.Float64,
            "Neff": pl.Float64,
        },
        ignore_errors=True,
    )

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        df = df.with_columns(
            pl.lit(chr_int).cast(pl.Int32).alias("#chr"),
            pl.col("Position_b38").str.strip_chars().cast(pl.Int32).alias("pos"),
            # mlog10p with underflow handling
            pl.when(pl.col("P-value").eq(0))
            .then(
                ((-log_ndtr(-(pl.col("BETA") / pl.col("SE")).abs()) - np.log(2)) / np.log(10)).round(4)
            )
            .otherwise((-np.log10(pl.col("P-value"))).round(4))
            .alias("mlog10p"),
        )

    # format numeric columns as scientific notation
    for col_name, orig in [("beta", "BETA"), ("se", "SE"), ("af", "avgA2FREQ"),
                           ("af_cases", "avgA2FREQ_CASES"), ("af_controls", "avgA2FREQ_CONTROLS")]:
        df = df.with_columns(
            pl.col(orig).map_elements(lambda x: f"{x:.3e}" if x is not None else None, return_dtype=pl.Utf8).alias(col_name),
        )

    # rename columns
    df = df.rename({"A1": "ref", "A2": "alt", "N_CASES": "n_cases", "N_CONTROLS": "n_controls", "Neff": "neff"})

    # round neff
    df = df.with_columns(pl.col("neff").round(2))

    # select output columns in order, dropping redundant ones
    output_cols = DERIVED_COLS + RENAMED_COLS + KEPT_COLS
    df = df.select(output_cols).sort("pos", "ref", "alt")

    return df


def _write_to_pipe(df: pl.DataFrame, pipe, include_header: bool) -> None:
    """Write DataFrame to a bgzip pipe as TSV."""
    csv_str = df.write_csv(None, separator="\t", null_value="NA")
    if not include_header:
        csv_str = csv_str[csv_str.index("\n") + 1:]
    pipe.write(csv_str.encode())


def munge_chromosomes(input_dir: str, output_path: str, phenotype: str, build_cpra: bool = False) -> set[tuple[int, int]] | None:
    """Process all chromosomes and write full + filtered output.
    If build_cpra is True, accumulates and returns a set of (chr, pos) tuples."""
    is_gcs = output_path.startswith("gs://")
    if is_gcs:
        tmpdir = tempfile.mkdtemp()
        local_full = f"{tmpdir}/full.tsv.gz"
        local_filt = f"{tmpdir}/filtered.tsv.gz"
    else:
        local_full = output_path
        local_filt = output_path.replace(".tsv.gz", ".mlog10p_gt4.tsv.gz")

    full_proc = subprocess.Popen(
        ["bgzip", "-c"],
        stdin=subprocess.PIPE,
        stdout=open(local_full, "wb"),
    )
    filt_proc = subprocess.Popen(
        ["bgzip", "-c"],
        stdin=subprocess.PIPE,
        stdout=open(local_filt, "wb"),
    )

    total_rows = 0
    total_filtered = 0
    cpra_set = set() if build_cpra else None

    for i, chrom in enumerate(CHROMOSOMES):
        chr_int = CHR_MAP[chrom] if chrom in CHR_MAP else int(chrom)
        filename_suffix = FILENAME_SUFFIX_TEMPLATE.format(phenotype=phenotype)
        filepath = f"{input_dir}/{chrom}{filename_suffix}"
        path = Path(filepath)
        if not path.exists():
            print(f"  warning: {filepath} not found, skipping", file=sys.stderr)
            continue

        print(f"  chr{chrom}...", end=" ", flush=True)
        df = read_and_transform(filepath, chr_int)
        include_header = (i == 0)

        _write_to_pipe(df, full_proc.stdin, include_header)
        total_rows += df.height

        if cpra_set is not None:
            cpra_set.update(zip(df["#chr"].to_list(), df["pos"].to_list()))

        filtered = df.filter(pl.col("mlog10p") > 4)
        _write_to_pipe(filtered, filt_proc.stdin, include_header)
        total_filtered += filtered.height

        print(f"{df.height} variants ({filtered.height} with mlog10p > 4)")

    full_proc.stdin.close()
    filt_proc.stdin.close()
    full_proc.wait()
    filt_proc.wait()

    # tabix index
    subprocess.run(["tabix", "-f", "-s1", "-b2", "-e2", local_full], check=True)
    subprocess.run(["tabix", "-f", "-s1", "-b2", "-e2", local_filt], check=True)

    if is_gcs:
        filtered_gcs = output_path.replace(".tsv.gz", ".mlog10p_gt4.tsv.gz")
        upload_to_gcs(local_full, output_path)
        upload_to_gcs(local_filt, filtered_gcs)

    print(f"  wrote {total_rows} variants to {output_path}")
    print(f"  wrote {total_filtered} variants with mlog10p > 4")

    return cpra_set


def read_munged(filepath: str, columns: list[str] | None = None) -> pl.DataFrame:
    """Re-read the munged output. If columns is set, only read those columns.
    Uses zcat pipe for bgzipped files to avoid loading the full decompressed file into memory."""
    if filepath.endswith(".gz"):
        proc = subprocess.Popen(["zcat", filepath], stdout=subprocess.PIPE)
        df = pl.read_csv(
            proc.stdout, separator="\t", null_values=["NA"],
            columns=columns,
            schema_overrides={"#chr": pl.Int32, "pos": pl.Int32},
            ignore_errors=True,
        )
        proc.wait()
        return df
    return pl.read_csv(
        filepath, separator="\t", null_values=["NA"],
        columns=columns,
        schema_overrides={"#chr": pl.Int32, "pos": pl.Int32},
        ignore_errors=True,
    )


def read_munged_sampled(filepath: str, columns: list[str], sample_n: int, seed: int = 42) -> pl.DataFrame:
    """Reservoir-sample rows from a bgzipped munged file without loading it all into memory."""
    rng = random.Random(seed)
    proc = subprocess.Popen(["zcat", filepath], stdout=subprocess.PIPE, text=True, bufsize=4_000_000)

    header = proc.stdout.readline().strip().split("\t")
    col_indices = [header.index(c) for c in columns]

    reservoir: list[list[str]] = []
    n_seen = 0
    for line in proc.stdout:
        n_seen += 1
        fields = line.rstrip("\n").split("\t")
        row = [fields[i] for i in col_indices]
        if n_seen <= sample_n:
            reservoir.append(row)
        else:
            j = rng.randint(0, n_seen - 1)
            if j < sample_n:
                reservoir[j] = row
    proc.wait()

    df = pl.DataFrame({col: [row[i] for row in reservoir] for i, col in enumerate(columns)})
    df = df.with_columns(
        pl.col("#chr").cast(pl.Int32),
        pl.col("pos").cast(pl.Int32),
    )
    print(f"  {n_seen} total variants, sampled {df.height} for plot")
    return df


def build_cpra_set(df: pl.DataFrame) -> set[tuple[int, int]]:
    """Build set of (chr, pos) tuples for fast gnomAD filtering."""
    return set(zip(df["#chr"].to_list(), df["pos"].to_list()))


def stream_gnomad(
    gnomad_path: str,
    cpra_set: set[tuple[int, int]],
    save_path: str,
) -> pl.DataFrame:
    """Stream gnomAD from GCS, keeping only rows with matching (chr, pos)."""
    proc = subprocess.Popen(
        f"gsutil cat {gnomad_path} | zcat",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=4_000_000,
    )

    header_line = proc.stdout.readline().strip()
    header = header_line.split("\t")
    col_idx = {}
    for col in GNOMAD_KEEP_COLS:
        if col in header:
            col_idx[col] = header.index(col)

    chr_idx = col_idx["#chr"]
    pos_idx = col_idx["pos"]
    max_idx = max(col_idx.values())

    out_header = "\t".join(col_idx.keys())
    bgzip_proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=open(save_path, "wb"))
    bgzip_proc.stdin.write((out_header + "\n").encode())

    n_lines = 0
    n_matches = 0
    for line in proc.stdout:
        n_lines += 1
        if n_lines % 10_000_000 == 0:
            print(f"  gnomAD: {n_lines / 1e6:.0f}M lines scanned, {n_matches} matches", flush=True)

        fields = line.split("\t")
        if len(fields) <= max_idx:
            continue

        try:
            chrom = int(fields[chr_idx])
            pos = int(fields[pos_idx])
        except ValueError:
            continue

        if (chrom, pos) in cpra_set:
            out_fields = [fields[idx].rstrip("\n") for idx in col_idx.values()]
            bgzip_proc.stdin.write(("\t".join(out_fields) + "\n").encode())
            n_matches += 1

    bgzip_proc.stdin.close()
    bgzip_proc.wait()
    proc.wait()
    print(f"  gnomAD: done. {n_lines / 1e6:.1f}M lines scanned, {n_matches} matches", flush=True)
    print(f"  saved filtered gnomAD to {save_path}")

    return read_gnomad_filtered(save_path)


def join_with_gnomad(
    ibd: pl.DataFrame, gnomad: pl.DataFrame, af_col: str,
) -> pl.DataFrame:
    """Join study data with gnomAD on chr/pos, classify allele orientation, flip AF if needed.
    Expects gnomAD already filtered to genomes/exomes.
    Returns DataFrame with af_study (Float64) and gnomAD af_col for plotting."""
    ibd = ibd.with_columns(
        pl.col("af").replace("NA", None).cast(pl.Float64).alias("af_study"),
        pl.col("ref").str.to_uppercase().alias("ref_upper"),
        pl.col("alt").str.to_uppercase().alias("alt_upper"),
    )
    gnomad = gnomad.with_columns(
        pl.col("ref").str.to_uppercase().alias("ref_upper"),
        pl.col("alt").str.to_uppercase().alias("alt_upper"),
    )

    # exact allele match: IBD ref/alt == gnomAD ref/alt
    matched = ibd.join(
        gnomad,
        left_on=["#chr", "pos", "ref_upper", "alt_upper"],
        right_on=["chr", "pos", "ref_upper", "alt_upper"],
        how="inner",
        suffix="_gnomad",
    ).with_columns(
        pl.lit(False).alias("flipped"),
    ).select("#chr", "pos", "ref", "alt", "af_study", pl.col(af_col).alias("af_gnomad"), "flipped")

    # swapped alleles: study ref/alt == gnomAD alt/ref — don't flip AF, just flag
    swapped = ibd.join(
        gnomad,
        left_on=["#chr", "pos", "ref_upper", "alt_upper"],
        right_on=["chr", "pos", "alt_upper", "ref_upper"],
        how="inner",
        suffix="_gnomad",
    ).with_columns(
        pl.lit(True).alias("flipped"),
    ).select("#chr", "pos", "ref", "alt", "af_study", pl.col(af_col).alias("af_gnomad"), "flipped")

    result = pl.concat([matched, swapped])

    # classify SNP vs indel
    result = result.with_columns(
        ((pl.col("ref").str.len_chars() == 1) & (pl.col("alt").str.len_chars() == 1)).alias("is_snp"),
    )

    # print counts per category
    n_total = ibd.height
    snp_match = result.filter(pl.col("is_snp") & ~pl.col("flipped")).height
    snp_swap = result.filter(pl.col("is_snp") & pl.col("flipped")).height
    indel_match = result.filter(~pl.col("is_snp") & ~pl.col("flipped")).height
    indel_swap = result.filter(~pl.col("is_snp") & pl.col("flipped")).height
    n_in_gnomad = result.height
    n_no_gnomad = n_total - n_in_gnomad
    print(f"  SNP match: {snp_match}, SNP swap: {snp_swap}")
    print(f"  indel match: {indel_match}, indel swap: {indel_swap}")
    print(f"  not in gnomAD: {n_no_gnomad} (of {n_total} sampled)")

    return result


def create_af_af_plot(
    df: pl.DataFrame,
    plot_path: str,
    tsv_path: str,
    phenotype: str = "ibd",
) -> None:
    """Create AF-AF scatter plot with four categories and write plot data as TSV."""
    sample_n = min(200_000, df.height)
    sample = df.sample(n=sample_n, seed=42)

    categories = [
        (~sample["flipped"] & sample["is_snp"], "SNP match", "blue"),
        (sample["flipped"] & sample["is_snp"], "SNP swap", "orange"),
        (~sample["flipped"] & ~sample["is_snp"], "indel match", "green"),
        (sample["flipped"] & ~sample["is_snp"], "indel swap", "red"),
    ]

    fig, ax = plt.subplots(figsize=(8, 8))

    for mask, label, color in categories:
        subset = sample.filter(mask)
        if subset.height == 0:
            continue
        ax.scatter(
            subset["af_study"].to_numpy(),
            subset["af_gnomad"].to_numpy(),
            alpha=0.3, s=1, c=color, label=f"{label} ({subset.height:,})",
            rasterized=True,
        )

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=0.5)
    pheno_upper = phenotype.upper()
    ax.set_xlabel(f"{pheno_upper} AF (original, not flipped)")
    ax.set_ylabel("gnomAD AF")
    ax.set_title(f"AF-AF: {pheno_upper} meta (n_variants_in_plot={sample.height:,})")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.legend(markerscale=10, loc="upper left")
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  plot saved to {plot_path}")

    df.write_csv(tsv_path, separator="\t", null_value="NA")
    print(f"  wrote {df.height} plot variants to {tsv_path}")


def main():
    args = parse_args()

    max_bytes = args.max_memory_gb * 1024 ** 3
    resource.setrlimit(resource.RLIMIT_AS, (max_bytes, max_bytes))

    input_dir = args.input_dir.rstrip("/")
    phenotype = args.phenotype
    output_path = args.output or f"{input_dir}/{phenotype}_meta.munged.tsv.gz"
    plot_path = args.plot or output_path.replace(".tsv.gz", ".af_af.png")
    plot_tsv_path = plot_path.replace(".png", ".tsv")
    gnomad_save_path = output_path.replace(".tsv.gz", ".gnomad_filtered.tsv.gz")

    cpra_set = None
    if args.skip_munge:
        print(f"Skipping munging, using existing {output_path}")
    else:
        # build cpra_set during munging if we need to stream gnomAD and don't have a pre-filtered file
        need_gnomad_stream = (args.gnomad_af_plot or args.gnomad_filter_only) and not args.gnomad_filtered
        print(f"Munging {phenotype.upper()} chromosomes...")
        cpra_set = munge_chromosomes(input_dir, output_path, phenotype, build_cpra=need_gnomad_stream)

    if args.gnomad_filter_only:
        if args.gnomad_filtered:
            print("--gnomad-filtered already provided, nothing to do for --gnomad-filter-only")
        else:
            if cpra_set is None:
                print(f"Reading munged output to build (chr, pos) set...")
                ibd = read_munged(output_path, columns=["#chr", "pos"])
                cpra_set = build_cpra_set(ibd)
                del ibd
            print(f"  {len(cpra_set)} unique positions in cpra set")
            print(f"Streaming gnomAD from {args.gnomad}...")
            stream_gnomad(args.gnomad, cpra_set, gnomad_save_path)

    elif args.gnomad_af_plot:
        gnomad_cols = ["#chr", "pos", "ref", "alt", "genome_or_exome", args.af_col]
        if args.gnomad_filtered:
            print(f"Reading filtered gnomAD from {args.gnomad_filtered}...")
            gnomad = read_gnomad_filtered(args.gnomad_filtered, columns=gnomad_cols)
        else:
            print(f"  {len(cpra_set)} unique positions in cpra set")
            print(f"Streaming gnomAD from {args.gnomad}...")
            gnomad = stream_gnomad(args.gnomad, cpra_set, gnomad_save_path)
            del cpra_set
            gnomad = read_gnomad_filtered(gnomad_save_path, columns=gnomad_cols)

        # filter to genomes/exomes early and drop the column to free memory
        source_label = "genomes" if args.gnomad_source == "g" else "exomes"
        gnomad = gnomad.filter(pl.col("genome_or_exome") == args.gnomad_source).drop("genome_or_exome")
        print(f"  {gnomad.height} gnomAD {source_label} rows loaded")
        _log_memory()

        ibd_cols = ["#chr", "pos", "ref", "alt", "af"]
        print(f"Re-reading munged output from {output_path} (reservoir sampling)...")
        ibd = read_munged_sampled(output_path, ibd_cols, sample_n=500_000)
        _log_memory()

        print("Joining with gnomAD...")
        joined = join_with_gnomad(ibd, gnomad, args.af_col)
        del ibd, gnomad
        _log_memory()

        print("Creating AF-AF plot...")
        create_af_af_plot(joined, plot_path, plot_tsv_path, phenotype)

    print("Done.")


if __name__ == "__main__":
    main()
