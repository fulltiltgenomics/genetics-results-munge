#!/usr/bin/env python3
"""Munge ALS exome results and create AF-AF plot against gnomAD.

Input: 41588_2026_2535_MOESM4_ESM.csv (build 38 exome association results)
- Converts OR to beta = log(OR), computes SE from confidence interval
- Computes -log10(p) from beta/se when p-value underflows
- Matches to gnomAD by position to determine ref/alt alleles
- Flips beta and frequencies for SNPs where effect_allele == gnomAD ref
- Does not flip indels (allele representation may differ)
- Outputs only variants with p < 1e-4
"""

import argparse
import subprocess
import tempfile
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from scipy.special import log_ndtr

GNOMAD_PATH = "/mnt/disks/data/gnomad/gnomad.genomes.exomes.v4.0.sites.v2.tsv.bgz"

# p-value floor in the input data
P_FLOOR = 1e-16


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to ALS exome CSV")
    parser.add_argument("--gnomad", default=GNOMAD_PATH, help="Path to gnomAD tabix-indexed file")
    parser.add_argument("--gnomad-source", default="e", choices=["g", "e"],
                        help="Use gnomAD genomes (g) or exomes (e) (default: e)")
    parser.add_argument("--output", help="Output munged TSV path (default: derived from input)")
    parser.add_argument("--plot", help="Output AF-AF plot PNG path (default: derived from input)")
    return parser.parse_args()


def read_als(filepath: str) -> pl.DataFrame:
    """Read ALS exome CSV and compute beta, se, mlog10p."""
    df = pl.read_csv(
        filepath,
        null_values=["NA", ""],
        schema_overrides={
            "POS": pl.Int64,
            "OR": pl.Float64,
            "ORCIlower": pl.Float64,
            "ORCIupper": pl.Float64,
            "P": pl.Float64,
            "caseMAF": pl.Float64,
            "ctrlMAF": pl.Float64,
            "caseN": pl.Int64,
            "ctrlN": pl.Int64,
        },
    )

    # strip "chr" prefix and convert to int
    df = df.with_columns(
        pl.col("CHROM").str.replace("chr", "")
            .str.replace("X", "23").str.replace("Y", "24").str.replace("MT", "26")
            .cast(pl.Int32).alias("chr"),
    )

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        # beta = log(OR), se from 95% CI
        df = df.with_columns(
            np.log(pl.col("OR")).alias("beta"),
            ((np.log(pl.col("ORCIupper")) - np.log(pl.col("ORCIlower"))) / (2 * 1.96)).alias("se"),
        )

        # mlog10p: use beta/se when p-value is at the floor
        df = df.with_columns(
            pl.when(pl.col("P").is_not_null() & (pl.col("P") > P_FLOOR))
                .then((-np.log10(pl.col("P"))).round(4))
            .when(pl.col("P").is_not_null() & pl.col("se").is_not_null() & (pl.col("se") > 0))
                .then(
                    ((-log_ndtr(-(pl.col("beta") / pl.col("se")).abs()) - np.log(2)) / np.log(10)).round(4)
                )
            .otherwise(None)
            .alias("mlog10p"),
        )

    # classify as SNP or indel
    df = df.with_columns(
        ((pl.col("effect_allele").str.len_chars() == 1) & (pl.col("other_allele").str.len_chars() == 1))
            .alias("is_snp"),
    )

    return df


def query_gnomad(df: pl.DataFrame, gnomad_path: str, gnomad_source: str) -> pl.DataFrame:
    """Query gnomAD by position using tabix and return matched rows."""
    # build regions file for tabix -R
    positions = df.select(
        pl.col("chr").cast(pl.Utf8),
        (pl.col("POS") - 1).cast(pl.Utf8).alias("start"),
        pl.col("POS").cast(pl.Utf8).alias("end"),
    ).unique()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
        regions_path = f.name
        for row in positions.iter_rows():
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")

    print(f"  querying gnomAD for {positions.height} unique positions...")
    proc = subprocess.run(
        ["tabix", "-R", regions_path, gnomad_path],
        capture_output=True, text=True,
    )
    Path(regions_path).unlink()

    if not proc.stdout.strip():
        print("  warning: no gnomAD matches found")
        return pl.DataFrame()

    # parse tabix output
    header = ["#chr", "pos", "ref", "alt", "rsids", "filters", "AN", "AF",
              "AF_afr", "AF_amr", "AF_asj", "AF_eas", "AF_fin", "AF_mid",
              "AF_nfe", "AF_remaining", "AF_sas", "most_severe",
              "gene_most_severe", "consequences", "genome_or_exome"]
    lines = proc.stdout.strip().split("\n")
    gnomad = pl.DataFrame(
        [line.split("\t") for line in lines],
        schema={col: pl.Utf8 for col in header},
        orient="row",
    )
    gnomad = gnomad.with_columns(
        pl.col("#chr").cast(pl.Int32).alias("gnomad_chr"),
        pl.col("pos").cast(pl.Int64).alias("gnomad_pos"),
        pl.col("AF").cast(pl.Float64).alias("gnomad_af"),
        pl.col("AF_nfe").cast(pl.Float64).alias("gnomad_af_nfe"),
    ).rename({"ref": "gnomad_ref", "alt": "gnomad_alt",
              "most_severe": "gnomad_most_severe",
              "gene_most_severe": "gnomad_gene_most_severe"})

    # filter to requested source
    gnomad = gnomad.filter(pl.col("genome_or_exome") == gnomad_source)
    print(f"  {gnomad.height} gnomAD rows ({('exomes' if gnomad_source == 'e' else 'genomes')})")

    return gnomad


def match_alleles(df: pl.DataFrame, gnomad: pl.DataFrame) -> pl.DataFrame:
    """Join ALS data with gnomAD by position and determine allele orientation."""
    joined = df.join(
        gnomad.select("gnomad_chr", "gnomad_pos", "gnomad_ref", "gnomad_alt",
                       "gnomad_af", "gnomad_af_nfe",
                       "gnomad_most_severe", "gnomad_gene_most_severe"),
        left_on=["chr", "POS"],
        right_on=["gnomad_chr", "gnomad_pos"],
        how="left",
    )

    # for SNPs: classify allele match
    joined = joined.with_columns(
        pl.when(~pl.col("is_snp"))
            .then(pl.lit("indel"))
        .when(pl.col("gnomad_ref").is_null())
            .then(pl.lit("no_gnomad"))
        .when(
            (pl.col("effect_allele").str.to_uppercase() == pl.col("gnomad_alt").str.to_uppercase())
            & (pl.col("other_allele").str.to_uppercase() == pl.col("gnomad_ref").str.to_uppercase())
        ).then(pl.lit("effect_is_alt"))
        .when(
            (pl.col("effect_allele").str.to_uppercase() == pl.col("gnomad_ref").str.to_uppercase())
            & (pl.col("other_allele").str.to_uppercase() == pl.col("gnomad_alt").str.to_uppercase())
        ).then(pl.lit("effect_is_ref"))
        .otherwise(pl.lit("mismatch"))
        .alias("allele_match"),
    )

    n_alt = joined.filter(pl.col("allele_match") == "effect_is_alt").height
    n_ref = joined.filter(pl.col("allele_match") == "effect_is_ref").height
    n_mismatch = joined.filter(pl.col("allele_match") == "mismatch").height
    n_no_gnomad = joined.filter(pl.col("allele_match") == "no_gnomad").height
    n_indel = joined.filter(pl.col("allele_match") == "indel").height
    print(f"  allele matching: effect=alt {n_alt}, effect=ref (flip) {n_ref}, "
          f"mismatch {n_mismatch}, no gnomAD {n_no_gnomad}, indel {n_indel}")

    # drop mismatches
    result = joined.filter(pl.col("allele_match") != "mismatch")

    # set ref/alt based on gnomAD for SNPs, use other_allele/effect_allele for indels
    # for indels, treat other_allele as ref and effect_allele as alt (convention from the study)
    is_flip = pl.col("allele_match") == "effect_is_ref"
    result = result.with_columns(
        # ref allele
        pl.when(pl.col("is_snp") & pl.col("gnomad_ref").is_not_null())
            .then(pl.col("gnomad_ref"))
            .otherwise(pl.col("other_allele"))
            .alias("ref"),
        # alt allele
        pl.when(pl.col("is_snp") & pl.col("gnomad_alt").is_not_null())
            .then(pl.col("gnomad_alt"))
            .otherwise(pl.col("effect_allele"))
            .alias("alt"),
        # flip beta for SNPs where effect_allele is gnomAD ref
        pl.when(is_flip).then(-pl.col("beta")).otherwise(pl.col("beta")).alias("beta"),
        # flip case/control MAF
        pl.when(is_flip).then(1.0 - pl.col("caseMAF")).otherwise(pl.col("caseMAF")).alias("caseMAF"),
        pl.when(is_flip).then(1.0 - pl.col("ctrlMAF")).otherwise(pl.col("ctrlMAF")).alias("ctrlMAF"),
    )

    return result


def create_af_af_plot(df: pl.DataFrame, output_path: str) -> None:
    """Create AF-AF scatter plot: study alt AF vs gnomAD alt AF, indels and SNPs separate."""
    # only plot rows that have gnomAD AF
    plot_df = df.filter(pl.col("gnomad_af").is_not_null())

    # compute overall AF = weighted average of case and control MAF
    plot_df = plot_df.with_columns(
        ((pl.col("caseMAF") * pl.col("caseN") + pl.col("ctrlMAF") * pl.col("ctrlN"))
         / (pl.col("caseN") + pl.col("ctrlN"))).alias("study_af"),
    )

    snps = plot_df.filter(pl.col("is_snp"))
    indels = plot_df.filter(~pl.col("is_snp"))

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    for ax, subset, label in [(axes[0], snps, "SNPs"), (axes[1], indels, "Indels")]:
        if subset.height == 0:
            ax.set_title(f"{label} (n=0)")
            continue

        ax.scatter(
            subset["study_af"].to_numpy(),
            subset["gnomad_af"].to_numpy(),
            alpha=0.3, s=2, c="blue",
            rasterized=True,
        )
        ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=0.5)
        ax.set_xlabel("ALS study alt AF")
        ax.set_ylabel("gnomAD AF")
        ax.set_title(f"{label} (n={subset.height:,})")
        ax.set_xlim(0, max(0.01, subset["study_af"].max() * 1.1, subset["gnomad_af"].max() * 1.1))
        ax.set_ylim(0, max(0.01, subset["study_af"].max() * 1.1, subset["gnomad_af"].max() * 1.1))
        ax.set_aspect("equal")

    fig.suptitle(f"ALS exome AF vs gnomAD AF (n={plot_df.height:,})")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  plot saved to {output_path}")


def write_output(df: pl.DataFrame, output_path: str) -> None:
    """Write munged TSV through bgzip."""
    # compute overall AF
    df = df.with_columns(
        ((pl.col("caseMAF") * pl.col("caseN") + pl.col("ctrlMAF") * pl.col("ctrlN"))
         / (pl.col("caseN") + pl.col("ctrlN"))).alias("af_overall_raw"),
    )

    # format numeric columns
    for out_col, src_col in [("beta_fmt", "beta"), ("se_fmt", "se"),
                              ("af_overall", "af_overall_raw"),
                              ("af_cases", "caseMAF"), ("af_controls", "ctrlMAF")]:
        df = df.with_columns(
            pl.col(src_col).map_elements(
                lambda x: f"{x:.4e}" if x is not None and np.isfinite(x) else None,
                return_dtype=pl.Utf8,
            ).alias(out_col),
        )

    out = df.select(
        pl.lit("ALS_exome").alias("#dataset"),
        "chr", pl.col("POS").alias("pos"), "ref", "alt",
        pl.col("gene_name").alias("gene"),
        pl.col("HGVSp").alias("annotation"),
        pl.col("gnomad_most_severe").alias("most_severe"),
        pl.col("gnomad_gene_most_severe").alias("gene_most_severe"),
        "mlog10p",
        pl.col("beta_fmt").alias("beta"),
        pl.col("se_fmt").alias("se"),
        "af_overall", "af_cases", "af_controls",
        pl.lit(None, dtype=pl.Utf8).alias("ac"),
        pl.lit(None, dtype=pl.Utf8).alias("an"),
        pl.lit(None, dtype=pl.Utf8).alias("heritability"),
        pl.lit("amyotrophic_lateral_sclerosis").alias("trait"),
    ).sort("chr", "pos", "ref", "alt")

    with subprocess.Popen(
        ["bgzip", "-c"],
        stdin=subprocess.PIPE,
        stdout=open(output_path, "wb"),
    ) as proc:
        out.write_csv(proc.stdin, separator="\t", null_value="NA")

    subprocess.run(["tabix", "-f", "-s2", "-b3", "-e3", output_path], check=True)
    print(f"  wrote {out.height} rows to {output_path}")


def main():
    args = parse_args()

    input_path = Path(args.input)
    output_dir = input_path.parent
    output_path = args.output or str(output_dir / f"{input_path.stem}.munged.tsv.gz")
    plot_path = args.plot or str(output_dir / f"{input_path.stem}.af_af.png")

    print(f"Reading {args.input}...")
    df = read_als(args.input)
    print(f"  {df.height} variants loaded")
    print(f"  SNPs: {df.filter(pl.col('is_snp')).height}, indels: {df.filter(~pl.col('is_snp')).height}")

    p_underflow = df.filter(pl.col("P") <= P_FLOOR)
    print(f"  {p_underflow.height} variants with P <= {P_FLOOR} (underflow)")
    if p_underflow.height > 0:
        print(f"    mlog10p range: {p_underflow['mlog10p'].min():.1f} - {p_underflow['mlog10p'].max():.1f}")

    print("Querying gnomAD...")
    gnomad = query_gnomad(df, args.gnomad, args.gnomad_source)

    if gnomad.height > 0:
        print("Matching alleles...")
        df = match_alleles(df, gnomad)
    else:
        # no gnomAD data, keep original alleles
        df = df.with_columns(
            pl.col("other_allele").alias("ref"),
            pl.col("effect_allele").alias("alt"),
            pl.lit("no_gnomad").alias("allele_match"),
            pl.lit(None, dtype=pl.Utf8).alias("gnomad_af"),
            pl.lit(None, dtype=pl.Utf8).alias("gnomad_most_severe"),
            pl.lit(None, dtype=pl.Utf8).alias("gnomad_gene_most_severe"),
        )

    print("Creating AF-AF plot...")
    create_af_af_plot(df, plot_path)

    # filter to p < 1e-4 (mlog10p > 4)
    filtered = df.filter(pl.col("mlog10p").is_not_null() & (pl.col("mlog10p") > 4))
    print(f"  {filtered.height} variants with p < 1e-4")

    print("Writing output...")
    write_output(filtered, output_path)

    print("Done.")


if __name__ == "__main__":
    main()
