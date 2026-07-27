#!/usr/bin/env python3
"""Munge IBD exome sequencing gene burden and variant results into genebass-compatible format.

Input files contain results for IBD, UC, and CD in a 'group' column.
Outputs separate files per disease and per result type (gene/variant).
Gene annotations file: gencode v45 (GRCh38).
"""

import argparse
import resource
import warnings
from pathlib import Path

import numpy as np
import polars as pl
from scipy.special import log_ndtr

from sumstat_utils import write_exome_output

GROUPS = ["IBD", "CD", "UC"]

TRAIT_MAP = {
    "IBD": "inflammatory_bowel_disease",
    "CD": "crohns_disease",
    "UC": "ulcerative_colitis",
}

# gene burden annotation classes: prefix -> human-readable name
GENE_ANNOTATION_CLASSES = [
    {"annotation": "pLoF",          "p_col": "ptv_0_001_P_meta",   "beta_col": "ptv_0_001_beta_meta",   "het_p_col": "ptv_0_001_het_P_meta"},
    {"annotation": "nonsynonymous", "p_col": "nsyn_0_001_P_meta",  "beta_col": "nsyn_0_001_beta_meta",  "het_p_col": "nsyn_0_001_het_P_meta"},
]

MLOG10P_THRESHOLD = 4

CODING_CONSEQUENCES = {
    "missense_variant", "frameshift_variant", "inframe_deletion", "inframe_insertion",
    "stop_gained", "stop_lost", "start_lost",
    "splice_donor_variant", "splice_acceptor_variant", "splice_region_variant",
    "coding_sequence_variant", "protein_altering_variant",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-input", help="Path to IBD_gene_results.tsv.bgz")
    parser.add_argument("--variant-input", help="Path to IBD_variant_results.tsv.bgz")
    parser.add_argument("--gencode", default="/mnt/disks/data/gencode.v45.annotation.genes.tsv",
                        help="Path to gencode annotation genes TSV (default: v45)")
    parser.add_argument("--output-dir", help="Output directory (default: same as input)")
    parser.add_argument("--max-memory-gb", type=int, default=24, help="Max virtual memory in GB (default: 24)")
    return parser.parse_args()


def read_gencode(path: str) -> pl.DataFrame:
    """Read gencode annotation, strip version from gene_id, convert chr to numeric."""
    df = pl.read_csv(path, separator="\t")
    df = df.with_columns(
        pl.col("gene_id").str.split(".").list.first().alias("gene_id_base"),
        pl.col("chrom").cast(pl.Utf8)
            .str.replace("X", "23").str.replace("Y", "24").str.replace("M", "26")
            .cast(pl.Int32).alias("gene_chr"),
        pl.col("gene_start").cast(pl.Int32).alias("gene_start_pos"),
        pl.col("gene_end").cast(pl.Int32).alias("gene_end_pos"),
        pl.col("gene_name").alias("gene"),
    ).select("gene_id_base", "gene", "gene_chr", "gene_start_pos", "gene_end_pos")
    return df.unique(subset=["gene_id_base"], keep="first")


def read_gene_results(path: str) -> pl.DataFrame:
    """Read IBD gene burden results."""
    float_cols = {col for cls in GENE_ANNOTATION_CLASSES
                  for col in [cls["p_col"], cls["beta_col"], cls["het_p_col"]]}
    float_cols.update({"n_cases", "n_controls", "variant_p_meta", "variant_beta_meta",
                        "variant_af_case", "variant_af_ctrl", "variant_het_p_meta"})
    return pl.read_csv(
        path, separator="\t", null_values=["NA", ""],
        schema_overrides={c: pl.Float64 for c in float_cols},
    )


def melt_gene_to_long(df: pl.DataFrame, group: str) -> pl.DataFrame:
    """Filter to one disease group, melt annotation classes to long format."""
    df = df.filter(pl.col("group") == group)
    frames = []
    for cls in GENE_ANNOTATION_CLASSES:
        sub = df.select(
            "gene_id", "n_cases", "n_controls",
            pl.col(cls["p_col"]).alias("pvalue"),
            pl.col(cls["beta_col"]).alias("beta_raw"),
            pl.col(cls["het_p_col"]).alias("het_p"),
        ).with_columns(pl.lit(cls["annotation"]).alias("annotation"))
        sub = sub.filter(pl.col("pvalue").is_not_null())
        frames.append(sub)
    return pl.concat(frames)


def compute_gene_stats(df: pl.DataFrame) -> pl.DataFrame:
    """Compute se and mlog10p from beta and p-value."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        # se = |beta| / |z| where z = qnorm(p/2)
        pvals = df["pvalue"].to_numpy()
        betas = df["beta_raw"].to_numpy()
        from scipy.stats import norm
        z = np.abs(norm.ppf(pvals / 2))
        z[z == 0] = np.nan
        se = np.where(np.isfinite(betas) & np.isfinite(z), np.abs(betas) / z, np.nan)

        df = df.with_columns(
            pl.Series("se_raw", se, dtype=pl.Float64)
                .map_elements(lambda x: None if (x is None or np.isnan(x)) else x, return_dtype=pl.Float64),
        )

        # mlog10p
        df = df.with_columns(
            pl.when(pl.col("pvalue").is_not_null() & (pl.col("pvalue") > 0))
                .then(pl.max_horizontal((-np.log10(pl.col("pvalue"))).round(4), 0.0))
                .when(pl.col("pvalue").is_not_null() & pl.col("pvalue").eq(0))
                .then(324.0)
                .otherwise(None)
                .alias("mlog10p_burden"),
        )

        # filter out rows where beta is NA
        df = df.filter(pl.col("beta_raw").is_not_null())

        # format beta/se as scientific notation
        df = df.with_columns(
            pl.col("beta_raw").map_elements(lambda x: f"{x:.4e}", return_dtype=pl.Utf8).alias("beta"),
            pl.col("se_raw").map_elements(lambda x: f"{x:.4e}" if x is not None else None, return_dtype=pl.Utf8).alias("se"),
        )

    return df


def build_gene_output(df: pl.DataFrame, gencode: pl.DataFrame, group: str) -> pl.DataFrame:
    """Join with gencode and build final gene output DataFrame."""
    n_before = df.height
    df = df.join(gencode, left_on="gene_id", right_on="gene_id_base", how="inner")
    n_unmatched = n_before - df.height
    if n_unmatched > 0:
        print(f"  warning: {n_unmatched} genes not found in gencode")

    trait = TRAIT_MAP[group]
    out = df.select(
        pl.lit("IBD_exome").alias("#dataset"),
        pl.lit(trait).alias("trait"),
        "gene",
        "gene_id",
        "gene_chr",
        "gene_start_pos",
        "gene_end_pos",
        "annotation",
        "mlog10p_burden",

        "beta",
        "se",
        pl.lit(None, dtype=pl.Utf8).alias("total_variants"),
        pl.lit(None, dtype=pl.Utf8).alias("total_variants_pheno"),
        pl.col("n_cases").cast(pl.Int64),
        pl.col("n_controls").cast(pl.Int64),
        pl.lit(trait).alias("trait_original"),
        pl.lit("NA").alias("flags"),
    ).sort("gene_chr", "gene_start_pos", "gene_end_pos", "annotation")

    return out


def process_gene_results(gene_path: str, gencode: pl.DataFrame, output_dir: str) -> None:
    """Process gene burden results for all disease groups."""
    print(f"Reading gene results from {gene_path}...")
    gene_df = read_gene_results(gene_path)
    print(f"  {gene_df.height} total rows")

    for group in GROUPS:
        print(f"\nProcessing {group} gene results...")
        long = melt_gene_to_long(gene_df, group)
        print(f"  {long.height} rows across {len(GENE_ANNOTATION_CLASSES)} annotation classes")

        long = compute_gene_stats(long)

        out = build_gene_output(long, gencode, group)
        print(f"  {out.height} total output rows")

        write_exome_output(out, f"{output_dir}/IBD_exome_{group}_gene_results.munged.tsv.gz",
                           tabix_args=["-s5", "-b6", "-e6"], mlog10p_col="mlog10p_burden")


def process_variant_results(variant_path: str, gencode: pl.DataFrame, output_dir: str) -> None:
    """Process variant results for all disease groups."""
    # gencode for gene name lookup only
    gencode_names = gencode.select("gene_id_base", "gene")

    for group in GROUPS:
        print(f"\nProcessing {group} variant results...")
        print("  scanning file...")
        df = (
            pl.scan_csv(
                variant_path, separator="\t", null_values=["NA", ""],
                schema_overrides={
                    "P_meta": pl.Float64, "BETA_meta": pl.Float64, "HetP": pl.Float64,
                    "ac_case": pl.Int64, "an_case": pl.Int64,
                    "ac_ctrl": pl.Int64, "an_ctrl": pl.Int64,
                    "cadd": pl.Float64, "revel": pl.Float64,
                    "gnomADv4_1_genome_nfe_frq": pl.Float64,
                    "gnomad_v4_1_genome_nfe_freq": pl.Float64,
                    "splice_ai": pl.Float64,
                },
            )
            .filter(pl.col("group") == group)
            .filter(pl.col("consequence").is_in(CODING_CONSEQUENCES))
            .filter(pl.col("P_meta").is_not_null())
            .collect()
        )
        print(f"  {df.height} coding rows with P_meta")

        if df.height == 0:
            print(f"  no data for {group}, skipping")
            continue

        # parse locus "chrN:POS" -> chr, pos
        df = df.with_columns(
            pl.col("locus").str.replace("chr", "").str.split(":").list.first()
                .str.replace("X", "23").str.replace("Y", "24").str.replace("MT", "26")
                .cast(pl.Int32).alias("chr"),
            pl.col("locus").str.split(":").list.last().cast(pl.Int32).alias("pos"),
        )

        # parse alleles '["REF","ALT"]' -> ref, alt
        df = df.with_columns(
            pl.col("alleles").str.strip_chars("[]").str.replace_all('"', '').str.split(","),
        ).with_columns(
            pl.col("alleles").list.first().alias("ref"),
            pl.col("alleles").list.last().alias("alt"),
        )

        df = compute_variant_stats(df)

        n_before = df.height
        df = df.join(gencode_names, left_on="gene_id", right_on="gene_id_base", how="inner")
        n_unmatched = n_before - df.height
        if n_unmatched > 0:
            print(f"  warning: {n_unmatched} variants not matched to gencode")

        trait = TRAIT_MAP[group]
        out = df.select(
            pl.lit("IBD_exome").alias("#dataset"),
            "chr", "pos", "ref", "alt",
            "gene",
            pl.col("consequence").alias("annotation"),
            "mlog10p",
            pl.col("beta").alias("beta"),
            "se",
            "af_overall", "af_cases", "af_controls",
            "ac", "an",
            pl.lit(None, dtype=pl.Utf8).alias("n_cases"),
            pl.lit(None, dtype=pl.Utf8).alias("n_controls"),
            pl.lit(trait).alias("trait"),
            pl.lit(trait).alias("trait_original"),
        ).sort("chr", "pos", "ref", "alt")

        write_exome_output(out, f"{output_dir}/IBD_exome_{group}_variant_results.munged.tsv.gz",
                           tabix_args=["-s2", "-b3", "-e3"])


def compute_variant_stats(df: pl.DataFrame) -> pl.DataFrame:
    """Compute mlog10p, allele frequencies, and format columns for variant results."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        # allele frequencies
        df = df.with_columns(
            ((pl.col("ac_case") + pl.col("ac_ctrl")) / (pl.col("an_case") + pl.col("an_ctrl"))).alias("af_overall_raw"),
            (pl.col("ac_case") / pl.col("an_case")).alias("af_cases_raw"),
            (pl.col("ac_ctrl") / pl.col("an_ctrl")).alias("af_controls_raw"),
            (pl.col("ac_case") + pl.col("ac_ctrl")).alias("ac"),
            (pl.col("an_case") + pl.col("an_ctrl")).alias("an"),
        )

        # mlog10p from beta/se when available
        df = df.with_columns(
            pl.when(pl.col("BETA_meta").is_not_null() & pl.col("P_meta").is_not_null() & (pl.col("P_meta") > 0))
                .then(pl.max_horizontal((-np.log10(pl.col("P_meta"))).round(4), 0.0))
                .when(pl.col("P_meta").is_not_null() & pl.col("P_meta").eq(0))
                .then(324.0)
                .otherwise(None)
                .alias("mlog10p"),
        )

        # se from beta and p: se = |beta| / |qnorm(p/2)|
        pvals = df["P_meta"].to_numpy()
        betas = df["BETA_meta"].to_numpy()
        from scipy.stats import norm
        z = np.abs(norm.ppf(pvals / 2))
        z[z == 0] = np.nan
        se_arr = np.where(np.isfinite(betas) & np.isfinite(z) & (z > 0), np.abs(betas) / z, np.nan)
        df = df.with_columns(pl.Series("se_raw", se_arr, dtype=pl.Float64))

        # format numeric columns as scientific notation
        for out_col, src_col in [("beta", "BETA_meta"), ("se", "se_raw"),
                                  ("af_overall", "af_overall_raw"),
                                  ("af_cases", "af_cases_raw"),
                                  ("af_controls", "af_controls_raw")]:
            df = df.with_columns(
                pl.col(src_col).map_elements(
                    lambda x: f"{x:.4e}" if x is not None and np.isfinite(x) else None,
                    return_dtype=pl.Utf8,
                ).alias(out_col),
            )

    return df


def main():
    args = parse_args()

    max_bytes = args.max_memory_gb * 1024 ** 3
    resource.setrlimit(resource.RLIMIT_AS, (max_bytes, max_bytes))

    if not args.gene_input and not args.variant_input:
        print("Error: at least one of --gene-input or --variant-input is required")
        return

    input_path = args.gene_input or args.variant_input
    output_dir = args.output_dir or str(Path(input_path).parent)

    print(f"Reading gencode from {args.gencode}...")
    gencode = read_gencode(args.gencode)
    print(f"  {gencode.height} genes")

    if args.gene_input:
        process_gene_results(args.gene_input, gencode, output_dir)

    if args.variant_input:
        process_variant_results(args.variant_input, gencode, output_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
