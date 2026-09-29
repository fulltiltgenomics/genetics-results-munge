#!/usr/bin/env python3
"""Munge the Autism Sequencing Consortium (ASC) exome release into BigQuery-loadable TSVs.

Source: Satterstrom, Auwerx, Fu et al. 2026, "Rare variation illuminates the distinct and
pleiotropic genetic architecture of autism across neuropsychiatric traits" (medRxiv
10.64898/2026.08.24.26360398). Two files, obtained pre-publication:

  ASC2_gene_results.tsv.bgz     one row per gene (20,893, unversioned Ensembl ids), WIDE:
                                bayes_factor, false_discovery_rate, qc_flagged and 42 count
                                columns = 7 variant classes (ptv mis2 mis1 mis0 syn del dup)
                                x {de_novo_<cls>_proband/_sibling, transmitted_/untransmitted_
                                <cls>_proband, <cls>_case/_control}
  ASC2_variant_results.tsv.bgz  one row per SNV/indel (7,045,130, unique on locus+alleles):
                                Hail locus "chrN:POS" and alleles '["REF","ALT"]', gene_id,
                                transcript_id, VEP consequence, variant_class (PTV, Mis2,
                                Mis1, Mis0, synonymous - no CNVs), hgvsc (NA on every row),
                                hgvsp, mpc, alpha_missense, is_other_splice, gnomad_af, and
                                six allele counts (de_novo_ac_proband/_sibling,
                                transmitted_/untransmitted_ac_proband, ac_case, ac_ctrl)

What was established from the bytes, not the paper:
  - build GRCh38 (SAMD11 at chr1:930274), chr-prefixed, no MT, no scaffolds. The gene file
    is autosomal only (TADA was run on the autosomes); the variant file also carries
    chrX and chrY rows, which are kept as chr 23 and 24
  - `group` is "meta" on every row of both files and is dropped
  - the release carries NO p-value, effect size or allele frequency; nothing is derived
    here - every statistic and count is passed through as the source spells it. The
    numeric statistics are read and written as strings for that reason, so "1.2272e+73"
    reaches BigQuery as written
  - the FDR is NOT a function of the Bayes factor: the paper's Bayesian FDR is the running
    mean of (1 - PPA) over BF-ranked genes, so the 6,860 genes floored at BF = 1 carry FDRs
    from 0.766 to 0.816. Consumers rank on fdr, never on bayes_factor
  - gene-file counts are "independent" counts (one variant per person per gene, LOFTEE
    filtering on inherited and case-control PTVs) and do NOT equal sums over the variant
    file (de novo proband PTVs: 6,279 vs 6,470). The two outputs are independent products
  - gene ids match GENCODE v29 (VEP 95) completely; v39 loses 57 of the 20,893
  - 4 qc_flagged genes (MIB1, LZTR1, NABP2, AC097634.4) are in the file and excluded from the
    paper's 253 at FDR < 0.001; they are kept, flagged

Outputs (dataset ASC2, one binary trait ASD), each bgzipped + tabixed:
  ASC2_gene_counts.munged.tsv.gz          LONG: one row per gene x variant_class x
                                          inheritance_mode with n_affected / n_unaffected
                                          (de_novo: proband/sibling; inherited:
                                          transmitted/untransmitted; case_control:
                                          case/control); gene point index (-s5 -b6 -e6)
  ASC2_gene_bayes_results.munged.tsv.gz   one row per gene: bayes_factor, fdr, qc_flagged;
                                          same index
  ASC2_variant_counts.munged.tsv.gz       one row per variant, WIDE counts as in the source;
                                          variant index (-s2 -b3 -e3)

No mlog10p-filtered companion is written (there is no mlog10p), which is why the shared
writer is called with mlog10p_col=None.

Usage:
  python3 munge_asc.py --gene-input ASC2_gene_results.tsv.bgz \
      --variant-input ASC2_variant_results.tsv.bgz --output-dir <dir or gs://...>
"""

import argparse
import resource
from pathlib import Path

import polars as pl

from munge_ibd_exome import read_gencode
from sumstat_utils import write_exome_output

DATASET = "ASC2"
TRAIT = "ASD"

# source column stem -> the spelling the variant file uses, so the two outputs join on it;
# the CNV classes follow rcnv's cnv_type spelling
VARIANT_CLASSES = {
    "ptv": "PTV", "mis2": "Mis2", "mis1": "Mis1", "mis0": "Mis0", "syn": "synonymous",
    "del": "DEL", "dup": "DUP",
}

# inheritance_mode -> (affected column template, unaffected column template)
INHERITANCE_MODES = {
    "de_novo": ("de_novo_{cls}_proband", "de_novo_{cls}_sibling"),
    "inherited": ("transmitted_{cls}_proband", "untransmitted_{cls}_proband"),
    "case_control": ("{cls}_case", "{cls}_control"),
}

GENE_TABIX = ["-s5", "-b6", "-e6"]
VARIANT_TABIX = ["-s2", "-b3", "-e3"]

VARIANT_COUNT_COLS = [
    "de_novo_ac_proband", "de_novo_ac_sibling",
    "transmitted_ac_proband", "untransmitted_ac_proband",
    "ac_case", "ac_ctrl",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gene-input", help="Path to ASC2_gene_results.tsv.bgz")
    parser.add_argument("--variant-input", help="Path to ASC2_variant_results.tsv.bgz")
    parser.add_argument("--gencode", default="/mnt/disks/data/gencode.v29.annotation.genes.tsv",
                        help="Gencode annotation genes TSV (default: v29, the release's VEP 95 annotation)")
    parser.add_argument("--output-dir", help="Output directory, local or gs:// (default: next to the input)")
    parser.add_argument("--max-memory-gb", type=int, default=24, help="Max virtual memory in GB (default: 24)")
    return parser.parse_args()


def join_gencode(df: pl.DataFrame, gencode: pl.DataFrame, what: str) -> pl.DataFrame:
    n_before = df.height
    df = df.join(gencode, left_on="gene_id", right_on="gene_id_base", how="inner")
    if df.height != n_before:
        print(f"  warning: {n_before - df.height} {what} not found in gencode, dropped")
    return df


def read_gene_results(path: str) -> pl.DataFrame:
    count_cols = [
        tmpl.format(cls=cls)
        for cls in VARIANT_CLASSES
        for tmpls in INHERITANCE_MODES.values()
        for tmpl in tmpls
    ]
    df = pl.read_csv(
        path, separator="\t",
        schema_overrides={
            "gene_id": pl.Utf8, "group": pl.Utf8, "bayes_factor": pl.Utf8,
            "false_discovery_rate": pl.Utf8, "qc_flagged": pl.Utf8,
            **{c: pl.Int64 for c in count_cols},
        },
    )
    groups = df["group"].unique().to_list()
    if groups != ["meta"]:
        raise SystemExit(f"gene file has group values {groups}; only 'meta' is expected and the column is dropped")
    if df["gene_id"].n_unique() != df.height:
        raise SystemExit("gene file has repeated gene_id rows")
    return df


def melt_gene_counts(df: pl.DataFrame) -> pl.DataFrame:
    frames = []
    for cls, class_label in VARIANT_CLASSES.items():
        for mode, (affected_tmpl, unaffected_tmpl) in INHERITANCE_MODES.items():
            frames.append(df.select(
                "gene_id",
                pl.lit(class_label).alias("variant_class"),
                pl.lit(mode).alias("inheritance_mode"),
                pl.col(affected_tmpl.format(cls=cls)).alias("n_affected"),
                pl.col(unaffected_tmpl.format(cls=cls)).alias("n_unaffected"),
            ))
    return pl.concat(frames)


def process_gene_results(gene_path: str, gencode: pl.DataFrame, output_dir: str) -> None:
    print(f"Reading gene results from {gene_path}...")
    df = read_gene_results(gene_path)
    print(f"  {df.height} genes")
    df = join_gencode(df, gencode, "genes")

    counts = melt_gene_counts(df)
    counts = counts.join(df.select("gene_id", "gene", "gene_chr", "gene_start_pos", "gene_end_pos"), on="gene_id")
    counts_out = counts.select(
        pl.lit(DATASET).alias("#dataset"),
        pl.lit(TRAIT).alias("trait"),
        "gene", "gene_id",
        pl.col("gene_chr").alias("chr"),
        "gene_start_pos", "gene_end_pos",
        "variant_class", "inheritance_mode", "n_affected", "n_unaffected",
        pl.lit(TRAIT).alias("trait_original"),
    ).sort("chr", "gene_start_pos", "gene_end_pos", "gene_id", "variant_class", "inheritance_mode")
    print(f"  {counts_out.height} gene x class x mode rows")
    write_exome_output(counts_out, f"{output_dir}/{DATASET}_gene_counts.munged.tsv.gz",
                       tabix_args=GENE_TABIX, mlog10p_col=None)

    bayes_out = df.select(
        pl.lit(DATASET).alias("#dataset"),
        pl.lit(TRAIT).alias("trait"),
        "gene", "gene_id",
        pl.col("gene_chr").alias("chr"),
        "gene_start_pos", "gene_end_pos",
        "bayes_factor",
        pl.col("false_discovery_rate").alias("fdr"),
        "qc_flagged",
        pl.lit(TRAIT).alias("trait_original"),
    ).sort("chr", "gene_start_pos", "gene_end_pos", "gene_id")
    print(f"  {bayes_out.height} gene Bayes rows")
    write_exome_output(bayes_out, f"{output_dir}/{DATASET}_gene_bayes_results.munged.tsv.gz",
                       tabix_args=GENE_TABIX, mlog10p_col=None)


def process_variant_results(variant_path: str, gencode: pl.DataFrame, output_dir: str) -> None:
    print(f"Reading variant results from {variant_path}...")
    df = pl.read_csv(
        variant_path, separator="\t",
        # every column but the counts stays a string: the annotations and gnomad_af are
        # passed through as spelled, and the sentinel NA is kept as NA by the writer
        schema_overrides={c: pl.Int64 for c in VARIANT_COUNT_COLS},
        infer_schema_length=0,
        null_values=["NA"],
    )
    print(f"  {df.height} variants")
    groups = df["group"].unique().to_list()
    if groups != ["meta"]:
        raise SystemExit(f"variant file has group values {groups}; only 'meta' is expected and the column is dropped")
    if df["hgvsc"].null_count() != df.height:
        raise SystemExit("hgvsc is no longer all-NA; stop dropping it")

    df = df.with_columns(
        pl.col("locus").str.split(":").list.first().str.replace(r"(?i)^chr", "")
            .str.replace(r"^X$", "23").str.replace(r"^Y$", "24").str.replace(r"^MT?$", "26")
            .cast(pl.Int32, strict=False).alias("chr"),
        pl.col("locus").str.split(":").list.last().cast(pl.Int64).alias("pos"),
        pl.col("alleles").str.strip_chars("[]").str.replace_all('"', "").str.split(",").alias("allele_list"),
    ).with_columns(
        pl.col("allele_list").list.first().alias("ref"),
        pl.col("allele_list").list.last().alias("alt"),
    )
    if (df["allele_list"].list.len() != 2).any():
        raise SystemExit("a variant row does not have exactly two alleles")
    n_before = df.height
    df = df.filter(pl.col("chr").is_not_null())
    if df.height != n_before:
        print(f"  warning: {n_before - df.height} rows on non-canonical contigs dropped")
    if df.select(pl.struct("chr", "pos", "ref", "alt").is_duplicated().any()).item():
        raise SystemExit("variant rows are not unique on chr:pos:ref:alt")

    df = join_gencode(df, gencode.select("gene_id_base", "gene"), "variants")

    out = df.select(
        pl.lit(DATASET).alias("#dataset"),
        "chr", "pos", "ref", "alt",
        "gene", "gene_id", "transcript_id",
        "consequence", "variant_class",
        "hgvsp", "mpc", "alpha_missense", "is_other_splice", "gnomad_af",
        *VARIANT_COUNT_COLS,
        pl.lit(TRAIT).alias("trait"),
        pl.lit(TRAIT).alias("trait_original"),
    ).sort("chr", "pos", "ref", "alt")
    write_exome_output(out, f"{output_dir}/{DATASET}_variant_counts.munged.tsv.gz",
                       tabix_args=VARIANT_TABIX, mlog10p_col=None)


def main() -> None:
    args = parse_args()
    max_bytes = args.max_memory_gb * 1024 ** 3
    resource.setrlimit(resource.RLIMIT_AS, (max_bytes, max_bytes))

    if not args.gene_input and not args.variant_input:
        raise SystemExit("at least one of --gene-input or --variant-input is required")
    output_dir = args.output_dir or str(Path(args.gene_input or args.variant_input).parent)

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
