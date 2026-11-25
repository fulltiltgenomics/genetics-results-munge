#!/usr/bin/env python3
"""Munge the gnomAD v4.1.1 MPC scores and missense constraint regions for BigQuery.

Source: gs://gcp-public-data--gnomad/papers/2026-rmc/ (GRCh38, GENCODE v39 transcripts)
  gnomad_v4.1.1_mpc.tsv.bgz            locus, alleles, transcript, mpc: one row per possible
                                       missense SNV and transcript
  gnomad_v4.1.1_all_mcrs.tsv           one row per missense constraint region of a transcript
  gnomad_v4.1.1_non_outlier_mcrs.tsv   the same rows for the transcripts the release does not
                                       treat as outliers

Nothing is derived: scores and region statistics are written as the source spells them. Added
are the gene of each transcript, from a GENCODE v39 GTF, and on the regions `is_outlier_transcript`,
true where the transcript is absent from the non-outlier file. The script stops if a transcript
is not in the GTF, or if a non-outlier row is not a verbatim row of the all-transcripts file.

Output, columns in the order of genetics-results-db/schemas:
  --product mpc      gnomad_v4.1.1_mpc.munged.tsv.gz, bgzipped + tabixed (-s1 -b2 -e2), read
                     and written one chromosome at a time (the file is 70M rows)
  --product regions  gnomad_v4.1.1_missense_constraint_regions.munged.tsv.gz, bgzipped, no index

Usage:
  python3 munge_gnomad_missense_constraint.py --product mpc --input gnomad_v4.1.1_mpc.tsv.bgz \
      --gtf gencode.v39.annotation.gtf.gz --output-dir <dir or gs://...>
  python3 munge_gnomad_missense_constraint.py --product regions \
      --input gnomad_v4.1.1_all_mcrs.tsv --non-outlier gnomad_v4.1.1_non_outlier_mcrs.tsv \
      --gtf gencode.v39.annotation.gtf.gz --output-dir <dir or gs://...>
"""

import argparse
import gzip
import io
import re
import subprocess
import tempfile
from pathlib import Path

import polars as pl

from exome_utils import locus_contig_chunks
from sumstat_utils import chrom_code_expr, upload_to_gcs

RELEASE = "gnomad_v4.1.1"
ATTR = {name: re.compile(rf'{name} "([^"]+)"') for name in ("gene_id", "transcript_id", "gene_name")}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--product", required=True, choices=["mpc", "regions"])
    parser.add_argument("--input", required=True, help="The MPC .tsv.bgz, or the all-transcripts regions TSV")
    parser.add_argument("--non-outlier", help="The non-outlier regions TSV (--product regions)")
    parser.add_argument("--gtf", required=True, help="GENCODE v39 annotation GTF (.gtf.gz)")
    parser.add_argument("--output-dir", help="Output directory, local or gs:// (default: next to the input)")
    return parser.parse_args()


def read_transcript_genes(gtf: str) -> pl.DataFrame:
    """transcript_id -> gene_id, gene, both ids unversioned. The _PAR_Y copies repeat a
    transcript id and are skipped: the chrX row is the one the id names."""
    rows = []
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if fields[2] != "transcript":
                continue
            transcript = ATTR["transcript_id"].search(fields[8]).group(1)
            if transcript.endswith("_PAR_Y"):
                continue
            rows.append((
                transcript.split(".")[0],
                ATTR["gene_id"].search(fields[8]).group(1).split(".")[0],
                ATTR["gene_name"].search(fields[8]).group(1),
            ))
    df = pl.DataFrame(rows, schema=["transcript_id", "gene_id", "gene"], orient="row")
    if df["transcript_id"].is_duplicated().any():
        raise SystemExit("a transcript id appears twice in the GTF")
    return df


def with_genes(df: pl.DataFrame, genes: pl.DataFrame) -> pl.DataFrame:
    out = df.join(genes, on="transcript_id", how="left")
    missing = out.filter(pl.col("gene_id").is_null())["transcript_id"].unique()
    if missing.len():
        raise SystemExit(f"{missing.len()} transcripts are not in the GTF, e.g. {missing[0]}")
    return out


def chr_code(seqname: pl.Expr) -> pl.Expr:
    """Strict: a contig that is not a chromosome stops the run instead of becoming a null."""
    return chrom_code_expr(seqname).cast(pl.Int32)


def write_plain(df: pl.DataFrame, proc_stdin, header: bool) -> None:
    df.write_csv(proc_stdin, separator="\t", null_value="NA", include_header=header)


def munge_mpc(path: str, genes: pl.DataFrame, local_path: str) -> int:
    n = 0
    with open(local_path, "wb") as out, subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=out) as proc:
        for header, buf in locus_contig_chunks(path):
            df = pl.read_csv(io.BytesIO(header + buf.getvalue()), separator="\t", infer_schema_length=0)
            df = df.rename({"transcript": "transcript_id"}).with_columns(
                chr_code(pl.col("locus").str.split(":").list.first()).alias("chr"),
                pl.col("locus").str.split(":").list.last().cast(pl.Int64).alias("pos"),
                pl.col("alleles").str.strip_chars("[]").str.replace_all('"', "").str.split(",").alias("allele_list"),
            )
            if (df["allele_list"].list.len() != 2).any():
                raise SystemExit("a row does not have exactly two alleles")
            if df["mpc"].cast(pl.Float64).null_count():
                raise SystemExit("a row has no MPC value")
            df = with_genes(df, genes).select(
                pl.col("chr").alias("#chr"), "pos",
                pl.col("allele_list").list.first().alias("ref"),
                pl.col("allele_list").list.last().alias("alt"),
                "gene", "gene_id", "transcript_id", "mpc",
            ).sort("#chr", "pos", "ref", "alt", "transcript_id")
            if df.select(pl.struct("#chr", "pos", "ref", "alt", "transcript_id").is_duplicated().any()).item():
                raise SystemExit("rows are not unique on variant + transcript")
            write_plain(df, proc.stdin, header=n == 0)
            n += df.height
    if proc.returncode:
        raise SystemExit(f"bgzip exited {proc.returncode}")
    subprocess.run(["tabix", "-f", "-s1", "-b2", "-e2", local_path], check=True)
    return n


def munge_regions(path: str, non_outlier_path: str, genes: pl.DataFrame, local_path: str) -> int:
    all_lines = Path(path).read_text().splitlines()
    non_outlier_lines = Path(non_outlier_path).read_text().splitlines()
    if all_lines[0] != non_outlier_lines[0] or not set(non_outlier_lines) <= set(all_lines):
        raise SystemExit("the non-outlier file is not a verbatim subset of the all-transcripts file")
    non_outlier = {line.split("\t", 1)[0] for line in non_outlier_lines[1:]}

    df = pl.read_csv(path, separator="\t", infer_schema_length=0).rename({"transcript": "transcript_id"})
    counted = df.group_by("transcript_id").agg(pl.len().alias("n"), pl.col("n_total_transcript_regions").first())
    if counted.filter(pl.col("n") != pl.col("n_total_transcript_regions").cast(pl.Int64)).height:
        raise SystemExit("n_total_transcript_regions does not match a transcript's row count")
    df = with_genes(df, genes).select(
        chr_code(pl.col("chr")).alias("#chr"),
        pl.col("region_start").cast(pl.Int64), pl.col("region_end").cast(pl.Int64),
        "gene", "gene_id", "transcript_id",
        pl.col("n_total_transcript_regions").alias("n_transcript_regions"),
        "region_obs", "region_exp", "region_oe", "region_oe_chisq", "region_oe_chisq_p",
        "is_high_coverage",
        pl.col("transcript_id").is_in(list(non_outlier)).not_().alias("is_outlier_transcript"),
    ).sort("#chr", "region_start", "transcript_id")
    with open(local_path, "wb") as out, subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=out) as proc:
        write_plain(df, proc.stdin, header=True)
    if proc.returncode:
        raise SystemExit(f"bgzip exited {proc.returncode}")
    print(f"  {df['transcript_id'].n_unique()} transcripts, {df['is_outlier_transcript'].sum()} rows on outlier transcripts")
    return df.height


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir or str(Path(args.input).parent)
    name = "mpc" if args.product == "mpc" else "missense_constraint_regions"
    output_path = f"{output_dir}/{RELEASE}_{name}.munged.tsv.gz"
    local_path = f"{tempfile.mkdtemp()}/out.tsv.gz" if output_path.startswith("gs://") else output_path

    print(f"Reading transcripts from {args.gtf}...")
    genes = read_transcript_genes(args.gtf)
    print(f"  {genes.height} transcripts")
    print(f"Reading {args.input}...")
    if args.product == "mpc":
        n = munge_mpc(args.input, genes, local_path)
    else:
        if not args.non_outlier:
            raise SystemExit("--product regions needs --non-outlier")
        n = munge_regions(args.input, args.non_outlier, genes, local_path)
    print(f"  wrote {n} rows to {output_path}")
    if local_path != output_path:
        upload_to_gcs(local_path, output_path)
    print("\nDone.")


if __name__ == "__main__":
    main()
