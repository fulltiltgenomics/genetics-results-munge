#!/usr/bin/env python3
"""Munge the SCHEMA2 variant file into the per-variant allele COUNT layout ASC2 uses.

Source: the SCHEMA v2 variant file published by the release's browser under
https://storage.googleapis.com/exome-results-browsers-public/downloads/<date>/SCHEMA/
(SCHEMA_v2_variant_results.tsv.bgz; earlier releases name it SCHEMA2_variant_results.tsv.bgz),
one row per variant (32,481,236, unique on locus+alleles): Hail locus "chrN:POS" and alleles '["REF","ALT"]', gene_id, VEP consequence,
transcript_id, hgvsc, hgvsp, four missense scores, `group`, ac_case/an_case, ac_ctrl/an_ctrl,
n_de_novo and in_analysis.

The release carries NO per-variant p-value or effect size and none is derived here: the
counts are passed through as the source spells them.

What is kept, and why:
  - only rows with in_analysis = true, the variants the source's gene-level analysis used.
    The other rows (in_analysis false or NA, most of them intronic or UTR) are not loaded
  - rows whose gene_id is NA or absent from the gencode table are dropped with a count:
    the target table clusters on the gene symbol
  - `group` is "meta" on every row and is dropped
  - n_de_novo is written as de_novo_ac_proband. The source has no sibling, transmitted or
    untransmitted counts, no variant class and no LOFTEE other-splice flag, so those
    columns are NA
  - hgvsc, mpc and alpha_missense are kept; misrank_percentile, misfit_s and pop_eve have
    no column in the target table and are dropped
  - consequence keeps the source's VEP spelling, including its `_LC` suffix

Output (dataset SCHEMA2, trait Schizophrenia / SCZ as in the gene burden munge), bgzipped +
tabixed on the variant (-s2 -b3 -e3), columns in the order of
genetics-results-db/schemas/exome_variant_counts.sql:
  SCHEMA2_variant_counts.munged.tsv.gz

Usage:
  python3 munge_schema2_variant_counts.py --input SCHEMA_v2_variant_results.tsv.bgz \
      --output-dir <dir or gs://...>
"""

import argparse
import gzip
import io
import resource
from pathlib import Path

import polars as pl

from exome_utils import read_gencode
from sumstat_utils import chrom_code_expr, write_exome_output

DATASET = "SCHEMA2"
TRAIT = "Schizophrenia"
TRAIT_ORIGINAL = "SCZ"

COUNT_COLS = ["ac_case", "an_case", "ac_ctrl", "an_ctrl", "n_de_novo"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Path to the SCHEMA v2 variant results .tsv.bgz")
    parser.add_argument("--gencode", default="/mnt/disks/data/gencode.v39.annotation.genes.tsv",
                        help="Gencode annotation genes TSV (default: v39, as the SCHEMA2 gene burden munge)")
    parser.add_argument("--output-dir", help="Output directory, local or gs:// (default: next to the input)")
    parser.add_argument("--max-memory-gb", type=int, default=8, help="Max virtual memory in GB (default: 8)")
    return parser.parse_args()


def read_in_analysis_rows(path: str) -> pl.DataFrame:
    """The in_analysis = true rows. The file is filtered as a stream because decompressing
    all of it into a DataFrame needs more memory than the rows kept justify."""
    kept = io.BytesIO()
    n_total = 0
    with gzip.open(path, "rb") as f:
        header = f.readline()
        columns = header.rstrip(b"\n").split(b"\t")
        in_analysis = columns.index(b"in_analysis")
        group = columns.index(b"group")
        kept.write(header)
        for line in f:
            n_total += 1
            fields = line.rstrip(b"\n").split(b"\t")
            if fields[group] != b"meta":
                raise SystemExit(f"group value {fields[group]!r}; only 'meta' is expected and the column is dropped")
            if fields[in_analysis] == b"true":
                kept.write(line)
    kept.seek(0)
    df = pl.read_csv(
        kept, separator="\t", null_values=["NA"], infer_schema_length=0,
        schema_overrides={c: pl.Int64 for c in COUNT_COLS},
    )
    print(f"  {n_total} variants in the file, {df.height} with in_analysis = true")
    return df


def main() -> None:
    args = parse_args()
    max_bytes = args.max_memory_gb * 1024 ** 3
    resource.setrlimit(resource.RLIMIT_AS, (max_bytes, max_bytes))
    output_dir = args.output_dir or str(Path(args.input).parent)

    print(f"Reading gencode from {args.gencode}...")
    gencode = read_gencode(args.gencode, coordinates=False)
    print(f"  {gencode.height} genes")

    print(f"Reading variant results from {args.input}...")
    df = read_in_analysis_rows(args.input)

    df = df.with_columns(
        chrom_code_expr(pl.col("locus").str.split(":").list.first()).cast(pl.Int32, strict=False).alias("chr"),
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

    n_no_gene = df["gene_id"].null_count()
    df = df.filter(pl.col("gene_id").is_not_null())
    n_before = df.height
    df = df.join(gencode, left_on="gene_id", right_on="gene_id_base", how="inner")
    print(f"  dropped {n_no_gene} rows with no gene_id and {n_before - df.height} rows "
          f"whose gene_id is not in gencode")

    na_str = pl.lit(None, dtype=pl.Utf8)
    na_int = pl.lit(None, dtype=pl.Int64)
    out = df.select(
        pl.lit(DATASET).alias("#dataset"),
        "chr", "pos", "ref", "alt",
        "gene", "gene_id", "transcript_id",
        "consequence",
        na_str.alias("variant_class"),
        "hgvsc", "hgvsp", "mpc", "alpha_missense",
        na_str.alias("is_other_splice"),
        pl.col("n_de_novo").alias("de_novo_ac_proband"),
        na_int.alias("de_novo_ac_sibling"),
        na_int.alias("transmitted_ac_proband"),
        na_int.alias("untransmitted_ac_proband"),
        "ac_case", "an_case", "ac_ctrl", "an_ctrl",
        pl.lit(TRAIT).alias("trait"),
        pl.lit(TRAIT_ORIGINAL).alias("trait_original"),
    ).sort("chr", "pos", "ref", "alt")
    write_exome_output(out, f"{output_dir}/{DATASET}_variant_counts.munged.tsv.gz",
                       tabix_args=["-s2", "-b3", "-e3"], mlog10p_col=None)
    print("\nDone.")


if __name__ == "__main__":
    main()
