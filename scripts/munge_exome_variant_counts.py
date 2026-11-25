#!/usr/bin/env python3
"""Munge a case-control exome variant file into the per-variant allele COUNT layout ASC2 uses.

Sources, both published under
https://storage.googleapis.com/exome-results-browsers-public/downloads/<date>/:
  SCHEMA2  SCHEMA/SCHEMA_v2_variant_results.tsv.bgz (earlier releases: SCHEMA2/SCHEMA2_...)
  BipEx2   BipEx/BipEx_v2_variant_results.tsv.bgz   (earlier releases: BipEx2/BipEx2_...)
One row per variant, unique on locus+alleles: Hail locus "chrN:POS" and alleles '["REF","ALT"]',
gene_id, VEP consequence, hgvsc, hgvsp, missense scores, `group`, ac_case/an_case,
ac_ctrl/an_ctrl and in_analysis. SCHEMA2 also has transcript_id and n_de_novo; BipEx2 has
neither, and has mac, worst_csq_for_variant_canonical and missense_passing instead.

The releases carry NO per-variant p-value or effect size and none is derived here: the counts
are passed through as the source spells them.

What is kept, and why:
  - every variant, whatever its in_analysis value. The flag is written as the source has it:
    true for the variants the source's gene-level analysis used, false or NA for the rest
  - rows whose gene_id is NA or absent from the gencode table are dropped with a count: the
    target table clusters on the gene symbol. Every row without a consequence is one of them,
    and the script stops if that ever fails to hold
  - `group` is "meta" on every row and is dropped
  - n_de_novo is written as de_novo_ac_proband, NA where the source has no such column. No
    source here has sibling, transmitted or untransmitted counts, a variant class or a LOFTEE
    other-splice flag, so those columns are NA
  - hgvsc, mpc and alpha_missense are kept; the other missense scores, mac,
    worst_csq_for_variant_canonical and missense_passing have no column in the target table
    and are dropped
  - consequence keeps the source's VEP spelling, including its `_LC` suffix

The file is read and written one chromosome at a time: all of it in one DataFrame does not fit
the memory of the machine this runs on.

Output, bgzipped + tabixed on the variant (-s2 -b3 -e3), columns in the order of
genetics-results-db/schemas/exome_variant_counts.sql:
  <dataset>_variant_counts.munged.tsv.gz

Usage:
  python3 munge_exome_variant_counts.py --dataset SCHEMA2 \
      --input SCHEMA_v2_variant_results.tsv.bgz --output-dir <dir or gs://...>
"""

import argparse
import io
import resource
import subprocess
import tempfile
from pathlib import Path

import polars as pl

from exome_utils import locus_contig_chunks, read_gencode
from sumstat_utils import chrom_code_expr, upload_to_gcs

# columns a source may lack: the script refuses a file that lacks one its dataset should have
OPTIONAL_COLS = {"transcript_id", "n_de_novo"}
DATASETS = {
    "SCHEMA2": {"trait": "Schizophrenia", "trait_original": "SCZ", "lacks": set()},
    "BipEx2": {"trait": "Bipolar disorder", "trait_original": "BIP", "lacks": OPTIONAL_COLS},
}
COUNT_COLS = ["ac_case", "an_case", "ac_ctrl", "an_ctrl", "n_de_novo"]
TABIX_ARGS = ["-s2", "-b3", "-e3"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", required=True, choices=sorted(DATASETS))
    parser.add_argument("--input", required=True, help="Path to the source variant results .tsv.bgz")
    parser.add_argument("--gencode", default="/mnt/disks/data/gencode.v39.annotation.genes.tsv",
                        help="Gencode annotation genes TSV (default: v39, as the gene burden munges of both datasets)")
    parser.add_argument("--output-dir", help="Output directory, local or gs:// (default: next to the input)")
    parser.add_argument("--max-memory-gb", type=int, default=8, help="Max virtual memory in GB (default: 8)")
    return parser.parse_args()


def munge_chunk(header: bytes, buf: io.BytesIO, dataset: str, gencode: pl.DataFrame, counts: dict) -> pl.DataFrame:
    spec = DATASETS[dataset]
    columns = header.decode().rstrip("\n").split("\t")
    absent = OPTIONAL_COLS - set(columns)
    if absent != spec["lacks"]:
        raise SystemExit(f"{dataset}: expected the file to lack {sorted(spec['lacks'])}, it lacks {sorted(absent)}")
    df = pl.read_csv(
        io.BytesIO(header + buf.getvalue()), separator="\t", null_values=["NA"], infer_schema_length=0,
        schema_overrides={c: pl.Int64 for c in COUNT_COLS if c in columns},
    )
    counts["read"] += df.height
    for value, n in df.group_by("in_analysis").len().iter_rows():
        counts[f"in_analysis={value}"] = counts.get(f"in_analysis={value}", 0) + n
    groups = df["group"].unique().to_list()
    if groups != ["meta"]:
        raise SystemExit(f"group values {groups}; only 'meta' is expected and the column is dropped")

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
    counts["non_canonical"] += n_before - df.height
    if df.select(pl.struct("chr", "pos", "ref", "alt").is_duplicated().any()).item():
        raise SystemExit("variant rows are not unique on chr:pos:ref:alt")

    counts["no_gene_id"] += df["gene_id"].null_count()
    df = df.filter(pl.col("gene_id").is_not_null())
    n_before = df.height
    df = df.join(gencode, left_on="gene_id", right_on="gene_id_base", how="inner")
    counts["not_in_gencode"] += n_before - df.height
    if df["consequence"].null_count():
        raise SystemExit("a row with a gene has no consequence; the target column is NOT NULL")

    na_str = pl.lit(None, dtype=pl.Utf8)
    na_int = pl.lit(None, dtype=pl.Int64)
    return df.select(
        pl.lit(dataset).alias("#dataset"),
        "chr", "pos", "ref", "alt",
        "gene", "gene_id",
        pl.col("transcript_id") if "transcript_id" in columns else na_str.alias("transcript_id"),
        "consequence",
        na_str.alias("variant_class"),
        "hgvsc", "hgvsp", "mpc", "alpha_missense",
        na_str.alias("is_other_splice"),
        "in_analysis",
        (pl.col("n_de_novo") if "n_de_novo" in columns else na_int).alias("de_novo_ac_proband"),
        na_int.alias("de_novo_ac_sibling"),
        na_int.alias("transmitted_ac_proband"),
        na_int.alias("untransmitted_ac_proband"),
        "ac_case", "an_case", "ac_ctrl", "an_ctrl",
        pl.lit(spec["trait"]).alias("trait"),
        pl.lit(spec["trait_original"]).alias("trait_original"),
    ).sort("chr", "pos", "ref", "alt")


def main() -> None:
    args = parse_args()
    max_bytes = args.max_memory_gb * 1024 ** 3
    resource.setrlimit(resource.RLIMIT_AS, (max_bytes, max_bytes))
    output_dir = args.output_dir or str(Path(args.input).parent)
    output_path = f"{output_dir}/{args.dataset}_variant_counts.munged.tsv.gz"
    local_path = f"{tempfile.mkdtemp()}/out.tsv.gz" if output_path.startswith("gs://") else output_path

    print(f"Reading gencode from {args.gencode}...")
    gencode = read_gencode(args.gencode, coordinates=False)
    print(f"  {gencode.height} genes")

    print(f"Reading variant results from {args.input}...")
    counts = {"read": 0, "non_canonical": 0, "no_gene_id": 0, "not_in_gencode": 0}
    n_written = 0
    with open(local_path, "wb") as out, subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=out) as proc:
        for header, buf in locus_contig_chunks(args.input):
            df = munge_chunk(header, buf, args.dataset, gencode, counts)
            df.write_csv(proc.stdin, separator="\t", null_value="NA", include_header=n_written == 0)
            n_written += df.height
    if proc.returncode:
        raise SystemExit(f"bgzip exited {proc.returncode}")
    subprocess.run(["tabix", "-f", *TABIX_ARGS, local_path], check=True)

    for key, n in counts.items():
        print(f"  {key}: {n}")
    print(f"  wrote {n_written} rows to {output_path}")
    if local_path != output_path:
        upload_to_gcs(local_path, output_path)
    print("\nDone.")


if __name__ == "__main__":
    main()
