"""Shared pieces of the exome munges: gene and variant level burden results.

What lives here is what does not depend on the source paper: the gencode gene table the results
are joined to, the statistics derived from an odds ratio and a p-value, and the command line and
output layout of the gene burden and SCHEMA variant munges. Each script keeps its own reader and
its own column meanings. The write itself (`write_exome_output`) is in `sumstat_utils.py`.
"""

import argparse
import gzip
import io
import resource
import warnings
from pathlib import Path

import numpy as np
import polars as pl

from sumstat_utils import chrom_code_expr, write_exome_output


def read_gencode(
    path: str | Path,
    *,
    coordinates: bool = True,
    chrom_is_numeric: bool = False,
    unique_on: str = "gene_id_base",
) -> pl.DataFrame:
    """Gene names keyed by unversioned ENSG, with numeric-chromosome coordinates unless
    `coordinates=False`.

    `chrom_is_numeric` is for a mapping file that already numbers its chromosomes (X 23, Y 24,
    MT 26): a value that will not cast is then a non-canonical contig and the gene is dropped
    rather than carried into a BigQuery INT64 column. Without it the letters are mapped here and
    anything else that will not cast raises.
    """
    exprs = [
        pl.col("gene_id").str.split(".").list.first().alias("gene_id_base"),
        pl.col("gene_name").alias("gene"),
    ]
    columns = ["gene_id_base", "gene"]
    if coordinates:
        chrom = pl.col("chrom").cast(pl.Utf8)
        if chrom_is_numeric:
            gene_chr = chrom.cast(pl.Int32, strict=False)
        else:
            gene_chr = chrom_code_expr(chrom).cast(pl.Int32)
        exprs += [
            gene_chr.alias("gene_chr"),
            pl.col("gene_start").cast(pl.Int32).alias("gene_start_pos"),
            pl.col("gene_end").cast(pl.Int32).alias("gene_end_pos"),
        ]
        columns += ["gene_chr", "gene_start_pos", "gene_end_pos"]
    df = pl.read_csv(path, separator="\t").with_columns(exprs)
    if coordinates and chrom_is_numeric:
        df = df.filter(pl.col("gene_chr").is_not_null())
    return df.select(columns).unique(subset=[unique_on], keep="first")


def burden_stats_from_or_p(df: pl.DataFrame) -> pl.DataFrame:
    """Compute beta, se, and mlog10p_burden from the OR and pvalue columns, dropping rows whose
    OR has no finite log.

    se = |log(OR)| / |z|  where z = qnorm(p/2)
    """
    # scipy.stats is slow to import and nothing else here needs it
    from scipy.stats import norm

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        df = df.with_columns(
            pl.when(pl.col("OR").is_not_null() & pl.col("OR").is_finite() & (pl.col("OR") > 0))
                .then(np.log(pl.col("OR")))
                .otherwise(None)
                .alias("beta"),
        )

        pvals = df["pvalue"].to_numpy()
        betas = df["beta"].to_numpy()
        z = np.abs(norm.ppf(pvals / 2))
        # avoid division by zero
        z[z == 0] = np.nan
        se = np.where(np.isfinite(betas) & np.isfinite(z), np.abs(betas) / z, np.nan)
        df = df.with_columns(
            pl.Series("se", se, dtype=pl.Float64)
                .map_elements(lambda x: None if (x is None or np.isnan(x)) else x, return_dtype=pl.Float64),
        )

        df = df.with_columns(
            pl.when(pl.col("pvalue").is_not_null() & (pl.col("pvalue") > 0))
                .then(pl.max_horizontal((-np.log10(pl.col("pvalue"))).round(4), 0.0))
                .when(pl.col("pvalue").is_not_null() & pl.col("pvalue").eq(0))
                .then(324.0)
                .otherwise(None)
                .alias("mlog10p_burden"),
        )

        df = df.filter(pl.col("beta").is_not_null())

        df = df.with_columns(
            pl.col("beta").map_elements(lambda x: f"{x:.4e}", return_dtype=pl.Utf8).alias("beta"),
            pl.col("se").map_elements(lambda x: f"{x:.4e}", return_dtype=pl.Utf8).alias("se"),
        )

    return df


def build_gene_burden_output(
    df: pl.DataFrame,
    gencode: pl.DataFrame,
    *,
    dataset: str,
    trait: str,
    trait_original: str,
    total_variants: pl.Expr,
    n_cases: pl.Expr,
    n_controls: pl.Expr,
    flags: pl.Expr,
) -> pl.DataFrame:
    """Join long-format burden rows with gencode and select the gene burden output columns.

    `total_variants` fills both `total_variants` and `total_variants_pheno`.
    """
    df = df.join(gencode, left_on="gene_id", right_on="gene_id_base", how="inner")

    n_unmatched = df.filter(pl.col("gene_chr").is_null()).height
    if n_unmatched > 0:
        print(f"  warning: {n_unmatched} genes not found in gencode")

    return df.select(
        pl.lit(dataset).alias("#dataset"),
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
        total_variants.alias("total_variants"),
        total_variants.alias("total_variants_pheno"),
        n_cases.alias("n_cases"),
        n_controls.alias("n_controls"),
        pl.lit(trait_original).alias("trait_original"),
        flags.alias("flags"),
    ).sort("gene_chr", "gene_start_pos", "gene_end_pos", "annotation")


def run_gene_burden_munge(
    description: str | None,
    source: str,
    gencode_default: str,
    read,
    melt,
    n_classes: int,
    compute_stats,
    output: dict,
    unit: str = "rows",
) -> None:
    """The whole command line of a gene burden munge that melts a wide table, one column set per
    annotation class, to long format: `read(path)` gives the wide frame, `melt` the long one,
    `compute_stats` adds beta/se/mlog10p_burden, and `output` holds the keyword arguments of
    `build_gene_burden_output`. `unit` is what a wide row is called in the progress line.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--input", required=True, help=f"Path to {source} gene results .bgz file")
    parser.add_argument("--gencode", default=gencode_default,
                        help="Path to gencode annotation genes TSV")
    parser.add_argument("--output-dir", help="Output directory (default: same as input)")
    parser.add_argument("--per-trait-dir", help="If given, also write one unfiltered <trait>.tsv.gz there (local or gs://)")
    args = parser.parse_args()

    output_dir = args.output_dir or str(Path(args.input).parent)

    print(f"Reading gencode from {args.gencode}...")
    gencode = read_gencode(args.gencode)
    print(f"  {gencode.height} genes")

    print(f"Reading {source} from {args.input}...")
    wide = read(args.input)
    print(f"  {wide.height} {unit}")

    print("Melting to long format...")
    long = melt(wide)
    print(f"  {long.height} rows across {n_classes} annotation classes")

    print("Computing stats...")
    long = compute_stats(long)

    print("Building output...")
    out = build_gene_burden_output(long, gencode, **output)
    print(f"  {out.height} total rows")

    print("Writing results...")
    write_exome_output(out, f"{output_dir}/{source}_gene_results.munged.tsv.gz",
                       tabix_args=["-s5", "-b6", "-e6"], mlog10p_col="mlog10p_burden",
                       per_trait_dir=args.per_trait_dir)

    print("Done.")


def run_schema_variant_munge(
    description: str | None, dataset: str, gencode_default: str, read, compute_stats, beta_col: str,
) -> None:
    """The whole command line of a SCHEMA-format variant munge: `read(path)` gives one row per
    coding variant with chr/pos/ref/alt parsed, `compute_stats` adds mlog10p, the allele
    frequencies and counts, and `beta_col` names the column that holds the effect size by then.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--input", required=True, help=f"Path to {dataset} variant results .bgz file")
    parser.add_argument("--gencode", default=gencode_default,
                        help="Path to gencode annotation genes TSV")
    parser.add_argument("--output", help="Output munged TSV path (default: derived from input)")
    args = parser.parse_args()

    output_path = args.output or str(
        Path(args.input).parent / f"{dataset}_variant_results.munged.tsv.gz"
    )

    print(f"Reading gencode from {args.gencode}...")
    gencode = read_gencode(args.gencode, coordinates=False)
    print(f"  {gencode.height} genes")

    print(f"Reading {dataset} variants from {args.input}...")
    df = read(args.input)

    print("Computing stats...")
    df = compute_stats(df)

    print("Building output...")
    n_before = df.height
    df = df.join(gencode, left_on="gene_id", right_on="gene_id_base", how="inner")
    n_unmatched = n_before - df.height
    if n_unmatched > 0:
        print(f"  warning: {n_unmatched} variants not matched to gencode")

    out = df.select(
        pl.lit(dataset).alias("#dataset"),
        "chr", "pos", "ref", "alt",
        "gene",
        pl.col("consequence").alias("annotation"),
        "mlog10p",
        pl.col(beta_col).alias("beta"),
        "se",
        "af_overall", "af_cases", "af_controls",
        "ac", "an",
        pl.lit(None, dtype=pl.Utf8).alias("n_cases"),
        pl.lit(None, dtype=pl.Utf8).alias("n_controls"),
        pl.lit("Schizophrenia").alias("trait"),
        pl.lit("SCZ").alias("trait_original"),
    ).sort("chr", "pos", "ref", "alt")

    print("Writing results...")
    write_exome_output(out, output_path, tabix_args=["-s2", "-b3", "-e3"])

    print("Done.")


def run_gene_and_variant_munge(args, process_gene_results, process_variant_results) -> None:
    """The entry point of a munge that takes gene level results, variant level results or both.

    `args` carries `max_memory_gb`, `gene_input`, `variant_input`, `output_dir` and `gencode`.
    Each callback gets (input path, gencode table, output dir). Giving neither input is a usage
    error and exits non-zero: a wrapper that forgot both arguments must not look like a run
    that had nothing to do.
    """
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


def locus_contig_chunks(path: str):
    """(header, lines of one contig) per contig, in file order, for a gzipped TSV whose first
    column is a Hail locus ("chrN:POS"). A contig that comes back after another one would be
    written as two blocks and break the tabix index, so it stops the run."""
    with gzip.open(path, "rb") as f:
        header = f.readline()
        seen, current, buf = set(), None, io.BytesIO()
        for line in f:
            contig = line[:line.index(b":")]
            if contig != current:
                if current is not None:
                    yield header, buf
                    buf = io.BytesIO()
                if contig in seen:
                    raise SystemExit(f"contig {contig.decode()} is not contiguous in the file")
                seen.add(contig)
                current = contig
            buf.write(line)
        if current is not None:
            yield header, buf
