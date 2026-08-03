"""
Convert the published PGC schizophrenia FINEMAP 95 % credible sets (Trubetskoy et al. 2022,
supplementary table ST11a) to the credible set TSV format used by the rest of this repository.

Input
  --input        ST11a_95_perc_Credible_Sets.tsv, one row per credible set variant, GRCh37
  --sumstats     daner_PGC_SCZ_w3_90_0418b.munged.tsv.gz, the munged wave 3 summary statistics
  --annotation   FinnGen annotated variants (tabix indexed), for most_severe / gene_most_severe

Output
  <output_dir>/<dataset>_cs_95.tsv   unsorted, with a header; the shell driver sorts, bgzips
                                     and tabix indexes it

ST11a is on GRCh37 and carries no ref/alt, so GRCh38 coordinates, alleles and the alternative
allele frequency are taken from the munged wave 3 summary statistics by rsid — the same variant
space the PGC pseudo credible sets live in, which keeps the two directly comparable. Variants
whose rsid neither source resolves are dropped and cs_size is recomputed from the members that
survive.

The wave 3 daner file covers the autosomes only, so the rsids it does not carry — in practice the
five X-linked sets — are looked up in the variant annotation instead. Those rows get no aaf: the
Finnish allele frequency is not the case/control weighted quantity the rest of the column holds.

Effect sizes, standard errors and p-values are kept from ST11a rather than from the summary
statistics: most of these sets come from the *extended* GWAS, whose association statistics are
not the wave 3 ones. The odds ratio is flipped to the alternative allele where ST11a's effect
allele is the reference allele. ST11a prints the odds ratio to two decimals, so beta is quantised
in steps of roughly 0.005 and rounds to exactly 0 for the weakest members.

Two properties of the source have no column in the credible set schema and are therefore lost
here: `expected_causals_k` (FINEMAP was run allowing several causal variants per locus, so PIPs
within one set sum to roughly k, not to 1) and `extended_gwas`. cs_min_r2 has no source at all
and is written as NA.
"""

import argparse
import io
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import polars as pl

CHR_MAP = {"X": "23"}

ST11A_COLUMNS = [
    "index_snp",
    "rsid",
    "chromosome",
    "GWAS_effect_allele",
    "other_allele",
    "or",
    "se",
    "pval",
    "finemap_posterior_probability",
]

SUMSTAT_COLUMNS = ["#chr", "pos", "ref", "alt", "rsid", "af"]

OUTPUT_COLUMNS = [
    "dataset",
    "data_type",
    "trait",
    "trait_original",
    "cell_type",
    "chr",
    "pos",
    "ref",
    "alt",
    "mlog10p",
    "beta",
    "se",
    "pip",
    "cs_id",
    "cs_size",
    "cs_min_r2",
    "aaf",
    "most_severe",
    "gene_most_severe",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="ST11a credible set TSV (GRCh37)")
    parser.add_argument("--sumstats", required=True, help="Munged PGC SCZ wave 3 summary statistics")
    parser.add_argument("--annotation", required=True, help="Tabix indexed FinnGen annotated variants")
    parser.add_argument("--output-dir", default=".", help="Directory to write the unsorted TSV into")
    parser.add_argument("--dataset", default="PGC_SCZ_2022", help="Value of the dataset column")
    parser.add_argument("--trait", default="SCZ", help="Value of the trait and trait_original columns")
    return parser.parse_args()


def _sci(col: str) -> pl.Expr:
    """Format a float column as %.3e in one numpy call rather than per row, keeping nulls null."""
    formatted = pl.col(col).fill_null(0.0).map_batches(
        lambda s: pl.Series(np.char.mod("%.3e", s.to_numpy())), return_dtype=pl.String
    )
    return pl.when(pl.col(col).is_null()).then(None).otherwise(formatted).alias(col)


def read_st11a(path: str) -> pl.DataFrame:
    """Read ST11a, keeping the credible set columns and the GRCh37 chromosome for a sanity check."""
    df = pl.read_csv(
        path,
        separator="\t",
        columns=ST11A_COLUMNS,
        null_values=["NA", "-", ""],
        schema_overrides={"chromosome": pl.Utf8},
    )
    return df.with_columns(
        pl.col("chromosome").replace(CHR_MAP).cast(pl.Int32, strict=False).alias("chr_b37"),
        pl.col("GWAS_effect_allele").str.to_uppercase().alias("effect_allele"),
        pl.col("other_allele").str.to_uppercase().alias("other_allele"),
    ).drop("chromosome")


def read_sumstats(path: str, rsids: pl.DataFrame) -> pl.DataFrame:
    """Read the GRCh38 locus, alleles and alt allele frequency of the rsids we need."""
    return (
        pl.scan_csv(path, separator="\t", null_values=["NA"], schema_overrides={"#chr": pl.Int32})
        .select(SUMSTAT_COLUMNS)
        .rename({"#chr": "chr"})
        .join(rsids.lazy(), on="rsid", how="semi")
        .collect()
    )


def annotation_header(path: str) -> list[str]:
    return subprocess.run(
        ["tabix", "-H", path], capture_output=True, text=True, check=True
    ).stdout.strip().lstrip("#").split("\t")


def read_locus_by_rsid(path: str, rsids: pl.DataFrame) -> pl.DataFrame:
    """Look the GRCh38 locus and alleles of the given rsids up in the variant annotation.

    This is the fallback for rsids the summary statistics do not carry, which in practice means
    chromosome X: the wave 3 daner file is autosomes only, so the five X-linked sets would
    otherwise be lost entirely. The annotation is indexed on position, not rsid, so it has to be
    streamed; awk does the filtering because the file is ~700 MB gzipped and only a few hundred
    rows are wanted. No alternative allele frequency comes out of this path — the rows it produces
    keep a null aaf rather than borrowing the Finnish frequency, which is a different quantity
    from the case/control weighted one every other row carries.
    """
    header = annotation_header(path)
    rsid_field = header.index("rsid") + 1

    with tempfile.NamedTemporaryFile("w", suffix=".rsids") as wanted:
        wanted.write("\n".join(rsids["rsid"].to_list()) + "\n")
        wanted.flush()
        proc = subprocess.run(
            ["bash", "-c",
             f"zcat {shlex.quote(path)} | awk -F'\\t' 'NR==FNR {{ want[$1]; next }} "
             f"${rsid_field} in want' {shlex.quote(wanted.name)} -"],
            capture_output=True,
            text=True,
            check=True,
        )
    if not proc.stdout:
        return pl.DataFrame(schema={"chr": pl.Int32, "pos": pl.Int64, "ref": pl.Utf8,
                                    "alt": pl.Utf8, "rsid": pl.Utf8, "af": pl.Float64})

    return pl.read_csv(
        io.BytesIO(proc.stdout.encode()),
        separator="\t",
        has_header=False,
        new_columns=header,
        null_values=["NA"],
        schema_overrides={"chr": pl.Int32, "pos": pl.Int64},
    ).select("chr", "pos", "ref", "alt", "rsid", pl.lit(None, dtype=pl.Float64).alias("af"))


def read_annotation(path: str, variants: pl.DataFrame) -> pl.DataFrame:
    """Fetch the annotation of the given variants with tabix.

    Unlike the rsid fallback above this can use the index: the regions are the credible set
    positions themselves, so a few thousand single-base lookups replace a full scan.
    """
    regions = "\n".join(
        f"{chrom}\t{pos - 1}\t{pos}"
        for chrom, pos in variants.select("chr", "pos").unique().sort("chr", "pos").iter_rows()
    )
    proc = subprocess.run(
        ["tabix", "-R", "-", path],
        input=regions,
        capture_output=True,
        text=True,
        check=True,
    )
    anno = pl.read_csv(
        io.BytesIO(proc.stdout.encode()),
        separator="\t",
        has_header=False,
        new_columns=annotation_header(path),
        null_values=["NA"],
        schema_overrides={"variant": pl.Utf8},
    ).select(pl.col("variant").alias("variant_id"), "most_severe", "gene_most_severe")

    return anno.join(variants.select("variant_id"), on="variant_id", how="semi").unique(
        subset="variant_id"
    )


def harmonize(st11a: pl.DataFrame, locus: pl.DataFrame) -> pl.DataFrame:
    """Join ST11a to the resolved GRCh38 loci by rsid and orient the effect on the alt allele.

    Multi-allelic sites give an rsid more than one candidate locus, so the allele check is what
    picks the right one and the deduplication happens after it, not before.
    """
    joined = st11a.join(locus, on="rsid", how="inner")

    n_chr_mismatch = joined.filter(pl.col("chr") != pl.col("chr_b37")).height
    if n_chr_mismatch:
        sys.exit(f"{n_chr_mismatch} rsids sit on a different chromosome than ST11a reports")

    joined = joined.with_columns(
        pl.when((pl.col("effect_allele") == pl.col("alt")) & (pl.col("other_allele") == pl.col("ref")))
        .then(pl.lit(False))
        .when((pl.col("effect_allele") == pl.col("ref")) & (pl.col("other_allele") == pl.col("alt")))
        .then(pl.lit(True))
        .otherwise(None)
        .alias("flip")
    )

    n_incompatible = (
        joined.group_by("rsid").agg(pl.col("flip").is_null().all()).filter(pl.col("flip")).height
    )
    if n_incompatible:
        print(f"  {n_incompatible} rsids dropped: no allele-compatible locus")
    joined = joined.filter(pl.col("flip").is_not_null()).unique(subset=["index_snp", "rsid"])

    return joined.with_columns(
        pl.when(pl.col("flip"))
        .then(-pl.col("or").log())
        .otherwise(pl.col("or").log())
        .round(6)
        .alias("beta_num"),
        (-pl.col("pval").log10()).round(4).alias("mlog10p"),
        pl.col("finemap_posterior_probability").round(4).alias("pip"),
        pl.concat_str("chr", "pos", "ref", "alt", separator=":").alias("variant_id"),
    )


def assign_credible_sets(df: pl.DataFrame) -> pl.DataFrame:
    """Give every set a lead-variant cs_id and a size counted over the variants that survived.

    The index SNP of a set is normally one of its members, so its GRCh38 variant id names the set
    the same way the PGC pseudo credible sets are named. Where the index SNP did not survive the
    rsid join — or was never a member, which ST11a allows — the highest-PIP member names it
    instead.
    """
    leads = (
        df.filter(pl.col("rsid") == pl.col("index_snp"))
        .select("index_snp", pl.col("variant_id").alias("lead_variant"))
        .unique(subset="index_snp")
    )
    fallback = (
        df.sort("pip", "mlog10p", descending=True)
        .group_by("index_snp")
        .first()
        .select("index_snp", pl.col("variant_id").alias("fallback_variant"))
    )

    df = df.join(leads, on="index_snp", how="left").join(fallback, on="index_snp", how="left")

    n_fallback = (
        df.filter(pl.col("lead_variant").is_null()).select("index_snp").unique().height
    )
    if n_fallback:
        print(f"  {n_fallback} of {df['index_snp'].n_unique()} sets named after their top-PIP "
              "member because the index SNP is not among the surviving variants")

    return df.with_columns(
        pl.concat_str(
            pl.lit("chr"),
            pl.coalesce("lead_variant", "fallback_variant").str.replace_all(":", "_"),
        ).alias("cs_id"),
        pl.len().over("index_snp").cast(pl.Int32).alias("cs_size"),
    )


def main() -> None:
    args = parse_args()

    print(f"Reading {args.input}...")
    st11a = read_st11a(args.input)
    print(f"  {st11a.height} credible set variants in {st11a['index_snp'].n_unique()} sets")

    rsids = st11a.select("rsid").unique()
    print(f"Reading {args.sumstats}...")
    locus = read_sumstats(args.sumstats, rsids)
    print(f"  {locus['rsid'].n_unique()} of {rsids.height} rsids found")

    missing = rsids.join(locus.select("rsid"), on="rsid", how="anti")
    if missing.height:
        print(f"Looking {missing.height} remaining rsids up in {args.annotation}...")
        recovered = read_locus_by_rsid(args.annotation, missing)
        print(f"  {recovered['rsid'].n_unique()} recovered without an allele frequency")
        locus = pl.concat([locus, recovered])

    print("Harmonizing...")
    cs = harmonize(st11a, locus)
    print(f"  {cs.height} variants kept in {cs['index_snp'].n_unique()} sets")

    cs = assign_credible_sets(cs)

    print(f"Reading annotation from {args.annotation}...")
    anno = read_annotation(args.annotation, cs.select("chr", "pos", "variant_id").unique())
    print(f"  {anno.height} of {cs['variant_id'].n_unique()} variants annotated")

    output_path = Path(args.output_dir) / f"{args.dataset}_cs_95.tsv"
    cs.join(anno, on="variant_id", how="left").with_columns(
        pl.lit(args.dataset).alias("dataset"),
        pl.lit("GWAS").alias("data_type"),
        pl.lit(args.trait).alias("trait"),
        pl.lit(args.trait).alias("trait_original"),
        pl.lit(None, dtype=pl.String).alias("cell_type"),
        _sci("beta_num").alias("beta"),
        _sci("se").alias("se"),
        # wave 3's own alt allele frequency, the same source the PGC pseudo credible sets use
        _sci("af").alias("aaf"),
        pl.lit(None, dtype=pl.Float64).alias("cs_min_r2"),
    ).select(OUTPUT_COLUMNS).write_csv(output_path, separator="\t", null_value="NA")
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
