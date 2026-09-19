#!/usr/bin/env python3
"""Munge BRaVa cross-ancestry exome meta-analysis gene burden results into the genebass
gene-burden layout (medRxiv 2026.05.21.26353759).

Source: `gs://daly-genetics-results/raw/brava/gene/<code>_gene_meta_analysis_100_cutoff[.<STRATUM>].tsv.gz`,
one file per (phenotype, meta-analysis stratum), copied whole to disk before reading. Columns:
Region (ENSG, unversioned), Group, max_MAF, Pvalue, Stat, type, df, sum_weights, BETA_Burden,
chisq_het, Pvalue_het, SE_Burden, class.

What the source is, established by reading it rather than from the paper:

- **Only `class == "Burden"` and `type == "Inverse variance weighted"` rows are kept.** They are
  the only rows carrying both an effect size and its standard error, and `beta` is NOT NULL in
  BigQuery's `gene_burden_results`. Dropped: the Stouffer burden rows (a beta, no se, and the same
  gene x mask x MAF cell as the IVW row), and the SKAT and SKAT-O rows (a p-value only). The
  heterogeneity columns `chisq_het`/`Pvalue_het` have no place in the layout and go with them.
- the source carries no gene symbol, no coordinates, no sample sizes and no variant counts, so
  `total_variants`/`total_variants_pheno` are NA and everything else is joined in from elsewhere.
- coordinates and symbols come from GENCODE **v39**, which is what BRaVa annotated its masks
  against (VEP 105 + LOFTEE + SpliceAI over v39); genes whose ENSG does not join are reported
  and dropped. This is a binding constraint, not a version choice: the results-api entry for
  this dataset declares `gencode_version` 39 with a point index on `gene_start_pos`, so
  coordinates from any other version return nothing at query time. The ~34 ENSGs per file
  (~350 rows per trait) that don't join are accessions minted after v39, and are the cost of
  that constraint rather than a sign v39 fits the data best.
- Pvalue underflows to 0 on the strongest cells (APOB x LDLC among them), so the `log_ndtr`
  recovery from beta/se is a live path here rather than a precaution. A row that underflowed
  *and* has no BETA_Burden leaves nothing to recover from and is dropped with the rest of the
  rows that carry no effect size.
- `annotation` is the source's own mask spelling: `Group` with `;` replaced by `|`, then
  `|MAF<` and `max_MAF` **as written in the file** (`1e-04`, `0.001`), e.g.
  `pLoF|damaging_missense_or_protein_altering|MAF<1e-04`.

Trait naming (the epic's "strata in the trait code" decision): `trait_original` is the phenocode
from `brava_pheno.json` -- `<code>` for the cross-ancestry meta and `<code>|<STRATUM>` for a
stratum -- and `trait` is that item's phenostring. `n_cases`/`n_controls` come from the same
item: a binary phenotype's `num_cases`/`num_controls`, and for a quantitative one `num_samples`
as `n_cases` with `n_controls` NA, which is how genebass writes a quantitative trait.

Selectors: `--phenotypes` takes phenocodes, where a bare `CODE` is the cross-ancestry meta,
`CODE|STRATUM` one stratum and `CODE|all` the meta plus every stratum the JSON has for it;
`--strata` adds those strata to every bare code. Both are checked against the JSON.

    python3 scripts/munge_brava.py --phenotypes 'AFib|all' LDLC --output-dir out \\
        --per-trait-dir out/gene_burden_per_trait

Outputs: one combined `BRaVa_gene_results.tsv.gz` over every trait in the run plus the
`mlog10p_burden > 4` companion `write_exome_output` builds beside it, and one unfiltered
`gene_burden_per_trait/<trait_original>.tsv.gz` per trait. `--stage` uploads the filtered
combined file and the per-trait files (the unfiltered combined file is a local by-product of the
shared writer and is not served).
"""

import argparse
import json
import sys
from pathlib import Path

import polars as pl

from brava_phenotypes import source_file_code
from sumstat_utils import fetch_gs, upload_to_gcs, write_exome_output

GENE_PREFIX = "gs://daly-genetics-results/raw/brava/gene/"
PHENO_JSON = "gs://daly-genetics-results/mapping_files/brava_pheno.json"
GENCODE = "gs://daly-genetics-results/mapping_files/gencode.v39.annotation.genes.tsv"
STAGE_PREFIX = "gs://daly-genetics-results/exome_results/brava/"

COMBINED = "BRaVa_gene_results.tsv.gz"
FILTERED = "BRaVa_gene_results.mlog10p_gt4.tsv.gz"
TABIX_ARGS = ["-s5", "-b6", "-e6"]

# max_MAF is read as text because the annotation string carries it verbatim: `1e-04` parsed and
# reformatted would silently become `0.0001` and no longer match the published mask names
SOURCE_COLUMNS = {
    "Region": pl.Utf8,
    "Group": pl.Utf8,
    "max_MAF": pl.Utf8,
    "Pvalue": pl.Float64,
    "type": pl.Utf8,
    "BETA_Burden": pl.Float64,
    "SE_Burden": pl.Float64,
    "class": pl.Utf8,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phenotypes", nargs="+",
                        help="phenocodes to munge: CODE (meta), CODE|STRATUM, or CODE|all "
                             "(default: the cross-ancestry meta of every phenotype in the JSON)")
    parser.add_argument("--strata", nargs="+", default=[],
                        help="strata to add for every bare CODE in --phenotypes")
    parser.add_argument("--gene-prefix", default=GENE_PREFIX, help="prefix holding the gene result files")
    parser.add_argument("--pheno-json", default=PHENO_JSON, help="brava_pheno.json (gs:// or local)")
    parser.add_argument("--gencode", default=GENCODE, help="gencode annotation genes TSV (gs:// or local)")
    parser.add_argument("--output-dir", required=True, help="directory for the combined files")
    parser.add_argument("--per-trait-dir", help="if given, also write one unfiltered <trait_original>.tsv.gz there")
    parser.add_argument("--cache-dir", default=str(Path.home() / "brava_munge" / "cache"),
                        help="where the source files are downloaded to")
    parser.add_argument("--stage", action="store_true", help=f"upload the served outputs to {STAGE_PREFIX}")
    return parser.parse_args()


def read_gencode(path: Path) -> pl.DataFrame:
    """Gene coordinates keyed by unversioned ENSG.

    The mapping file already numbers chromosomes the way the output does (X 23, Y 24, MT 26);
    anything that will not cast is a non-canonical contig and is dropped rather than carried into
    a BigQuery INT64 column.
    """
    df = pl.read_csv(path, separator="\t")
    return (
        df.with_columns(
            pl.col("gene_id").str.split(".").list.first().alias("gene_id_base"),
            pl.col("chrom").cast(pl.Utf8).cast(pl.Int32, strict=False).alias("gene_chr"),
            pl.col("gene_start").cast(pl.Int32).alias("gene_start_pos"),
            pl.col("gene_end").cast(pl.Int32).alias("gene_end_pos"),
            pl.col("gene_name").alias("gene"),
        )
        .filter(pl.col("gene_chr").is_not_null())
        .select("gene_id_base", "gene", "gene_chr", "gene_start_pos", "gene_end_pos")
        .unique(subset=["gene_id_base"], keep="first")
    )


def resolve_traits(requested: list[str], strata: list[str], by_code: dict[str, dict]) -> list[str]:
    """Expand the selectors into phenocodes, in the order asked for, dropping repeats."""
    wanted: list[str] = []
    for entry in requested:
        base, _, stratum = entry.partition("|")
        if stratum == "all":
            wanted += [code for code in by_code if code == base or code.startswith(f"{base}|")]
        elif stratum:
            wanted.append(entry)
        else:
            wanted.append(base)
            wanted += [f"{base}|{one}" for one in strata]
    unknown = [code for code in wanted if code not in by_code]
    if unknown:
        raise SystemExit(f"no such phenocode in the phenotype JSON: {sorted(set(unknown))}")
    return list(dict.fromkeys(wanted))


def gene_file_uri(prefix: str, phenocode: str) -> str:
    base, _, stratum = phenocode.partition("|")
    name = f"{source_file_code(base)}_gene_meta_analysis_100_cutoff"
    return prefix.rstrip("/") + "/" + name + (f".{stratum}" if stratum else "") + ".tsv.gz"


def keep_burden_ivw(df: pl.DataFrame) -> pl.DataFrame:
    """The inverse-variance-weighted burden rows, the only ones with both a beta and an se."""
    return df.filter((pl.col("class") == "Burden") & (pl.col("type") == "Inverse variance weighted"))


def mlog10p_expr(df: pl.DataFrame) -> pl.Expr:
    """-log10(Pvalue), recovered from beta/se for any row where p underflowed to 0."""
    direct = pl.max_horizontal((-pl.col("Pvalue").log10()).round(4), 0.0)
    if df.filter(pl.col("Pvalue") == 0).height == 0:
        return direct
    # only reached by a file that underflows, which is also the only reason this run needs scipy
    from numpy import log
    from scipy.special import log_ndtr

    recovered = (
        (-log_ndtr(-(pl.col("BETA_Burden") / pl.col("SE_Burden")).abs()) - log(2)) / log(10)
    ).round(4)
    return pl.when(pl.col("Pvalue") > 0).then(direct).otherwise(recovered)


def build_output(df: pl.DataFrame, gencode: pl.DataFrame, item: dict) -> tuple[pl.DataFrame, int, list[str]]:
    """Join one file's kept rows to gencode and select the gene-burden columns."""
    n_cases = item.get("num_cases", item.get("num_samples"))
    n_controls = item.get("num_controls")

    df = df.with_columns(
        (pl.col("Group").str.replace_all(";", "|") + pl.lit("|MAF<") + pl.col("max_MAF")).alias("annotation"),
        mlog10p_expr(df).alias("mlog10p_burden"),
    )
    # a degenerate meta cell can carry Pvalue 0 with no BETA_Burden and se 0, leaving nothing to
    # serve and nothing to recover a p from; beta and mlog10p_burden are NOT NULL in BigQuery's
    # gene_burden_results, so such a row is dropped rather than written as NA. SE_Burden == 0
    # with a non-null beta and Pvalue == 0 gives mlog10p_burden = inf, which the finite check
    # catches and folds into the same drop rather than writing it literally
    is_usable = (
        pl.col("BETA_Burden").is_not_null()
        & pl.col("SE_Burden").is_not_null()
        & pl.col("mlog10p_burden").is_not_null()
        & pl.col("mlog10p_burden").is_finite()
    )
    usable = df.filter(is_usable)
    no_stats = df.height - usable.height

    # a p=0 row among the drops may be a real result (e.g. the strongest cell of a flagship
    # gene) losing its only recovery path, not junk -- worth a named line, not just a count
    for row in df.filter(~is_usable & (pl.col("Pvalue") == 0)).iter_rows(named=True):
        print(f"  WARNING: dropped p=0 row with no usable beta/se: trait={item['phenostring']} "
              f"gene_id={row['Region']} annotation={row['annotation']}", file=sys.stderr)

    unmatched = (
        usable.join(gencode, left_on="Region", right_on="gene_id_base", how="anti")
        .select("Region").unique().to_series().sort().to_list()
    )

    out = (
        usable
        .join(gencode, left_on="Region", right_on="gene_id_base", how="inner")
        .select(
            pl.lit("BRaVa").alias("#dataset"),
            pl.lit(item["phenostring"]).alias("trait"),
            "gene",
            pl.col("Region").alias("gene_id"),
            "gene_chr",
            "gene_start_pos",
            "gene_end_pos",
            "annotation",
            "mlog10p_burden",
            pl.col("BETA_Burden").map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8).alias("beta"),
            pl.col("SE_Burden").map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8).alias("se"),
            pl.lit(None, dtype=pl.Utf8).alias("total_variants"),
            pl.lit(None, dtype=pl.Utf8).alias("total_variants_pheno"),
            pl.lit(n_cases, dtype=pl.Int64).alias("n_cases"),
            pl.lit(n_controls, dtype=pl.Int64).alias("n_controls"),
            pl.lit(item["phenocode"]).alias("trait_original"),
            pl.lit("NA").alias("flags"),
        )
    )
    return out, no_stats, unmatched


def munge_trait(uri: str, cache: Path, gencode: pl.DataFrame, item: dict) -> tuple[pl.DataFrame, dict, list[str]]:
    path = fetch_gs(uri, cache)
    source = pl.read_csv(path, separator="\t", null_values=["", "NA"],
                         columns=list(SOURCE_COLUMNS), schema_overrides=SOURCE_COLUMNS)
    kept = keep_burden_ivw(source)
    out, no_stats, unmatched = build_output(kept, gencode, item)
    counts = {
        "rows": source.height,
        "burden_ivw": kept.height,
        "dropped_no_stats": no_stats,
        "dropped_no_gencode": kept.height - no_stats - out.height,
        "mlog10p_gt4": out.filter(pl.col("mlog10p_burden") > 4).height,
    }
    return out, counts, unmatched


def stage(output_dir: str, per_trait_dir: str, traits: list[str]) -> None:
    upload_to_gcs(f"{output_dir}/{FILTERED}", STAGE_PREFIX + FILTERED)
    for trait in traits:
        upload_to_gcs(f"{per_trait_dir}/{trait}.tsv.gz",
                      f"{STAGE_PREFIX}gene_burden_per_trait/{trait}.tsv.gz")


def main() -> None:
    args = parse_args()
    cache = Path(args.cache_dir)

    if args.stage and (not args.per_trait_dir or args.per_trait_dir.startswith("gs://")):
        raise SystemExit("--stage uploads the local outputs, so --per-trait-dir must be a local directory")

    items = json.loads(fetch_gs(args.pheno_json, cache).read_text())
    by_code = {item["phenocode"]: item for item in items}
    traits = resolve_traits(args.phenotypes or [c for c in by_code if "|" not in c], args.strata, by_code)
    print(f"{len(traits)} traits: {', '.join(traits)}")

    print(f"Reading gencode from {args.gencode}...")
    gencode = read_gencode(fetch_gs(args.gencode, cache))
    print(f"  {gencode.height} genes")

    frames = []
    counts = {}
    unmatched: dict[str, list[str]] = {}
    for trait in traits:
        uri = gene_file_uri(args.gene_prefix, trait)
        print(f"\n{trait} <- {uri}")
        out, counts[trait], unmatched[trait] = munge_trait(uri, cache, gencode, by_code[trait])
        print("  " + ", ".join(f"{k}={v}" for k, v in counts[trait].items()))
        if unmatched[trait]:
            print(f"  {len(unmatched[trait])} gene ids not in gencode, dropped: {' '.join(unmatched[trait])}")
        frames.append(out)

    print(f"\n{'trait':<20} {'rows':>10} {'burden IVW':>12} {'no beta':>9} {'no gencode':>12} {'mlog10p>4':>11}")
    for trait in traits:
        row = counts[trait]
        print(f"{trait:<20} {row['rows']:>10,} {row['burden_ivw']:>12,} {row['dropped_no_stats']:>9,} "
              f"{row['dropped_no_gencode']:>12,} {row['mlog10p_gt4']:>11,}")

    combined = pl.concat(frames).sort("gene_chr", "gene_start_pos", "gene_end_pos",
                                      "trait_original", "annotation")
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    write_exome_output(combined, f"{args.output_dir}/{COMBINED}", tabix_args=TABIX_ARGS,
                       mlog10p_col="mlog10p_burden", per_trait_dir=args.per_trait_dir)

    if args.stage:
        stage(args.output_dir, args.per_trait_dir, traits)

    print("Done.")


if __name__ == "__main__":
    main()
