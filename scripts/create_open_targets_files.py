"""
Convert an Open Targets credible set release to the credible set TSV format used by the rest of
this repository.

Input
  <data_dir>/credible_set/*.parquet   credible set parquet files of an Open Targets release
  <variant_annotation_file>           FinnGen annotated variants, needs the columns
                                      #variant (chr:pos:ref:alt, X as 23), AF, most_severe,
                                      gene_most_severe

Output
  <data_dir>/<dataset>_cs_95.tsv      unsorted, with a header; the shell driver sorts, bgzips
                                      and tabix indexes it

Only non-FinnGen GWAS credible sets fine-mapped with SuSiE are kept, and of those only the
variants flagged as belonging to the 95 % credible set. The 99 % credible sets are not written
because the release does not consistently distinguish them from the 95 % ones.

Each parquet file is converted on its own so that only one file's worth of exploded loci is held
at a time; the annotation is read once and reduced to the variants actually present in the
credible sets before the single join.
"""

import glob
import os
import sys

import numpy as np
import polars as pl

# variant ids spell the X chromosome out, everything downstream uses 23
CHR_MAP = {"X": "23"}

PARQUET_COLUMNS = [
    "studyLocusId",
    "studyId",
    "studyType",
    "finemappingMethod",
    "purityMinR2",
    "variantId",
    "pValueMantissa",
    "pValueExponent",
    "locus",
]

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


def _sci(col: str) -> pl.Expr:
    """Format a float column as %.3e in one numpy call rather than per row, keeping nulls null."""
    formatted = pl.col(col).fill_null(0.0).map_batches(
        lambda s: pl.Series(np.char.mod("%.3e", s.to_numpy())), return_dtype=pl.String
    )
    return pl.when(pl.col(col).is_null()).then(None).otherwise(formatted).alias(col)


def _mlog10p(mantissa: pl.Expr, exponent: pl.Expr) -> pl.Expr:
    return (-mantissa.cast(pl.Float64).log10() - exponent).round(4)


def convert_pq_to_df(parquet_path: str) -> pl.DataFrame:
    """Read one credible set parquet file into one row per 95 % credible set variant."""
    df = (
        pl.read_parquet(parquet_path, columns=PARQUET_COLUMNS)
        .filter(
            (pl.col("studyType") == "gwas")
            & ~pl.col("studyId").str.contains("FINNGEN", literal=True)
            & pl.col("finemappingMethod").str.to_lowercase().str.contains("susie", literal=True)
            & pl.col("locus").is_not_null()
        )
        .rename(
            {
                "variantId": "lead_variant_id",
                "pValueMantissa": "cs_p_mantissa",
                "pValueExponent": "cs_p_exponent",
            }
        )
        # the row index identifies the credible set a locus variant came from, which is what
        # cs_size counts over; studyLocusId is not relied on to be unique
        .with_row_index("_cs_row")
        .explode("locus")
        .unnest("locus")
        .filter(pl.col("is95CredibleSet"))
        .rename({"standardError": "se"})
    )

    df = df.with_columns(pl.col("variantId").str.split("_").alias("_cpra"))
    return df.with_columns(
        pl.col("studyId").alias("trait"),
        pl.col("_cpra").list.get(0).replace(CHR_MAP).cast(pl.Int32, strict=False).alias("chr"),
        pl.col("_cpra").list.get(1).cast(pl.Int64, strict=False).alias("pos"),
        pl.col("_cpra").list.get(2).alias("ref"),
        pl.col("_cpra").list.get(3).alias("alt"),
        # only the lead variant falls back to the credible set level p-value
        pl.when(pl.col("pValueMantissa").is_not_null() & pl.col("pValueExponent").is_not_null())
        .then(_mlog10p(pl.col("pValueMantissa"), pl.col("pValueExponent")))
        .when(pl.col("variantId") == pl.col("lead_variant_id"))
        .then(_mlog10p(pl.col("cs_p_mantissa"), pl.col("cs_p_exponent")))
        .alias("mlog10p"),
        _sci("beta"),
        _sci("se"),
        pl.col("posteriorProbability").round(4).alias("pip"),
        pl.col("studyLocusId").alias("cs_id"),
        pl.len().over("_cs_row").cast(pl.Int32).alias("cs_size"),
        pl.col("purityMinR2").round(4).alias("cs_min_r2"),
    ).select(
        "trait",
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
        pl.concat_str("chr", "pos", "ref", "alt", separator=":").alias("variant_id"),
    )


def read_annotation(annotation_path: str, variant_ids: pl.DataFrame) -> pl.DataFrame:
    """Read the variant annotation, keeping only the variants present in the credible sets."""
    anno = (
        pl.scan_csv(annotation_path, separator="\t", null_values=["NA"])
        .rename({"#variant": "variant_id"})
        .select("variant_id", "AF", "most_severe", "gene_most_severe")
        .join(variant_ids.lazy(), on="variant_id", how="semi")
        .collect()
    )
    return anno.with_columns(_sci("AF").alias("aaf")).drop("AF")


def main(dataset: str, data_dir: str, annotation_path: str) -> None:
    files = sorted(glob.glob(os.path.join(data_dir, "credible_set", "*.parquet")))
    if not files:
        sys.exit(f"no parquet files under {data_dir}/credible_set")

    frames = []
    for i, path in enumerate(files, 1):
        frames.append(convert_pq_to_df(path))
        print(f"[{i}/{len(files)}] {os.path.basename(path)}: {frames[-1].height} rows", flush=True)

    cs = pl.concat(frames)
    del frames
    print(
        f"{cs.height} credible set variants, "
        f"{cs['trait'].n_unique()} studies, "
        f"{cs['cs_id'].n_unique()} credible sets"
    )

    variant_ids = cs.select("variant_id").unique()
    anno = read_annotation(annotation_path, variant_ids)
    print(f"{anno.height} of {variant_ids.height} variants found in the annotation")

    output_path = os.path.join(data_dir, f"{dataset}_cs_95.tsv")
    cs.join(anno, on="variant_id", how="left").with_columns(
        pl.lit(dataset).alias("dataset"),
        pl.lit("GWAS").alias("data_type"),
        pl.col("trait").alias("trait_original"),
        pl.lit(None, dtype=pl.String).alias("cell_type"),
    ).select(OUTPUT_COLUMNS).write_csv(output_path, separator="\t", null_value="NA")
    print(f"wrote {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(
            "usage: python create_open_targets_files.py "
            "<dataset_name> <data_dir> <variant_annotation_file>"
        )
    main(sys.argv[1], sys.argv[2], sys.argv[3])
