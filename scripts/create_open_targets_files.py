"""
Convert an Open Targets credible set release to the credible set TSV format used by the rest of
this repository.

Input
  <data_dir>/credible_set/*.parquet    credible set parquet files of an Open Targets release
  <data_dir>/study_metadata/*.parquet  its study table, read for studyId and traitFromSource

Output
  <data_dir>/<dataset>_cs_95.tsv      unsorted, with a header; the shell driver sorts, bgzips
                                      and tabix indexes it

Only non-FinnGen GWAS credible sets fine-mapped with SuSiE are kept, and of those only the
variants flagged as belonging to the 95 % credible set. The 99 % credible sets are not written
because the release does not consistently distinguish them from the 95 % ones.

`trait` is the study's traitFromSource made whitespace-free (see `sanitize_trait_name`) followed
by `_(<studyId>)`, and `trait_original` is the bare studyId. A study's name is not unique - many
accessions share one - so the accession suffix is what makes `trait` identify a study, and
consumers cut exactly that suffix to get the name back.

`aaf`, `most_severe` and `gene_most_severe` are written as NA. The release has no allele
frequency, and consequence is stamped afterwards by `annotate_resource.sh`, so that an
annotation refresh needs no re-munge.

Each parquet file is converted on its own so that only one file's worth of exploded loci is held
at a time.
"""

import glob
import os
import sys
import unicodedata

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
        pl.col("studyId").alias("trait_original"),
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
        "trait_original",
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
    )


def sanitize_trait_name(name: str | None) -> str:
    """Make a study's trait name safe as one TSV field without changing how it reads.

    Whitespace runs become one underscore, control and format characters are dropped, and the
    result is NFC so that the same name always compares equal. Punctuation and non-ASCII letters
    are kept: the value is matched by exact string equality downstream, never used as a path.
    The double quote is the one exception - it is the quote character of the TSV readers and of
    the BigQuery load, which would take a field starting with it as a quoted field.
    """
    if not name:
        return ""
    kept = "".join(
        ch for ch in unicodedata.normalize("NFC", name)
        # tabs and newlines are control characters, so they go here rather than becoming
        # underscores below
        if not unicodedata.category(ch).startswith("C")
    )
    return "_".join(kept.split()).replace('"', "'")


def trait_label(study_id: str, name: str | None) -> str:
    """`<sanitized name>_(<studyId>)`, or the bare studyId when the study has no usable name."""
    sanitized = sanitize_trait_name(name)
    return f"{sanitized}_({study_id})" if sanitized else study_id


def read_trait_labels(data_dir: str, study_ids: pl.Series) -> pl.DataFrame:
    """Map each studyId present in the credible sets to its `trait` value."""
    files = sorted(glob.glob(os.path.join(data_dir, "study_metadata", "*.parquet")))
    if not files:
        sys.exit(f"no parquet files under {data_dir}/study_metadata")
    names = dict(
        pl.read_parquet(files, columns=["studyId", "traitFromSource"])
        .filter(pl.col("studyId").is_in(study_ids.implode()))
        .iter_rows()
    )
    missing = [s for s in study_ids if s not in names]
    if missing:
        print(f"{len(missing)} studies are not in the study table, e.g. {missing[:5]}")
    labels = pl.DataFrame(
        {
            "trait_original": study_ids,
            "trait": [trait_label(s, names.get(s)) for s in study_ids],
        }
    )
    # per-study files are named by trait_original while consumers select by trait, so the two
    # have to name the same study
    assert labels["trait_original"].n_unique() == labels.height, "duplicate studyId"
    assert labels["trait"].n_unique() == labels.height, "trait is not unique per studyId"
    assert not labels["trait"].str.contains(r"\s").any(), "whitespace left in trait"
    return labels


def main(dataset: str, data_dir: str) -> None:
    files = sorted(glob.glob(os.path.join(data_dir, "credible_set", "*.parquet")))
    if not files:
        sys.exit(f"no parquet files under {data_dir}/credible_set")

    frames = []
    for i, path in enumerate(files, 1):
        frames.append(convert_pq_to_df(path))
        print(f"[{i}/{len(files)}] {os.path.basename(path)}: {frames[-1].height} rows", flush=True)

    cs = pl.concat(frames)
    del frames
    labels = read_trait_labels(data_dir, cs["trait_original"].unique())
    print(
        f"{cs.height} credible set variants, "
        f"{labels.height} studies, "
        f"{cs['cs_id'].n_unique()} credible sets"
    )

    output_path = os.path.join(data_dir, f"{dataset}_cs_95.tsv")
    na = pl.lit(None, dtype=pl.String)
    cs.join(labels, on="trait_original", how="left").with_columns(
        pl.lit(dataset).alias("dataset"),
        pl.lit("GWAS").alias("data_type"),
        na.alias("cell_type"),
        na.alias("aaf"),
        na.alias("most_severe"),
        na.alias("gene_most_severe"),
    ).select(OUTPUT_COLUMNS).write_csv(output_path, separator="\t", null_value="NA")
    print(f"wrote {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: python create_open_targets_files.py <dataset_name> <data_dir>")
    main(sys.argv[1], sys.argv[2])
