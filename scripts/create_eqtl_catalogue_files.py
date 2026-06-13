import warnings
import sys
import os
import glob

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from scipy.special import log_ndtr
import numpy as np
import polars as pl
from credible_set_stats import (
    calculate_stats,
    write_stats_json,
    get_tsv_header,
    stats_to_tsv_row,
)

# eQTL Catalogue R8 pilot munging.
# R8 SUSIE credible sets are distributed as parquet (R7 were tsv.gz) and every record
# carries an Ensembl gene_id, so the gene name is mapped directly from gene_id via the
# gene-level metadata (gene_counts) for ALL quantification methods -- including the new
# majiq splicing method and leafcutter, neither of which needs a per-quant phenotype
# metadata file. trait_original keeps the molecular_trait_id|quant_method form.


class NoDataException(Exception):
    pass


QUANT_DATA_TYPE = {
    "ge": "eQTL",
    "exon": "eQTL",
    "tx": "eQTL",
    "txrev": "eQTL",
    "microarray": "eQTL",
    "leafcutter": "sQTL",
    "majiq": "sQTL",
    "aptamer": "pQTL",
    "protein": "pQTL",
}

# columns needed from the credible-set file (parquet or tsv)
_NEEDED = [
    "molecular_trait_id",
    "gene_id",
    "cs_id",
    "variant",
    "cs_size",
    "pip",
    "pvalue",
    "beta",
    "se",
    "cs_min_r2",
]


def _read_cs(file_path: str) -> pl.DataFrame:
    if file_path.endswith(".parquet"):
        return pl.read_parquet(file_path, columns=_NEEDED)
    return pl.read_csv(
        file_path,
        separator="\t",
        columns=_NEEDED,
        schema_overrides={
            "molecular_trait_id": pl.Utf8,
            "gene_id": pl.Utf8,
            "cs_id": pl.Utf8,
            "variant": pl.Utf8,
            "cs_size": pl.Int32,
            "pip": pl.Float64,
            "pvalue": pl.Float64,
            "beta": pl.Float64,
            "se": pl.Float64,
            "cs_min_r2": pl.Float64,
        },
        infer_schema_length=100000,
        null_values=["NA"],
    )


def merge(
    study: str,
    data_type: str,
    file_path: str,
    cell_type: str,
    quant_method: str,
    gene_lookup: pl.DataFrame,
    pheno_lookup: pl.DataFrame,
    variant_annotation: pl.DataFrame | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    df = _read_cs(file_path)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="divide by zero encountered in log10")
        merged_df = (
            # gene name: prefer the authoritative phenotype->gene metadata (needed for the
            # splicing quants leafcutter/majiq), fall back to the Ensembl gene_id mapping
            df.with_columns(
                pl.col("gene_id").str.replace(r"\..*$", "").alias("_gene_key")
            )
            .join(gene_lookup, left_on="_gene_key", right_on="gene_key", how="left")
            .join(pheno_lookup, left_on="molecular_trait_id", right_on="phenotype_id", how="left")
            .with_columns(
                pl.lit(study).alias("dataset"),
                pl.lit(data_type).alias("data_type"),
                pl.lit(cell_type).alias("cell_type"),
                pl.col("variant").str.split("_").list.get(0).str.replace("chr", "").str.replace("X", "23").cast(pl.Int8).alias("chr"),
                pl.col("variant").str.split("_").list.get(1).cast(pl.Int32).alias("pos"),
                pl.col("variant").str.split("_").list.get(2).alias("ref"),
                pl.col("variant").str.split("_").list.get(3).alias("alt"),
                pl.coalesce(pl.col("gene_name_pheno"), pl.col("gene_name")).alias("trait"),
                pl.concat_str([pl.col("molecular_trait_id"), pl.lit(quant_method)], separator="|").alias("trait_original"),
                # compute -log10(p) from beta and se when p-value underflows
                pl.when(pl.col("pvalue").eq(0))
                .then(((-log_ndtr(-(pl.col("beta") / pl.col("se")).abs()) - np.log(2)) / np.log(10)).round(4))
                .otherwise((-np.log10(pl.col("pvalue")).round(4)))
                .alias("mlog10p"),
                # scientific notation strings for beta/se to avoid rounding se to zero
                pl.col("beta").map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8).alias("beta"),
                pl.col("se").map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8).alias("se"),
                pl.col("pip").round(4),
                pl.lit("NA").alias("aaf"),  # filled from variant annotation when available
                pl.col("cs_min_r2").round(4),
            )
            .select(
                "dataset", "data_type", "trait", "trait_original", "cell_type",
                "chr", "pos", "ref", "alt", "mlog10p", "beta", "se", "pip",
                "cs_id", "cs_size", "cs_min_r2", "aaf",
            )
        )

    null_traits = merged_df.filter(pl.col("trait").is_null())

    if variant_annotation is not None:
        merged_df = merged_df.drop(
            [c for c in merged_df.columns if c.endswith("most_severe") or c.startswith("aaf")]
        )
        merged_df = (
            merged_df.with_columns(
                pl.concat_str(
                    [pl.col("chr"), pl.col("pos").cast(pl.Utf8), pl.col("ref"), pl.col("alt")],
                    separator=":",
                ).alias("variant_id")
            )
            .join(variant_annotation, on="variant_id", how="left")
            .drop("variant_id")
        )

    return (
        merged_df.sort("chr", "pos", "ref", "alt", "trait"),
        null_traits.sort("chr", "pos", "ref", "alt", "trait"),
    )


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python create_eqtl_catalogue_files.py <data_dir> <variant_annotation_file>")
        sys.exit(1)
    data_dir = sys.argv[1]
    variant_annotation_file = sys.argv[2]

    metadata = pl.read_csv("metadata/eqtl_catalogue_files.tsv", separator="\t", has_header=False)
    study_metadata = pl.read_csv("metadata/eqtl_catalogue_studies.tsv", separator="\t")

    # gene_id -> gene_name lookup (Ensembl gene metadata); symbol-less genes map to their gene_id
    gene_lookup = (
        pl.read_csv(
            f"{data_dir}/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz",
            separator="\t",
            schema_overrides={"chromosome": pl.Utf8},
        )
        .select(
            pl.col("phenotype_id").str.replace(r"\..*$", "").alias("gene_key"),
            pl.col("gene_name"),
        )
        .unique(subset="gene_key")
    )

    # authoritative phenotype_id -> gene_name lookup from per-quant phenotype metadata
    # (e.g. leafcutter/majiq splicing-feature metadata); gene_counts excluded as it is the
    # gene_id fallback above. Used as the primary gene-name source in merge().
    pheno_files = [
        f for f in glob.glob(f"{data_dir}/*_phenotype_metadata.tsv.gz")
        if not os.path.basename(f).startswith("gene_counts_")
    ]
    if pheno_files:
        pheno_lookup = (
            pl.concat([
                pl.read_csv(f, separator="\t", schema_overrides={"chromosome": pl.Utf8})
                .select("phenotype_id", pl.col("gene_name").alias("gene_name_pheno"))
                for f in pheno_files
            ])
            .unique(subset="phenotype_id")
        )
    else:
        pheno_lookup = pl.DataFrame(schema={"phenotype_id": pl.Utf8, "gene_name_pheno": pl.Utf8})

    if variant_annotation_file is not None:
        variant_annotation = (
            pl.scan_csv(variant_annotation_file, separator="\t")
            .rename({"#variant": "variant_id"})
            .with_columns(pl.col("AF").map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8).alias("aaf"))
            .select("variant_id", "aaf", "most_severe", "gene_most_severe")
            .collect()
        )
    else:
        variant_annotation = None

    studies = metadata["column_1"].to_list()
    files = metadata["column_2"].to_list()
    out_dir = f"{data_dir}/eqtlcat_per_study"

    all_stats = []
    for i, study in enumerate(studies):
        out_file = f"{out_dir}/{study}.SUSIE.munged.tsv"
        if os.path.exists(out_file):
            print(f"{i+1}/{len(studies)}: {study}: already done, skipping")
            continue

        row = study_metadata.filter(pl.col("dataset_id").eq(study))
        tissue_label = row["tissue_label"][0]
        condition_label = row["condition_label"][0]
        cell_type = f"{tissue_label}|{condition_label}".replace(" ", "_")
        quant_method = row["quant_method"][0]
        if quant_method not in QUANT_DATA_TYPE:
            raise ValueError(f"Unknown quant_method for {study}: {quant_method}")
        data_type = QUANT_DATA_TYPE[quant_method]

        try:
            merged_df, null_traits = merge(
                study, data_type, files[i], cell_type, quant_method, gene_lookup, pheno_lookup, variant_annotation
            )
            if len(null_traits) > 0:
                null_traits.write_csv(
                    f"{out_dir}/{study}.SUSIE.munged.null_traits.tsv", separator="\t", null_value="NA"
                )
                print(
                    f"{i+1}/{len(studies)}: {study}: {len(set(null_traits['trait'].to_list()))} unmapped traits"
                )

            merged_df.write_csv(out_file, separator="\t", null_value="NA")

            # partition once (O(n)) instead of filtering per trait (O(traits*n))
            for trait_df in merged_df.filter(pl.col("trait").is_not_null()).partition_by("trait"):
                trait = trait_df["trait"][0]
                stats = calculate_stats(trait_df)
                write_stats_json(stats, f"{out_dir}/{study}.{trait}.SUSIE.munged.stats.json")
                all_stats.append(stats)

            print(f"{i+1}/{len(studies)}: {study}: OK")
        except NoDataException as e:
            print(f"{i+1}/{len(studies)}: {study}: {e}")

    with open(f"{out_dir}/credible_set_stats.tsv", "w") as f:
        f.write(get_tsv_header() + "\n")
        for s in all_stats:
            f.write(stats_to_tsv_row(s) + "\n")
    print(f"Wrote aggregate stats for {len(all_stats)} traits")
