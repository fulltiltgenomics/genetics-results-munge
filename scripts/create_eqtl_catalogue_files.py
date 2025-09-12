import warnings
import sys
from scipy.special import log_ndtr
from typing import Literal
import numpy as np
import polars as pl

Datatype = Literal["eQTL", "sQTL", "pQTL"]


class NoDataException(Exception):
    pass


def merge(
    study: str,
    data_type: Datatype,
    file_path: str,
    cell_type: str,
    quant_method: str,
    trait_mapping_file_paths: list[str],
    variant_annotation: pl.DataFrame | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    df = pl.read_csv(
        file_path,
        separator="\t",
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
        infer_schema_length=100000,  # needs to be large because some float fields may not otherwise be detected correctly
        null_values=["NA"],
    )

    # here i read the same files over and over again but it saves ram
    trait_df = pl.concat(
        [
            pl.read_csv(
                trait_mapping_file_path.replace("[STUDY]", study),
                separator="\t",
                schema_overrides={
                    "chromosome": pl.Utf8,
                },
            ).select("phenotype_id", "gene_name")
            for trait_mapping_file_path in trait_mapping_file_paths
        ]
    )

    with warnings.catch_warnings():
        # zero division warnings are given although we handle them
        warnings.filterwarnings("ignore", message="divide by zero encountered in log10")
        merged_df = (
            df.join(
                trait_df,
                left_on="molecular_trait_id",
                right_on="phenotype_id",
                how="left",
            )
            .with_columns(
                pl.lit(study).alias("dataset"),
                pl.lit(data_type).alias("data_type"),
                pl.lit(cell_type).alias("cell_type"),
                pl.col("variant")
                .str.split("_")
                .list.get(0)
                .str.replace("chr", "")
                .str.replace("X", "23")
                .cast(pl.Int8)
                .alias("chr"),
                pl.col("variant")
                .str.split("_")
                .list.get(1)
                .cast(pl.Int32)
                .alias("pos"),
                pl.col("variant").str.split("_").list.get(2).alias("ref"),
                pl.col("variant").str.split("_").list.get(3).alias("alt"),
                pl.col("gene_name").alias("trait"),
                pl.concat_str(
                    [pl.col("molecular_trait_id"), pl.lit(quant_method)], separator="|"
                ).alias("trait_original"),
                # compute -log10(p) from beta and se when p-value underflows
                pl.when(pl.col("pvalue").eq(0))
                .then(
                    (
                        (-log_ndtr(-(pl.col("beta") / pl.col("se")).abs()) - np.log(2))
                        / np.log(10)
                    ).round(4)
                )
                .otherwise((-np.log10(pl.col("pvalue")).round(4)))
                .alias("mlog10p"),
                # use scientific notation strings for beta and se to avoid rounding them (se) to zero
                pl.col("beta")
                .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
                .alias("beta"),
                pl.col("se")
                .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
                .alias("se"),
                pl.col("pip").round(4),
                pl.col("cs_min_r2").round(4),
            )
            .select(
                pl.col("dataset"),
                pl.col("data_type"),
                pl.col("trait"),
                pl.col("trait_original"),
                pl.col("cell_type"),
                pl.col("chr"),
                pl.col("pos"),
                pl.col("ref"),
                pl.col("alt"),
                pl.col("mlog10p"),
                pl.col("beta"),
                pl.col("se"),
                pl.col("pip"),
                pl.col("cs_id"),
                pl.col("cs_size"),
                pl.col("cs_min_r2"),
            )
        )

    null_traits = merged_df.filter(pl.col("trait").is_null())

    if variant_annotation is not None:
        # drop existing annotation columns if any
        merged_df = merged_df.drop(
            [col for col in merged_df.columns if col.endswith("most_severe")]
        )
        merged_df = (
            merged_df.with_columns(
                pl.concat_str(
                    [
                        pl.col("chr"),
                        pl.col("pos").cast(pl.Utf8),
                        pl.col("ref"),
                        pl.col("alt"),
                    ],
                    separator=":",
                ).alias("variant_id")
            )
            .join(
                variant_annotation,
                on="variant_id",
                how="left",
            )
            .drop("variant_id")
        )

    return (
        merged_df.sort("chr", "pos", "ref", "alt", "trait"),
        null_traits.sort("chr", "pos", "ref", "alt", "trait"),
    )


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(
            "Usage: python create_open_targets_files.py <data_dir> <variant_annotation_file>"
        )
        sys.exit(1)
    data_dir = sys.argv[1]
    variant_annotation_file = sys.argv[2]

    metadata = pl.read_csv(
        "metadata/eqtl_catalogue_files.tsv",
        separator="\t",
        has_header=False,
    )

    study_metadata = pl.read_csv(
        "metadata/eqtl_catalogue_studies.tsv",
        separator="\t",
    )

    trait_mapping_files = {
        "exon": [f"{data_dir}/exon_counts_Ensembl_105_phenotype_metadata.tsv.gz"],
        "ge": [f"{data_dir}/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz"],
        "microarray": [
            f"{data_dir}/Affy_Human_Gene_1_0_ST_Ensembl_96_phenotype_metadata.tsv.gz",  # Affy
            f"{data_dir}/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz",  # Illumina
        ],
        "aptamer": [f"{data_dir}/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz"],
        "tx": [f"{data_dir}/transcript_usage_Ensembl_105_phenotype_metadata.tsv.gz"],
        "txrev": [f"{data_dir}/txrevise_Ensembl_105_phenotype_metadata.tsv.gz"],
        "leafcutter": [
            f"{data_dir}/leafcutter_[STUDY]_Ensembl_105_phenotype_metadata.tsv.gz"
        ],
    }

    studies = metadata["column_1"].to_list()
    files = metadata["column_2"].to_list()

    if variant_annotation_file is not None:
        variant_annotation = (
            pl.scan_csv(variant_annotation_file, separator="\t")
            .rename({"#variant": "variant_id"})
            .select("variant_id", "most_severe", "gene_most_severe")
            .collect()
        )
    else:
        variant_annotation = None

    for i, study in enumerate(studies):
        tissue_label = study_metadata.filter(pl.col("dataset_id").eq(study))[
            "tissue_label"
        ][0]
        condition_label = study_metadata.filter(pl.col("dataset_id").eq(study))[
            "condition_label"
        ][0]
        cell_type = f"{tissue_label}|{condition_label}".replace(" ", "_")
        quant_method = study_metadata.filter(pl.col("dataset_id").eq(study))[
            "quant_method"
        ][0]
        if quant_method in ["exon", "ge", "tx", "txrev", "microarray"]:
            data_type = "eQTL"
        elif quant_method == "leafcutter":
            data_type = "sQTL"
        elif quant_method == "aptamer":
            data_type = "pQTL"
        else:
            raise ValueError(f"Unknown quant_method for {study}: {quant_method}")

        trait_mapping_file = trait_mapping_files[quant_method]

        try:
            merged_df, null_traits = merge(
                study,
                data_type,
                files[i],
                cell_type,
                quant_method,
                trait_mapping_file,
                variant_annotation,
            )
            if len(null_traits) > 0:
                null_traits.write_csv(
                    f"{data_dir}/eqtlcat_per_study/{study}.SUSIE.munged.null_traits.tsv",
                    separator="\t",
                    null_value="NA",
                )
                print(
                    f"{i+1}/{len(studies)}: {study}: {len(set(null_traits.select('trait').to_series().to_list()))} unmapped traits"
                )

            merged_df.write_csv(
                f"{data_dir}/eqtlcat_per_study/{study}.SUSIE.munged.tsv",
                separator="\t",
                null_value="NA",
            )
            print(f"{i+1}/{len(studies)}: {study}: OK")
        except NoDataException as e:
            print(f"{i+1}/{len(studies)}: {study}: {e}")
