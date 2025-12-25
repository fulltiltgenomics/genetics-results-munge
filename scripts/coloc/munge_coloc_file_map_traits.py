import polars as pl
import os
import sys

### THIS SCRIPT ASSUMES THAT EQTL CATALOGUE IS ALWAYS TRAIT2
### AND THAT EQTL CATALOGUE HAS NOT BEEN COLOCALIZED WITH OTHER QTL DATA (THOSE OTHER QTL TRAIT NAMES ARE NOT MAPPED)

# TODO FinnLiver mapping, get mapping file for that
# TODO |ge etc. for FinnLiver

DATA_TYPES = pl.Enum(["eQTL", "caQTL", "sQTL", "pQTL", "GWAS", "metaboQTL"])


def munge_non_eqtl_cat_data(filename: str, run_name: str):
    trait_mapping_files = {
        "FinnGen_snRNAseq": "/mnt/disks/dist_data/features.gex.tsv.gz",
        "FinnGen_Olink": "/mnt/disks/dist_data/manual_gene_mapping.20250916.tsv",
        "FinnGen_SomaScan": "/mnt/disks/dist_data/soma_aptamer_gene_mapping.tsv",
        "UKB_PPP": "/mnt/disks/dist_data/manual_gene_mapping.20250916.tsv",
    }
    soma_extra_mapping = "/mnt/disks/dist_data/manual_gene_mapping.20250916.tsv"
    # somascan mapping file may have "wrong" gene names so remap those using the manual mapping file
    # TODO do this up front..
    if "FinnGen_SomaScan" in trait_mapping_files.keys():
        mapping = pl.read_csv(
            trait_mapping_files["FinnGen_SomaScan"], separator="\t", has_header=False
        ).rename({"column_1": "trait", "column_2": "mapped_trait"})
        remapping = pl.read_csv(
            soma_extra_mapping, separator="\t", has_header=False
        ).rename({"column_1": "mapped_trait", "column_2": "remapped_trait"})
        mapping.join(remapping, on="mapped_trait", how="left").with_columns(
            pl.when(~pl.col("remapped_trait").is_null())
            .then(pl.col("remapped_trait"))
            .otherwise("mapped_trait")
            .alias("mapped_trait")
        ).select(["trait", "mapped_trait"]).write_csv(
            "soma_mapping.tsv", separator="\t", include_header=False
        )
        trait_mapping_files["FinnGen_SomaScan"] = "soma_mapping.tsv"
    trait_mapping = None
    for dataset, trait_mapping_file in trait_mapping_files.items():
        mapping = (
            pl.read_csv(
                trait_mapping_file,
                separator="\t",
                has_header=False,
            )
            .rename({"column_1": "trait", "column_2": "mapped_trait"})
            .with_columns(
                pl.lit(dataset).alias("dataset"),
            )
            .select("dataset", "trait", "mapped_trait")
        )
        if trait_mapping is None:
            trait_mapping = mapping
        else:
            trait_mapping = trait_mapping.vstack(mapping)

    data = (
        pl.read_csv(
            filename,
            separator="\t",
            schema_overrides={
                "data_type1": DATA_TYPES,
                "data_type2": DATA_TYPES,
                "chr": pl.UInt8,
                "pos": pl.UInt32,
                "cs1_size": pl.UInt16,
                "cs2_size": pl.UInt16,
                "cs_overlap": pl.UInt16,
            },
            infer_schema_length=100000,
        )
        .filter(~(pl.col("dataset2").str.starts_with("QTD")))
        .join(
            trait_mapping,
            left_on=["#dataset1", "trait1"],
            right_on=["dataset", "trait"],
            how="left",
        )
        .rename(
            {
                "trait1": "trait1_original",
                "mapped_trait": "trait1",
            }
        )
        .with_columns(
            pl.when(pl.col("trait1").is_null())
            .then(pl.col("trait1_original"))
            .otherwise(pl.col("trait1"))
            .alias("trait1"),
        )
        .join(
            trait_mapping,
            left_on=["dataset2", "trait2"],
            right_on=["dataset", "trait"],
            how="left",
        )
        .rename(
            {
                "trait2": "trait2_original",
                "mapped_trait": "trait2",
            }
        )
        .with_columns(
            pl.when(pl.col("trait2").is_null())
            .then(pl.col("trait2_original"))
            .otherwise(pl.col("trait2"))
            .alias("trait2"),
        )
        .with_columns(pl.col("trait1").str.replace_all("\|", "_").alias("trait1"))
        .with_columns(pl.col("trait2").str.replace_all("\|", "_").alias("trait2"))
    )

    data = data.select(
        "#dataset1",
        "dataset2",
        "data_type1",
        "data_type2",
        "trait1",
        "trait1_original",
        "trait2",
        "trait2_original",
        *[
            col
            for col in data.columns
            if col
            not in [
                "#dataset1",
                "dataset2",
                "data_type1",
                "data_type2",
                "trait1",
                "trait1_original",
                "trait2",
                "trait2_original",
            ]
        ],
    )

    data.write_csv(f"{run_name}.non_eqtl_cat_data.tsv", separator="\t", null_value="NA")

    # write individual trait data
    traits1 = (
        data.select("#dataset1", "trait1")
        .with_columns(
            pl.concat_str(pl.col("#dataset1"), pl.col("trait1"), separator="|").alias(
                "dataset_trait"
            )
        )
        .select("dataset_trait")
        .unique()
        .to_series()
        .to_list()
    )
    traits2 = (
        data.select("dataset2", "trait2")
        .with_columns(
            pl.concat_str(pl.col("dataset2"), pl.col("trait2"), separator="|").alias(
                "dataset_trait"
            )
        )
        .select("dataset_trait")
        .unique()
        .to_series()
        .to_list()
    )
    all_traits = sorted(list(set(traits1 + traits2)))
    os.makedirs(f"{run_name}", exist_ok=True)
    for dataset_trait in all_traits:
        dataset, trait = dataset_trait.split("|")
        os.makedirs(f"{run_name}/{dataset}", exist_ok=True)
        data.filter(
            (pl.col("data_type1") == "GWAS")
            & (
                ((pl.col("#dataset1") == dataset) & (pl.col("trait1") == trait))
                | ((pl.col("dataset2") == dataset) & (pl.col("trait2") == trait))
            )
        ).write_csv(
            f"{run_name}/{dataset}/{trait.replace('/', '_')}.tsv",  # replace / with _ to avoid directories (there is "THRA1/BTR" gene in features.gex.tsv.gz)
            separator="\t",
            null_value="NA",
        )
        print(f"{dataset} {trait} OK")


def munge_eqtl_cat_data(filename: str, eqtl_cat_metadata_filename: str, run_name: str):
    eqtl_cat_study_metadata = (
        pl.read_csv(
            eqtl_cat_metadata_filename,
            separator="\t",
        )
        .with_columns(
            pl.col("sample_group")
            .str.replace("+", "_", literal=True)
            .alias("sample_group")
        )
        .with_columns(
            pl.concat_str(
                pl.col("tissue_label"),
                pl.col("condition_label"),
                separator="|",
            )
            .str.replace_all(" ", "_")
            .alias("cell_type"),
            pl.when(pl.col("quant_method").str.contains("leafcutter"))
            .then(pl.lit("sQTL"))
            .when(pl.col("quant_method").str.contains("aptamer"))
            .then(pl.lit("pQTL"))
            .otherwise(pl.lit("eQTL"))
            .alias("data_type"),
        )
        .select("dataset_id", "cell_type", "quant_method", "data_type")
    )
    eqtl_cat_study_metadata_dicts = eqtl_cat_study_metadata.to_dict(as_series=False)
    eqtl_cat_quant_method_mapping = dict(
        zip(
            eqtl_cat_study_metadata_dicts["dataset_id"],
            eqtl_cat_study_metadata_dicts["quant_method"],
        )
    )

    # create a mapping from qtd study id to dataframe of molecular trait mapping
    eqtl_cat_trait_mapping_files = {
        "exon": [
            "/mnt/disks/data/eqtl_catalogue_r7/metadata/exon_counts_Ensembl_105_phenotype_metadata.tsv.gz"
        ],
        "ge": [
            "/mnt/disks/data/eqtl_catalogue_r7/metadata/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz"
        ],
        "microarray": [
            "/mnt/disks/data/eqtl_catalogue_r7/metadata/Affy_Human_Gene_1_0_ST_Ensembl_96_phenotype_metadata.tsv.gz",  # Affy
            "/mnt/disks/data/eqtl_catalogue_r7/metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz",  # Illumina
        ],
        "aptamer": [
            "/mnt/disks/data/eqtl_catalogue_r7/metadata/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz"
        ],
        "tx": [
            "/mnt/disks/data/eqtl_catalogue_r7/metadata/transcript_usage_Ensembl_105_phenotype_metadata.tsv.gz"
        ],
        "txrev": [
            "/mnt/disks/data/eqtl_catalogue_r7/metadata/txrevise_Ensembl_105_phenotype_metadata.tsv.gz"
        ],
        "leafcutter": [
            "/mnt/disks/data/eqtl_catalogue_r7/metadata/leafcutter_[STUDY]_Ensembl_105_phenotype_metadata.tsv.gz"
        ],
    }
    eqtl_cat_study_trait_mapping = {}
    for quant_method, trait_mapping_file_paths in eqtl_cat_trait_mapping_files.items():
        study_ids = eqtl_cat_study_metadata.filter(
            pl.col("quant_method").eq(quant_method)
        )["dataset_id"].to_list()
        # separate mapping files for each leafcutter study
        if quant_method == "leafcutter":
            for study in study_ids:
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
                eqtl_cat_study_trait_mapping[study] = trait_df
        else:
            trait_df = pl.concat(
                [
                    pl.read_csv(
                        trait_mapping_file_path,
                        separator="\t",
                        schema_overrides={
                            "chromosome": pl.Utf8,
                        },
                    ).select("phenotype_id", "gene_name")
                    for trait_mapping_file_path in trait_mapping_file_paths
                ]
            )
            for study in study_ids:
                eqtl_cat_study_trait_mapping[study] = trait_df

    data = pl.read_csv(
        filename,
        separator="\t",
        schema_overrides={
            "data_type1": DATA_TYPES,
            "data_type2": DATA_TYPES,
            "chr": pl.UInt8,
            "pos": pl.UInt32,
            "cs1_size": pl.UInt16,
            "cs2_size": pl.UInt16,
            "cs_overlap": pl.UInt16,
        },
        infer_schema_length=100000,
    ).filter(
        (pl.col("dataset2").str.starts_with("QTD"))
    )  # filter to eQTL Catalogue datasets
    if len(data) == 0:
        print("No eQTL Catalogue data found")
        return

    # map eQTL Catalogue study by study due to leafcutter trait mapping per study
    eqtl_cat_studies = data.select(pl.col("dataset2")).unique().to_series().to_list()
    eqtl_cat_data = None
    os.makedirs(f"{run_name}/eQTL_Catalogue", exist_ok=True)
    for i, study in enumerate(eqtl_cat_studies):
        quant_method = eqtl_cat_quant_method_mapping[study]
        study_data = (
            data.filter(pl.col("dataset2") == study)
            .join(
                eqtl_cat_study_trait_mapping[study],
                left_on="trait2",
                right_on="phenotype_id",
                how="left",
            )
            .with_columns(
                pl.col("trait1").alias("trait1_original"),
            )
            .with_columns(
                pl.concat_str(
                    [pl.col("trait2"), pl.lit(quant_method)], separator="|"
                ).alias("trait2_original"),
            )
            .with_columns(
                pl.col("gene_name").alias("trait2"),
            )
        )
        study_data = study_data.select(
            "#dataset1",
            "dataset2",
            "data_type1",
            "data_type2",
            "trait1",
            "trait1_original",
            "trait2",
            "trait2_original",
            *[
                col
                for col in study_data.columns
                if col
                not in [
                    "#dataset1",
                    "dataset2",
                    "data_type1",
                    "data_type2",
                    "trait1",
                    "trait1_original",
                    "trait2",
                    "trait2_original",
                    "gene_name",
                ]
            ],
        )
        study_data.write_csv(
            f"{run_name}/eQTL_Catalogue/{study}.tsv",
            separator="\t",
            null_value="NA",
        )
        if eqtl_cat_data is None:
            eqtl_cat_data = study_data
        else:
            eqtl_cat_data = eqtl_cat_data.vstack(study_data)
        print(f"{i+1}/{len(eqtl_cat_studies)}: {study} OK")

    eqtl_cat_data.sort(
        "chr", "region_start_min", "region_end_max", "trait1", "trait2"
    ).write_csv(f"{run_name}.eqtl_cat_data.tsv", separator="\t", null_value="NA")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python3 munge_coloc_file_map_traits.py <filename> <eqtl_cat_metadata_filename> <run_name>"
        )
        sys.exit(1)
    filename = sys.argv[1]
    eqtl_cat_metadata_filename = sys.argv[2]
    run_name = sys.argv[3]
    munge_non_eqtl_cat_data(filename, run_name)
    munge_eqtl_cat_data(filename, eqtl_cat_metadata_filename, run_name)
