import polars as pl
from scipy.special import log_ndtr
import numpy as np
import sys
import warnings

# time for file in colocQC.tsv.gz FinnGen-R12.eQTL.colocQC.tsv.gz FinnGen-KANTA.eQTL.colocQC.tsv.gz; do time scripts/munge_coloc_file.sh /mnt/disks/dist_data/coloc/$file /mnt/disks/dist_data/eqtl_catalogue_r7/dataset_metadata.tsv /mnt/disks/dist_data/coloc/$(basename $file .tsv.gz).munged.tsv; done

# convert dataset columns in data file to these names
# TODO include eQTL Catalogue in this same mapping
dataset_mapping = {
    "FIN-R12-Metabolism--Plasma-Metabolism": "FinnGen_NMR",
    "FIN-R12-Somascan--Plasma-Proteomics": "FinnGen_SomaScan",
    "FinnGen-Olink-pQTL-meta--Plasma-Proteomics": "FinnGen_Olink",
    "FinnGen-R13--GWAS": "FinnGen_R13",
    "FinnLiver_exon_Liver--FinnLiver": "FinnLiver",
    "FinnLiver_ge_Liver--FinnLiver": "FinnLiver",
    "FinnLiver_leafcutter_Liver--FinnLiver": "FinnLiver",
    "FinnLiver_tx_Liver--FinnLiver": "FinnLiver",
    "FinnLiver_txrev_Liver--FinnLiver": "FinnLiver",
    "GeneRisk--GWAS": "GeneRisk",
    "INTERVAL--Plasma-Proteomics": "INTERVAL",
    "Kanta_labs--GWAS": "FinnGen_kanta",
    "UKB-PPP--Plasma-Proteomics": "UKB_PPP",
    "UKB-finucane--GWAS": "UKB_Finucane",
    "FinnGen-R12--GWAS": "FinnGen_R12",
    "FinnGen-KANTA--GWAS": "FinnGen_kanta",
}

data_type_mapping = {
    "FIN-R12-Metabolism--Plasma-Metabolism": "metaboQTL",
    "FIN-R12-Somascan--Plasma-Proteomics": "pQTL",
    "FinnGen-Olink-pQTL-meta--Plasma-Proteomics": "pQTL",
    "FinnGen-R13--GWAS": "GWAS",
    "FinnLiver_exon_Liver--FinnLiver": "eQTL",
    "FinnLiver_ge_Liver--FinnLiver": "eQTL",
    "FinnLiver_leafcutter_Liver--FinnLiver": "sQTL",
    "FinnLiver_tx_Liver--FinnLiver": "eQTL",
    "FinnLiver_txrev_Liver--FinnLiver": "eQTL",
    "GeneRisk--GWAS": "GWAS",
    "INTERVAL--Plasma-Proteomics": "pQTL",
    "Kanta_labs--GWAS": "GWAS",
    "UKB-PPP--Plasma-Proteomics": "pQTL",
    "UKB-finucane--GWAS": "GWAS",
    "FinnGen-R12--GWAS": "GWAS",
    "FinnGen-KANTA--GWAS": "GWAS",
}

cell_type_mapping = {
    "FIN-R12-Metabolism--Plasma-Metabolism": "plasma",
    "FIN-R12-Somascan--Plasma-Proteomics": "plasma",
    "FinnGen-Olink-pQTL-meta--Plasma-Proteomics": "plasma",
    "FinnGen-R13--GWAS": None,
    "FinnLiver_exon_Liver--FinnLiver": "liver",
    "FinnLiver_ge_Liver--FinnLiver": "liver",
    "FinnLiver_leafcutter_Liver--FinnLiver": "liver",
    "FinnLiver_tx_Liver--FinnLiver": "liver",
    "FinnLiver_txrev_Liver--FinnLiver": "liver",
    "GeneRisk--GWAS": None,
    "INTERVAL--Plasma-Proteomics": "plasma",
    "Kanta_labs--GWAS": None,
    "UKB-PPP--Plasma-Proteomics": "plasma",
    "UKB-finucane--GWAS": None,
    "FinnGen-R12--GWAS": None,
    "FinnGen-KANTA--GWAS": None,
}

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


def main(filename, eqtl_cat_metadata_filename, output_filename):

    eqtl_cat_metadata = (
        pl.read_csv(eqtl_cat_metadata_filename, separator="\t")
        .with_columns(
            pl.col("sample_group")
            .str.replace("+", "_", literal=True)
            .alias("sample_group")
        )
        .with_columns(
            pl.concat_str(
                pl.col("study_label"),
                pl.col("sample_group"),
                pl.col("quant_method"),
                pl.lit("eQTL_Catalogue"),
                separator="--",
            ).alias("dataset"),
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
        .select("dataset_id", "dataset", "cell_type", "quant_method", "data_type")
        .to_dict(as_series=False)
    )
    eqtl_cat_dataset_mapping = dict(
        zip(eqtl_cat_metadata["dataset"], eqtl_cat_metadata["dataset_id"])
    )
    eqtl_cat_quant_method_mapping = dict(
        zip(eqtl_cat_metadata["dataset"], eqtl_cat_metadata["quant_method"])
    )
    eqtl_cat_cell_type_mapping = dict(
        zip(eqtl_cat_metadata["dataset"], eqtl_cat_metadata["cell_type"])
    )
    eqtl_cat_data_type_mapping = dict(
        zip(eqtl_cat_metadata["dataset"], eqtl_cat_metadata["data_type"])
    )

    with warnings.catch_warnings():
        # zero division warnings are given although we handle them
        warnings.filterwarnings("ignore", message="divide by zero encountered in log10")
        data = (
            pl.scan_csv(
                filename,
                separator="\t",
                infer_schema_length=100000,
                schema_overrides={
                    "trait1": pl.Utf8,
                    "trait2": pl.Utf8,
                    "region1": pl.Utf8,
                    "region2": pl.Utf8,
                    "cs1_log10bf": pl.Float64,
                    "cs2_log10bf": pl.Float64,
                    "clpp": pl.Float64,
                    "clpa": pl.Float64,
                    "probmass_1": pl.Float64,
                    "probmass_2": pl.Float64,
                    "nsnps": pl.Int32,
                    "nsnps1": pl.Int32,
                    "nsnps2": pl.Int32,
                    "cs1_size": pl.Int32,
                    "cs2_size": pl.Int32,
                    "cs_overlap": pl.Int32,
                    "hit1_info": pl.Utf8,
                    "hit2_info": pl.Utf8,
                },
                null_values=["NA"],
            )
            # standard FinnGen filtering
            .filter((pl.col("low_purity1") == 0) & (pl.col("low_purity2") == 0))
            .filter((pl.col("probmass_1") > 0.9) & (pl.col("probmass_2") > 0.9))
            .filter((pl.col("cs1_log10bf") > 0.9) & (pl.col("cs2_log10bf") > 0.9))
            # filter out Plasma-SingleCell because newer snRNAseq data are in separate files
            .filter(
                ~(pl.col("dataset1").str.ends_with("--Plasma-SingleCell"))
                & ~(pl.col("dataset2").str.ends_with("--Plasma-SingleCell"))
            )
            # filter out old Olink
            .filter(
                ~(pl.col("dataset1").str.ends_with("FIN-R12-Olink--Plasma-Proteomics"))
                & ~(
                    pl.col("dataset2").str.ends_with("FIN-R12-Olink--Plasma-Proteomics")
                )
            )
            .with_columns(
                pl.concat_str(pl.col("region1"), pl.col("cs1"), separator="_").alias(
                    "cs1_id"
                ),
                pl.concat_str(pl.col("region2"), pl.col("cs2"), separator="_").alias(
                    "cs2_id"
                ),
            )
            .with_columns(
                pl.col("region1")
                .str.split(":")
                .list.get(0)
                .str.replace("chr", "")
                .str.replace("X", "23")
                .cast(pl.Int8)
                .alias("chr"),
                # ugly regex because splitting by : and - is not enough as there can be extra - from negative region starts
                pl.col("region1")
                .str.replace(r"X", "23")
                .str.extract(r"^(chr\d+):(-?\d+)-(\d+)$", 2)
                .cast(pl.Int32)
                .alias("region1_start"),
                pl.col("region1")
                .str.replace(r"X", "23")
                .str.extract(r"^(chr\d+):(-?\d+)-(\d+)$", 3)
                .cast(pl.Int32)
                .alias("region1_end"),
                pl.col("region2")
                .str.replace(r"X", "23")
                .str.extract(r"^(chr\d+):(-?\d+)-(\d+)$", 2)
                .cast(pl.Int32)
                .alias("region2_start"),
                pl.col("region2")
                .str.replace(r"X", "23")
                .str.extract(r"^(chr\d+):(-?\d+)-(\d+)$", 3)
                .cast(pl.Int32)
                .alias("region2_end"),
            )
            .with_columns(
                pl.max_horizontal([pl.col("region1_start"), pl.lit(1)]).alias(
                    "region1_start"
                ),
                pl.max_horizontal([pl.col("region2_start"), pl.lit(1)]).alias(
                    "region2_start"
                ),
            )
            .with_columns(
                pl.min_horizontal(
                    [pl.col("region1_start"), pl.col("region2_start")]
                ).alias("region_start_min"),
                pl.max_horizontal([pl.col("region1_end"), pl.col("region2_end")]).alias(
                    "region_end_max"
                ),
            )
            .with_columns(
                pl.col("hit1")
                .str.replace_all("_", ":")
                .str.replace("chr", "")
                .str.replace("X", "23")
                .alias("hit1"),
                pl.col("hit2")
                .str.replace_all("_", ":")
                .str.replace("chr", "")
                .str.replace("X", "23")
                .alias("hit2"),
            )
            .with_columns(
                pl.when((pl.col("hit1_info").is_null()) | (pl.col("hit1_info") == "NA"))
                .then(pl.lit(None))
                .when(pl.col("hit1_info").str.split(",").list.len() >= 1)
                .then(
                    pl.col("hit1_info")
                    .str.split(",")
                    .list.get(
                        0, null_on_oob=True
                    )  # these list.get are evaluated regardless of condition and we need null_on_oob to avoid errors
                    .cast(pl.Float64, strict=False)
                )
                .otherwise(pl.lit(None))
                .alias("hit1_beta"),
                pl.when((pl.col("hit1_info").is_null()) | (pl.col("hit1_info") == "NA"))
                .then(pl.lit(None))
                .when(pl.col("hit1_info").str.split(",").list.len() >= 2)
                .then(
                    pl.col("hit1_info")
                    .str.split(",")
                    .list.get(1, null_on_oob=True)
                    .cast(pl.Float64, strict=False)
                )
                .otherwise(pl.lit(None))
                .alias("hit1_p"),
                pl.when((pl.col("hit2_info").is_null()) | (pl.col("hit2_info") == "NA"))
                .then(pl.lit(None))
                .when(pl.col("hit2_info").str.split(",").list.len() >= 1)
                .then(
                    pl.col("hit2_info")
                    .str.split(",")
                    .list.get(0, null_on_oob=True)
                    .cast(pl.Float64, strict=False)
                )
                .otherwise(pl.lit(None))
                .alias("hit2_beta"),
                pl.when((pl.col("hit2_info").is_null()) | (pl.col("hit2_info") == "NA"))
                .then(pl.lit(None))
                .when(pl.col("hit2_info").str.split(",").list.len() >= 2)
                .then(
                    pl.col("hit2_info")
                    .str.split(",")
                    .list.get(1, null_on_oob=True)
                    .cast(pl.Float64, strict=False)
                )
                .otherwise(pl.lit(None))
                .alias("hit2_p"),
            )
            .with_columns(  # TODO we can't currently handle underflowing p-values here because SE is missing
                pl.when((pl.col("hit1_p").is_null()))
                .then(None)
                .otherwise(
                    pl.when((pl.col("hit1_p").eq(0)))
                    .then(308)
                    .otherwise((-np.log10(pl.col("hit1_p")).round(4)))
                )
                .alias("hit1_mlog10p"),
                pl.when((pl.col("hit2_p").is_null()))
                .then(None)
                .otherwise(
                    pl.when((pl.col("hit2_p").eq(0)))
                    .then(308)
                    .otherwise((-np.log10(pl.col("hit2_p")).round(4)))
                )
                .alias("hit2_mlog10p"),
            )
            .with_columns(
                pl.col("PP.H0.abf")
                .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
                .alias("PP.H0.abf"),
                pl.col("PP.H1.abf")
                .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
                .alias("PP.H1.abf"),
                pl.col("PP.H2.abf")
                .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
                .alias("PP.H2.abf"),
                pl.col("PP.H3.abf")
                .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
                .alias("PP.H3.abf"),
                pl.col("PP.H4.abf")
                .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
                .alias("PP.H4.abf"),
                pl.when((pl.col("cs1_log10bf").is_infinite()))
                .then(1000)
                .otherwise(pl.col("cs1_log10bf").round(4))
                .alias("cs1_log10bf"),
                pl.when((pl.col("cs2_log10bf").is_infinite()))
                .then(1000)
                .otherwise(pl.col("cs2_log10bf").round(4))
                .alias("cs2_log10bf"),
                pl.col("clpp").map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8),
                pl.col("clpa").map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8),
                pl.col("probmass_1").map_elements(
                    lambda x: f"{x:.5e}", return_dtype=pl.Utf8
                ),
                pl.col("probmass_2").map_elements(
                    lambda x: f"{x:.5e}", return_dtype=pl.Utf8
                ),
            )
            .with_columns(
                pl.when(pl.col("dataset1").str.ends_with("--eQTL_Catalogue"))
                .then(
                    pl.col("dataset1").replace_strict(
                        eqtl_cat_data_type_mapping, default=None
                    )
                )
                .when(pl.col("dataset1").str.ends_with("--FinnGen_eQTL"))
                .then(pl.lit("eQTL"))
                .when(pl.col("dataset1").str.ends_with("--FinnGen_caQTL"))
                .then(pl.lit("caQTL"))
                .otherwise(
                    pl.col("dataset1").replace_strict(data_type_mapping, default=None)
                )
                .alias("data_type1")
            )
            .with_columns(
                pl.when(pl.col("dataset2").str.ends_with("--eQTL_Catalogue"))
                .then(
                    pl.col("dataset2").replace_strict(
                        eqtl_cat_data_type_mapping, default=None
                    )
                )
                .when(pl.col("dataset2").str.ends_with("--FinnGen_eQTL"))
                .then(pl.lit("eQTL"))
                .when(pl.col("dataset2").str.ends_with("--FinnGen_caQTL"))
                .then(pl.lit("caQTL"))
                .otherwise(
                    pl.col("dataset2").replace_strict(data_type_mapping, default=None)
                )
                .alias("data_type2")
            )
            .with_columns(
                pl.when(pl.col("dataset1").str.ends_with("--eQTL_Catalogue"))
                .then(
                    pl.col("dataset1").replace_strict(
                        eqtl_cat_cell_type_mapping, default=None
                    )
                )
                .when(pl.col("dataset1").str.ends_with("--FinnGen_eQTL"))
                .then(
                    pl.col("dataset1")
                    .str.replace(r"--FinnGen_eQTL$", "")
                    .str.replace(r"^.*predicted\.celltype\.", "")
                    .str.replace(r"\.chr\d+$|\.chrX$", "")
                    .str.replace(r"\.mean\.inv\.SAIGE", "")
                )
                .when(pl.col("dataset1").str.ends_with("--FinnGen_caQTL"))
                .then(
                    pl.col("dataset1")
                    .str.replace(r"--FinnGen_caQTL$", "")
                    .str.replace(r"^.*predicted\.celltype\.", "")
                    .str.replace(r"\.chr\d+$|\.chrX$", "")
                    .str.replace(r"\.mean\.inv\.SAIGE", "")
                    .str.replace(r"\.sum\.inv", "")
                )
                .otherwise(
                    pl.col("dataset1").replace_strict(cell_type_mapping, default=None)
                )
                .alias("cell_type1")
            )
            .with_columns(
                pl.when(pl.col("dataset2").str.ends_with("--eQTL_Catalogue"))
                .then(
                    pl.col("dataset2").replace_strict(
                        eqtl_cat_cell_type_mapping, default=None
                    )
                )
                .when(pl.col("dataset2").str.ends_with("--FinnGen_eQTL"))
                .then(
                    pl.col("dataset2")
                    .str.replace(r"--FinnGen_eQTL$", "")
                    .str.replace(r"^.*predicted\.celltype\.", "")
                    .str.replace(r"\.chr\d+$|\.chrX$", "")
                    .str.replace(r"\.mean\.inv\.SAIGE", "")
                )
                .when(pl.col("dataset2").str.ends_with("--FinnGen_caQTL"))
                .then(
                    pl.col("dataset2")
                    .str.replace(r"--FinnGen_caQTL$", "")
                    .str.replace(r"^.*predicted\.celltype\.", "")
                    .str.replace(r"\.chr\d+$|\.chrX$", "")
                    .str.replace(r"\.mean\.inv\.SAIGE", "")
                    .str.replace(r"\.sum\.inv", "")
                )
                .otherwise(
                    pl.col("dataset2").replace_strict(cell_type_mapping, default=None)
                )
                .alias("cell_type2")
            )
            .with_columns(
                pl.when(pl.col("dataset1").str.ends_with("--eQTL_Catalogue"))
                .then(
                    pl.col("dataset1").replace_strict(
                        eqtl_cat_dataset_mapping, default=None
                    )
                )
                .when(pl.col("dataset1").str.ends_with("--FinnGen_eQTL"))
                .then(pl.lit("FinnGen_snRNAseq"))
                .when(pl.col("dataset1").str.ends_with("--FinnGen_caQTL"))
                .then(pl.lit("FinnGen_ATACseq"))
                .otherwise(
                    pl.col("dataset1").replace_strict(dataset_mapping, default=None)
                )
                .alias("#dataset1")
            )
            .with_columns(
                pl.when(pl.col("dataset2").str.ends_with("--eQTL_Catalogue"))
                .then(
                    pl.col("dataset2").replace_strict(
                        eqtl_cat_dataset_mapping, default=None
                    )
                )
                .when(pl.col("dataset2").str.ends_with("--FinnGen_eQTL"))
                .then(pl.lit("FinnGen_snRNAseq"))
                .when(pl.col("dataset2").str.ends_with("--FinnGen_caQTL"))
                .then(pl.lit("FinnGen_ATACseq"))
                .otherwise(
                    pl.col("dataset2").replace_strict(dataset_mapping, default=None)
                )
                .alias("dataset2")
            )
        )

        data.select(
            [
                "#dataset1",
                "dataset2",
                "data_type1",
                "data_type2",
                "trait1",
                "trait2",
                "cell_type1",
                "cell_type2",
                "cs1_id",
                "cs2_id",
                "hit1",
                "hit2",
                "hit1_beta",
                "hit1_mlog10p",
                "hit2_beta",
                "hit2_mlog10p",
                "chr",
                "region_start_min",
                "region_end_max",
                "PP.H0.abf",
                "PP.H1.abf",
                "PP.H2.abf",
                "PP.H3.abf",
                "PP.H4.abf",
                "nsnps",
                "nsnps1",
                "nsnps2",
                "cs1_log10bf",
                "cs2_log10bf",
                "clpp",
                "clpa",
                "cs1_size",
                "cs2_size",
                "cs_overlap",
                "topInOverlap",
            ]
        ).sink_csv(output_filename, separator="\t", null_value="NA")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python munge_coloc_file.py <filename> <eqtl_cat_metadata_filename> <output_filename>"
        )
        sys.exit(1)
    filename = sys.argv[1]
    eqtl_cat_metadata_filename = sys.argv[2]
    output_filename = sys.argv[3]
    main(filename, eqtl_cat_metadata_filename, output_filename)
