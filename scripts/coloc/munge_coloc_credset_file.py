import polars as pl
from scipy.special import log_ndtr
import numpy as np
import warnings
import sys

warnings.filterwarnings(
    "ignore", category=RuntimeWarning, message="divide by zero encountered in log10"
)

# convert dataset column in data file to these names
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
    "UKB-PPP--Plasma-Proteomics": "UKB-PPP",
    "UKB-finucane--GWAS": "UKB-Finucane",
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
    eqtl_cat_cell_type_mapping = dict(
        zip(eqtl_cat_metadata["dataset"], eqtl_cat_metadata["cell_type"])
    )
    eqtl_cat_data_type_mapping = dict(
        zip(eqtl_cat_metadata["dataset"], eqtl_cat_metadata["data_type"])
    )
    eqtl_cat_quant_method_mapping = dict(
        zip(eqtl_cat_metadata["dataset"], eqtl_cat_metadata["quant_method"])
    )

    data = (
        pl.scan_csv(
            filename,
            separator="\t",
            infer_schema_length=100000,
            schema_overrides={
                "trait": pl.Utf8,
                "p": pl.Float64 | None,
                "beta": pl.Utf8,  # some data files are broken so first read beta as string and then filter non-numeric values out
                "se": pl.Utf8,
            },
            null_values=["NA", "T"],  # T allele...
        )
        .with_columns(
            pl.col("beta").cast(pl.Float64, strict=False),
            pl.col("se").cast(pl.Float64, strict=False),
        )
        .filter(pl.col("beta").is_not_null())
        .filter(pl.col("low_purity") == 0)
        # filter out Plasma-SingleCell because newer snRNAseq data are in separate files
        .filter(~(pl.col("dataset").str.ends_with("--Plasma-SingleCell")))
        # filter out old Olink
        .filter(~(pl.col("dataset").str.ends_with("FIN-R12-Olink--Plasma-Proteomics")))
        # TODO filter out SomaScan?
        .with_columns(
            pl.when((pl.col("p").eq(0)) | (pl.col("p").is_null()))
            .then(
                pl.when((pl.col("se").eq(0)) | (pl.col("se").is_null()))
                .then(None)
                .otherwise(
                    (
                        (-log_ndtr(-(pl.col("beta") / pl.col("se")).abs()) - np.log(2))
                        / np.log(10)
                    ).round(4)
                )
            )
            .otherwise((-np.log10(pl.col("p")).round(4)))
            .alias("mlog10p")
        )
        .with_columns(
            pl.concat_str(pl.col("region"), pl.col("cs"), separator="_").alias("cs_id")
        )
        .with_columns(
            pl.col("rsid")
            .str.split("_")
            .list.get(0)
            .str.replace("chr", "")
            .str.replace("X", "23")
            .cast(pl.Int8)
            .alias("chr"),
            pl.col("rsid").str.split("_").list.get(1).cast(pl.Int32).alias("pos"),
            pl.col("rsid").str.split("_").list.get(2).alias("ref"),
            pl.col("rsid").str.split("_").list.get(3).alias("alt"),
        )
        .with_columns(
            pl.col("beta")
            .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
            .alias("beta"),
            pl.col("se")
            .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
            .alias("se"),
            pl.col("cs_specific_prob").round(4),
        )
        .with_columns(
            pl.when(pl.col("dataset").str.ends_with("--eQTL_Catalogue"))
            .then(
                pl.col("dataset").replace_strict(eqtl_cat_dataset_mapping, default=None)
            )
            .when(pl.col("dataset").str.ends_with("--FinnGen_eQTL"))
            .then(pl.lit("FinnGen_snRNAseq"))
            .when(pl.col("dataset").str.ends_with("--FinnGen_caQTL"))
            .then(pl.lit("FinnGen_ATACseq"))
            .otherwise(pl.col("dataset").replace_strict(dataset_mapping, default=None))
            .alias("#dataset")
        )
        .with_columns(
            pl.when(pl.col("dataset").str.ends_with("--eQTL_Catalogue"))
            .then(
                pl.col("dataset").replace_strict(
                    eqtl_cat_data_type_mapping, default=None
                )
            )
            .when(pl.col("dataset").str.ends_with("--FinnGen_eQTL"))
            .then(pl.lit("eQTL"))
            .when(pl.col("dataset").str.ends_with("--FinnGen_caQTL"))
            .then(pl.lit("caQTL"))
            .otherwise(
                pl.col("dataset").replace_strict(data_type_mapping, default=None)
            )
            .alias("data_type")
        )
        .with_columns(
            pl.when(pl.col("dataset").str.ends_with("--eQTL_Catalogue"))
            .then(
                pl.col("dataset").replace_strict(
                    eqtl_cat_cell_type_mapping, default=None
                )
            )
            .when(pl.col("dataset").str.ends_with("--FinnGen_caQTL"))
            .then(
                pl.col("dataset")
                .str.replace(r"--FinnGen_caQTL$", "")
                .str.replace(r"^.*predicted\.celltype\.", "")
                .str.replace(r"\.chr\d+$|\.chrX$", "")
                .str.replace(r"\.mean\.inv\.SAIGE", "")
                .str.replace(r"\.sum\.inv", "")
            )
            .when(pl.col("dataset").str.ends_with("--FinnGen_eQTL"))
            .then(
                pl.col("dataset")
                .str.replace(r"--FinnGen_eQTL$", "")
                .str.replace(r"^.*predicted\.celltype\.", "")
                .str.replace(r"\.chr\d+$|\.chrX$", "")
                .str.replace(r"\.mean\.inv\.SAIGE", "")
            )
            .otherwise(
                pl.col("dataset").replace_strict(cell_type_mapping, default=None)
            )
            .alias("cell_type")
        )
        # TODO map to trait_original like we do in credible set munging
        .with_columns(
            pl.when(
                (pl.col("dataset").str.contains("_ge_"))
                | (pl.col("dataset").str.contains("--ge--"))
            )
            .then(pl.concat_str(pl.col("trait"), pl.lit("ge"), separator="|"))
            .when(
                (pl.col("dataset").str.contains("_tx_"))
                | (pl.col("dataset").str.contains("--tx--"))
            )
            .then(pl.concat_str(pl.col("trait"), pl.lit("tx"), separator="|"))
            .when(
                (pl.col("dataset").str.contains("_txrev_"))
                | (pl.col("dataset").str.contains("--txrev--"))
            )
            .then(pl.concat_str(pl.col("trait"), pl.lit("txrev"), separator="|"))
            .when(
                (pl.col("dataset").str.contains("_exon_"))
                | (pl.col("dataset").str.contains("--exon--"))
            )
            .then(pl.concat_str(pl.col("trait"), pl.lit("exon"), separator="|"))
            .when(
                (pl.col("dataset").str.contains("_leafcutter_"))
                | (pl.col("dataset").str.contains("--leafcutter--"))
            )
            .then(pl.concat_str(pl.col("trait"), pl.lit("leafcutter"), separator="|"))
            .when(pl.col("dataset").str.ends_with("--eQTL_Catalogue"))
            .then(
                pl.concat_str(
                    pl.col("trait"),
                    pl.col("dataset").replace_strict(
                        eqtl_cat_quant_method_mapping, default=None
                    ),
                    separator="|",
                )
            )
            .otherwise(pl.col("trait"))
            .alias("trait_original")
        )
        .rename({"cs_specific_prob": "pip"})
    )

    # if data.select("#dataset").is_null().any():
    #     raise ValueError(f"Could not map dataset in {filename} for some rows")

    # if data.select("data_type").is_null().any():
    #     raise ValueError(f"Could not map data type in {filename} for some rows")

    # if data.filter(pl.col("data_type").ne("GWAS")).select("cell_type").is_null().any():
    #     raise ValueError(
    #         f"Could not determine cell type for some non-GWAS rows in {filename}"
    #     )

    data.select(
        [
            "#dataset",
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
        ]
    ).sink_csv(output_filename, separator="\t", null_value="NA")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python munge_coloc_credset_file.py <filename> <eqtl_cat_metadata_filename> <output_filename>"
        )
        sys.exit(1)
    filename = sys.argv[1]
    eqtl_cat_metadata_filename = sys.argv[2]
    output_filename = sys.argv[3]
    main(filename, eqtl_cat_metadata_filename, output_filename)
