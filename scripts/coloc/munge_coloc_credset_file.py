import polars as pl
import sys

# R14 coloc credible set munging.
# Input columns: trait region rsid cs low_purity p beta se mlogp cs_specific_prob dataset tissue quant
# dataset/tissue/quant are now explicit columns (in R13 they were packed into one string).

# quant method -> data_type (applies to all rows; in R13 this was derived from dataset suffixes)
QUANT_DATA_TYPE = {
    "genome": "GWAS",
    "metabolite": "metaboQTL",
    "protein": "pQTL",
    "aptamer": "pQTL",
    "ge": "eQTL",
    "exon": "eQTL",
    "tx": "eQTL",
    "txrev": "eQTL",
    "microarray": "eQTL",
    "leafcutter": "sQTL",
    "majiq": "sQTL",
}

# FinnGen-internal / external datasets: bare dataset name -> (#dataset display name, cell_type).
# Everything not listed here is treated as eQTL Catalogue and resolved via the catalogue metadata.
# cell_type is None for GWAS datasets (no cell type), constant otherwise.
FINNGEN_DATASETS = {
    "FIN-R12-Metabolism": ("FinnGen_NMR", "plasma"),
    "FIN-R12-Somascan": ("FinnGen_SomaScan", "plasma"),
    "FIN-R12-Olink-3K": ("FinnGen_Olink_3K", "plasma"),
    "FIN-R14-Olink-5K-b1-b2": ("FinnGen_Olink_5K", "plasma"),
    "FinnGen-Olink-pQTL-meta": ("FinnGen_Olink", "plasma"),
    "FinnGen-R14": ("FinnGen_R14", None),
    "GeneRisk": ("GeneRisk", None),
    "Kanta_labs_v3": ("FinnGen_kanta", None),
    "UKB-finucane": ("UKB_Finucane", None),
    "UKB-PPP": ("UKB_PPP", "plasma"),
    "INTERVAL": ("INTERVAL", "plasma"),
    "FIN-SingleCell-b1": ("FinnGen_snRNAseq", "PBMC"),
    "FinnLiver": ("FinnLiver", "liver"),
}
# FinnLiver is FinnGen-internal but a multi-quant molecular QTL, so its trait_original
# also gets the quant suffix (like eQTL Catalogue traits)
QUANT_TAG_FINNGEN = {"FinnLiver"}

finngen_name_map = {k: v[0] for k, v in FINNGEN_DATASETS.items()}
finngen_cell_type_map = {k: v[1] for k, v in FINNGEN_DATASETS.items()}


def main(filename, eqtl_cat_metadata_filename, output_filename):
    # catalogue metadata: (study_label, sample_group, quant_method) -> (dataset_id, cell_type)
    # cell_type = tissue_label|condition_label with spaces replaced by underscores
    meta = (
        pl.read_csv(eqtl_cat_metadata_filename, separator="\t")
        .with_columns(
            pl.concat_str(pl.col("tissue_label"), pl.col("condition_label"), separator="|")
            .str.replace_all(" ", "_")
            .alias("cat_cell_type")
        )
        .select(
            pl.col("study_label").alias("dataset"),
            pl.col("sample_group").alias("tissue"),
            pl.col("quant_method").alias("quant"),
            pl.col("dataset_id").alias("cat_dataset"),
            "cat_cell_type",
        )
        .unique(subset=["dataset", "tissue", "quant"])
        .lazy()
    )

    data = (
        pl.scan_csv(
            filename,
            separator="\t",
            infer_schema_length=100000,
            schema_overrides={
                "trait": pl.Utf8,
                "p": pl.Float64,
                "beta": pl.Float64,
                "se": pl.Float64,
                "mlogp": pl.Float64,
                "low_purity": pl.Utf8,  # 0 / 1 / NA (NA for catalogue datasets)
                "dataset": pl.Utf8,
                "tissue": pl.Utf8,
                "quant": pl.Utf8,
            },
            null_values=["NA"],
        )
        # keep good and unknown-purity credible sets, drop only explicitly low purity (=1)
        .filter(pl.col("low_purity").ne("1").fill_null(True))
        .filter(pl.col("beta").is_not_null())
        .with_columns(
            pl.col("mlogp").round(4).alias("mlog10p"),
            pl.concat_str(pl.col("region"), pl.col("cs"), separator="_").alias("cs_id"),
            pl.col("rsid").str.split("_").list.get(0).str.replace("chr", "").str.replace("X", "23").cast(pl.Int8).alias("chr"),
            pl.col("rsid").str.split("_").list.get(1).cast(pl.Int32).alias("pos"),
            pl.col("rsid").str.split("_").list.get(2).alias("ref"),
            pl.col("rsid").str.split("_").list.get(3).alias("alt"),
            pl.col("quant").replace_strict(QUANT_DATA_TYPE, default=None).alias("data_type"),
            pl.col("dataset").is_in(list(FINNGEN_DATASETS.keys())).alias("_is_finngen"),
        )
        .join(meta, on=["dataset", "tissue", "quant"], how="left")
        .with_columns(
            pl.col("beta").map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8).alias("beta"),
            pl.col("se").map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8).alias("se"),
            pl.col("cs_specific_prob").round(4).alias("pip"),
            pl.when(pl.col("_is_finngen"))
            .then(pl.col("dataset").replace_strict(finngen_name_map, default=None))
            .otherwise(pl.col("cat_dataset"))
            .alias("#dataset"),
            pl.when(pl.col("_is_finngen"))
            .then(pl.col("dataset").replace_strict(finngen_cell_type_map, default=None))
            .otherwise(pl.col("cat_cell_type"))
            .alias("cell_type"),
            # quant suffix on trait_original for molecular QTLs (catalogue + FinnLiver)
            pl.when(~pl.col("_is_finngen") | pl.col("dataset").is_in(list(QUANT_TAG_FINNGEN)))
            .then(pl.concat_str(pl.col("trait"), pl.col("quant"), separator="|"))
            .otherwise(pl.col("trait"))
            .alias("trait_original"),
        )
    )

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
        print("Usage: python munge_coloc_credset_file.py <filename> <eqtl_cat_metadata_filename> <output_filename>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
