import polars as pl
import sys

# R14 colocQC (H4 / QC table) munging.
# dataset/tissue/quant are explicit per side (dataset1/tissue1/quant1, dataset2/tissue2/quant2);
# hit stats are explicit columns (hit1_beta/hit1_se/hit1_p/hit1_mlogp, ...) instead of the
# old comma-packed hit1_info/hit2_info. Trait gene-name mapping happens in the second pass
# (munge_coloc_file_map_traits.py).

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
finngen_name_map = {k: v[0] for k, v in FINNGEN_DATASETS.items()}
finngen_cell_type_map = {k: v[1] for k, v in FINNGEN_DATASETS.items()}
finngen_keys = list(FINNGEN_DATASETS.keys())

LOW_PURITY_TRUE = ["1", "TRUE", "True", "true"]


def _sci(x):
    return f"{x:.3e}" if x is not None else None


def _dataset_cols(side, meta):
    """expressions to derive #dataset/cell_type/data_type for one side (1 or 2)."""
    ds = pl.col(f"dataset{side}")
    is_fg = ds.is_in(finngen_keys)
    return [
        pl.col(f"quant{side}").replace_strict(QUANT_DATA_TYPE, default=None).alias(f"data_type{side}"),
        pl.when(is_fg)
        .then(ds.replace_strict(finngen_name_map, default=None))
        .otherwise(pl.col(f"cat_dataset{side}"))
        .alias(f"_disp_dataset{side}"),
        pl.when(is_fg)
        .then(ds.replace_strict(finngen_cell_type_map, default=None))
        .otherwise(pl.col(f"cat_cell_type{side}"))
        .alias(f"cell_type{side}"),
    ]


def main(filename, eqtl_cat_metadata_filename, output_filename):
    meta_base = (
        pl.read_csv(eqtl_cat_metadata_filename, separator="\t")
        .with_columns(
            pl.concat_str(pl.col("tissue_label"), pl.col("condition_label"), separator="|")
            .str.replace_all(" ", "_")
            .alias("cat_cell_type")
        )
        .select("study_label", "sample_group", "quant_method", "dataset_id", "cat_cell_type")
        .unique(subset=["study_label", "sample_group", "quant_method"])
    )

    def meta_for(side):
        return meta_base.rename(
            {
                "study_label": f"dataset{side}",
                "sample_group": f"tissue{side}",
                "quant_method": f"quant{side}",
                "dataset_id": f"cat_dataset{side}",
                "cat_cell_type": f"cat_cell_type{side}",
            }
        ).lazy()

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
                "dataset1": pl.Utf8,
                "dataset2": pl.Utf8,
                "tissue1": pl.Utf8,
                "tissue2": pl.Utf8,
                "quant1": pl.Utf8,
                "quant2": pl.Utf8,
                "low_purity1": pl.Utf8,
                "low_purity2": pl.Utf8,
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
                "hit1_beta": pl.Float64,
                "hit2_beta": pl.Float64,
                "hit1_mlogp": pl.Float64,
                "hit2_mlogp": pl.Float64,
            },
            null_values=["NA"],
        )
        # standard FinnGen filtering; drop only explicitly low purity (1/TRUE), keep 0/FALSE/NA
        # (low_purity is null for eQTL Catalogue datasets; is_in on null -> null, so fill_null)
        .filter(
            ~pl.col("low_purity1").is_in(LOW_PURITY_TRUE).fill_null(False)
            & ~pl.col("low_purity2").is_in(LOW_PURITY_TRUE).fill_null(False)
        )
        .filter((pl.col("probmass_1") > 0.9) & (pl.col("probmass_2") > 0.9))
        .filter((pl.col("cs1_log10bf") > 0.9) & (pl.col("cs2_log10bf") > 0.9))
        .with_columns(
            pl.concat_str(pl.col("region1"), pl.col("cs1"), separator="_").alias("cs1_id"),
            pl.concat_str(pl.col("region2"), pl.col("cs2"), separator="_").alias("cs2_id"),
            pl.col("region1").str.split(":").list.get(0).str.replace("chr", "").str.replace("X", "23").cast(pl.Int8).alias("chr"),
            pl.col("region1").str.replace(r"X", "23").str.extract(r"^(chr\d+):(-?\d+)-(\d+)$", 2).cast(pl.Int32).alias("region1_start"),
            pl.col("region1").str.replace(r"X", "23").str.extract(r"^(chr\d+):(-?\d+)-(\d+)$", 3).cast(pl.Int32).alias("region1_end"),
            pl.col("region2").str.replace(r"X", "23").str.extract(r"^(chr\d+):(-?\d+)-(\d+)$", 2).cast(pl.Int32).alias("region2_start"),
            pl.col("region2").str.replace(r"X", "23").str.extract(r"^(chr\d+):(-?\d+)-(\d+)$", 3).cast(pl.Int32).alias("region2_end"),
            pl.col("hit1").str.replace_all("_", ":").str.replace("chr", "").str.replace("X", "23").alias("hit1"),
            pl.col("hit2").str.replace_all("_", ":").str.replace("chr", "").str.replace("X", "23").alias("hit2"),
        )
        .with_columns(
            pl.max_horizontal([pl.col("region1_start"), pl.lit(1)]).alias("region1_start"),
            pl.max_horizontal([pl.col("region2_start"), pl.lit(1)]).alias("region2_start"),
        )
        .with_columns(
            pl.min_horizontal([pl.col("region1_start"), pl.col("region2_start")]).alias("region_start_min"),
            pl.max_horizontal([pl.col("region1_end"), pl.col("region2_end")]).alias("region_end_max"),
        )
        .with_columns(
            pl.col("hit1_mlogp").round(4).alias("hit1_mlog10p"),
            pl.col("hit2_mlogp").round(4).alias("hit2_mlog10p"),
            pl.col("hit1_beta").alias("hit1_beta"),
            pl.col("hit2_beta").alias("hit2_beta"),
            pl.col("PP.H0.abf").map_elements(_sci, return_dtype=pl.Utf8),
            pl.col("PP.H1.abf").map_elements(_sci, return_dtype=pl.Utf8),
            pl.col("PP.H2.abf").map_elements(_sci, return_dtype=pl.Utf8),
            pl.col("PP.H3.abf").map_elements(_sci, return_dtype=pl.Utf8),
            pl.col("PP.H4.abf").map_elements(_sci, return_dtype=pl.Utf8),
            pl.when(pl.col("cs1_log10bf").is_infinite()).then(1000).otherwise(pl.col("cs1_log10bf").round(4)).alias("cs1_log10bf"),
            pl.when(pl.col("cs2_log10bf").is_infinite()).then(1000).otherwise(pl.col("cs2_log10bf").round(4)).alias("cs2_log10bf"),
            pl.col("clpp").map_elements(_sci, return_dtype=pl.Utf8),
            pl.col("clpa").map_elements(_sci, return_dtype=pl.Utf8),
        )
        .join(meta_for(1), on=["dataset1", "tissue1", "quant1"], how="left")
        .join(meta_for(2), on=["dataset2", "tissue2", "quant2"], how="left")
        .with_columns(_dataset_cols(1, meta_base) + _dataset_cols(2, meta_base))
        .rename({"_disp_dataset1": "#dataset1", "_disp_dataset2": "dataset2_disp"})
    )

    data.select(
        [
            "#dataset1",
            pl.col("dataset2_disp").alias("dataset2"),
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
        print("Usage: python munge_coloc_file.py <filename> <eqtl_cat_metadata_filename> <output_filename>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
