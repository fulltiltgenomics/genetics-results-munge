import polars as pl


def merge(
    data_file: str,
    gene_mapping: pl.DataFrame,
    map_from_column: str,
    map_to_column: str,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    df = (
        pl.scan_csv(
            data_file,
            separator="\t",
            schema_overrides={
                "#dataset": pl.Utf8,
                "data_type": pl.Utf8,
                "trait": pl.Utf8,
                "trait_original": pl.Utf8,
                "cell_type": pl.Utf8,
                "chr": pl.Int8,
                "pos": pl.Int32,
                "ref": pl.Utf8,
                "alt": pl.Utf8,
                "mlog10p": pl.Float64,
                "beta": pl.Float64,
                "se": pl.Float64,
                "pip": pl.Float64,
                "cs_id": pl.Utf8,
                "cs_size": pl.Int32,
                "cs_min_r2": pl.Float64,
                "aaf": pl.Float64,
                "most_severe": pl.Utf8,
                "gene_most_severe": pl.Utf8,
            },
            infer_schema_length=100000,  # make sure the columns are correct
            null_values=["NA"],
        )
        .join(
            gene_mapping,
            left_on=map_from_column,
            right_on=map_to_column,
            how="left",
            maintain_order="right",
        )
        .drop(["ensg", "gene_name"], strict=False)
        # .sort(
        #     "trait_chr", "trait_start", "trait_end", "chr", "pos", "ref", "alt", "trait"
        # )
        .collect()
    )
    unmapped = df.filter(pl.col("trait_chr").is_null())
    return df, unmapped


# data_file = "/mnt/disks/dist_data/FinnGen_Olink_1-4_credible_sets.tsv.gz"
# gene_mapping_file = "/mnt/disks/dist_data/gencode.v49.annotation.genes.tsv"
# map_from_column = "trait"
# map_to_column = "gene_name"

# data_file = "/mnt/disks/dist_data/FinnGen_snRNAseq_202509_credible_sets.tsv.gz"
# gene_mapping_file = "/mnt/disks/dist_data/gencode.v32.annotation.genes.tsv"
# map_from_column = "trait_original"
# map_to_column = "ensg"

# data_file = "/mnt/disks/dist_data/UKB_PPP_credible_sets.tsv.gz"
# gene_mapping_file = "/mnt/disks/dist_data/gencode.v49.annotation.genes.tsv"
# map_from_column = "trait"
# map_to_column = "gene_name"

data_file = "/mnt/disks/dist_data/eQTL_Catalogue_R7.tsv.gz"
gene_mapping_file = "/mnt/disks/dist_data/gencode.v39.annotation.genes.tsv"
map_from_column = "trait"
map_to_column = "gene_name"

gene_mapping = (
    pl.scan_csv(
        gene_mapping_file,
        separator="\t",
        # gene_id chrom   gene_start      gene_end        gene_strand     gene_name       gene_type
        schema_overrides={
            "gene_id": pl.Utf8,
            "chrom": pl.Utf8,
            "gene_start": pl.Int32,
            "gene_end": pl.Int32,
            "gene_strand": pl.Utf8,
            "gene_name": pl.Utf8,
            "gene_type": pl.Utf8,
        },
    )
    .with_columns(
        pl.col("chrom")
        .str.replace(r"^X$", "23")
        .str.replace(r"^Y$", "24")
        .str.replace(r"^M$", "26")
        .cast(pl.Int8)
        .alias("trait_chr")
    )
    .filter(pl.col("trait_chr").lt(24))
    .rename(
        {
            "gene_start": "trait_start",
            "gene_end": "trait_end",
        }
    )
    .with_columns(pl.col("gene_id").str.split(".").list.first().alias("ensg"))
    .select("ensg", "gene_name", "trait_chr", "trait_start", "trait_end")
)

merged_df, unmapped = merge(
    data_file,
    gene_mapping,
    map_from_column,
    map_to_column,
)

if len(unmapped) > 0:
    unmapped.write_csv(
        f"{data_file.removesuffix('.tsv.gz').removesuffix('.gz').removesuffix('.tsv')}.qtl.unmapped.tsv",
        separator="\t",
        null_value="NA",
    )
print(f"{len(set(unmapped.select('trait').to_series().to_list()))} unmapped genes")


merged_df.filter(~pl.col("trait_chr").is_null()).write_csv(
    f"{data_file.removesuffix('.tsv.gz').removesuffix('.gz').removesuffix('.tsv')}.qtl.tsv",
    separator="\t",
    null_value="NA",
)

# bgzip -f FinnGen_Olink.qtl.tsv && tabix -f -s19 -b20 -e20 FinnGen_Olink.qtl.tsv.gz
