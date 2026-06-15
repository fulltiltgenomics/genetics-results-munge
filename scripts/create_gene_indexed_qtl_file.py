import polars as pl


def merge(
    data_file: str,
    gene_mapping: pl.LazyFrame,
    map_from_column: str,
    map_to_column: str,
) -> pl.LazyFrame:
    # lazy plan; sinked with the streaming engine (low memory) below. order does not
    # matter because the output is sorted by gene position afterwards.
    return (
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
        )
        .drop(["ensg", "gene_name"], strict=False)
    )


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

data_file = "/mnt/disks/data/eqtl_cat_r8_run/eQTL_Catalogue_R8.tsv"
gene_mapping_file = "/mnt/disks/data/gencode.v39.annotation.genes.tsv"
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

merged_lf = merge(
    data_file,
    gene_mapping,
    map_from_column,
    map_to_column,
)

out_base = data_file.removesuffix(".tsv.gz").removesuffix(".gz").removesuffix(".tsv")

# stream (low memory) the gene-mapped rows and, separately, the unmapped ones
merged_lf.filter(pl.col("trait_chr").is_not_null()).sink_csv(
    f"{out_base}.qtl.tsv", separator="\t", null_value="NA"
)
merged_lf.filter(pl.col("trait_chr").is_null()).sink_csv(
    f"{out_base}.qtl.unmapped.tsv", separator="\t", null_value="NA"
)

n_unmapped = (
    merged_lf.filter(pl.col("trait_chr").is_null())
    .select(pl.col("trait").n_unique())
    .collect()
    .item()
)
print(f"{n_unmapped} unmapped genes")

# eQTL Catalogue qtl schema is 19 data cols + trait_chr/trait_start/trait_end (cols 20-22):
# sort -T . -k20,20g -k21,21g -k22,22g -k6,6g -k7,7g -k8,8 -k9,9 -k3,3 <out>.qtl.tsv \
#   | bgzip -@4 > <out>.qtl.tsv.gz && tabix -f -s20 -b21 -e21 <out>.qtl.tsv.gz
