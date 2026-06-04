import polars as pl

median_tpm_file = (
    "/mnt/disks/dist_data/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"
)
annotation_file = "/mnt/disks/dist_data/gencode.v39.annotation.genes.tsv"

median_tpm = pl.read_csv(median_tpm_file, separator="\t", skip_lines=2)
annotation = pl.read_csv(annotation_file, separator="\t")

df = (
    annotation.join(median_tpm, left_on="gene_id", right_on="Name", how="right")
    .drop(["Description", "gene_strand", "gene_type"])
    .rename({"Name": "gene_id"})
    .with_columns(pl.lit("GTEx_v10").alias("#dataset"))
)

# write wide format (original)
df.write_csv("gtex_v10_median_tpm.tsv", separator="\t", null_value="NA")

# convert to longitudinal format: melt tissue columns into `tissue_cell`
id_cols = ["#dataset", "chrom", "gene_start", "gene_end", "gene_name", "gene_id"]
long_df = df.unpivot(index=id_cols, variable_name="tissue_cell", value_name="level")

long_df = long_df.with_columns(pl.col("tissue_cell").str.to_lowercase())

# write long format
long_df.select(
    [
        "#dataset",
        "chrom",
        "gene_start",
        "gene_end",
        "gene_name",
        "gene_id",
        "tissue_cell",
        "level",
    ]
).filter(pl.col("chrom").is_not_null()).sort(
    ["chrom", "gene_start", "gene_end", "gene_id"]
).write_csv(
    "gtex_v10_median_tpm.long.tsv", separator="\t", null_value="NA"
)

# compute z-scores for each gene
# (TPM_tissue - mean_TPM_all_tissues) / SD_TPM_all_tissues
# then z > e.g. 2 (as long as TPM > 1) is "highly expressed"

# TODO make tissue names lowercase
