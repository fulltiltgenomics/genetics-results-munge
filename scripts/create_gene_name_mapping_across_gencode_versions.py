### this script creates a gene name mapping file across all gencode versions
### output file:
### ensg    gene_name_49    gene_name_43    gene_name_39    gene_name_32
### ENSG00000290825 DDX11L16        DDX11L2 NA      NA
### ENSG00000223972 DDX11L1 DDX11L1 DDX11L1 DDX11L1
### ...
### ENSG00000234290 NA      NA      NA      AC116366.1

import polars as pl

gencode_versions = ["49", "43", "39", "35", "32", "19"]
data = None
for version in gencode_versions:
    df = (
        pl.read_csv(
            f"/mnt/disks/dist_data/gencode.v{version}.annotation.genes.tsv",
            separator="\t",
        )
        .with_columns(pl.col("gene_id").str.split(".").list.first().alias("ensg"))
        .select(["ensg", "gene_name"])
        .rename({"gene_name": f"gene_name_{version}"})
    )
    if data is None:
        data = df
    else:
        data = data.join(df, on="ensg", how="full", coalesce=True)

data.write_csv(
    f"/mnt/disks/dist_data/gencode_gene_name_mapping_{'-'.join(gencode_versions)}.tsv",
    separator="\t",
    null_value="NA",
)
