import sys

import polars as pl

# usage: munge_hpa.py <normal_ihc_data.tsv> <HPA version> <gencode genes tsv>
# the gencode release must be the one HPA's Ensembl version corresponds to (109 -> v43)
ihc_file, ihc_version, annotation_file = sys.argv[1:4]

# HPA scores each gene's IHC annotation once, best to worst: Enhanced, Supported, Approved,
# Uncertain. Only Uncertain is left out
RELIABLE = ["Enhanced", "Supported", "Approved"]

ihc = (
    pl.read_csv(ihc_file, separator="\t")
    .filter(pl.col("Reliability").is_in(RELIABLE))
    # one source row has no level at all
    .filter(pl.col("Level") != "N/A")
    .with_columns(
        pl.concat_str(
            pl.col("Tissue"),
            pl.col("IHC tissue name"),
            pl.col("Cell type"),
            separator="|",
        ).alias("tissue_cell")
    )
    .drop(["Tissue", "IHC tissue name", "Cell type", "Reliability"])
    .with_columns(pl.lit(f"HPA_{ihc_version}").alias("#dataset"))
)
annotation = pl.read_csv(annotation_file, separator="\t").with_columns(
    pl.col("gene_id").str.split(".").list.first().alias("ensg")
)

data = annotation.join(ihc, left_on="ensg", right_on="Gene", how="right").drop(
    ["Gene name", "gene_strand", "gene_type", "Gene"]
)

data = data.rename({c: c.replace(" ", "_").lower() for c in data.columns})

# replace spaces with underscores in all string columns
data = data.with_columns(pl.col(pl.String).str.replace_all(" ", "_"))

data = data.with_columns(pl.col("tissue_cell").str.to_lowercase())

# some ENSG ids are not present in gencode annotation, remove them
data.filter(pl.col("chrom").is_not_null()).select(
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
).sort(["chrom", "gene_start", "gene_end", "gene_id"]).write_csv(
    f"hpa_normal_ihc_data_{ihc_version}.long.tsv", separator="\t", null_value="NA"
)
