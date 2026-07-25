"""Build the gene-indexed caQTL credible set file (`*.qtl.tsv.gz`) served by the API's
`/credible_sets_by_qtl_gene/{gene}` endpoint.

Unlike eQTL/pQTL, the caQTL molecular trait is a chromatin peak, not a gene, so the trait
cannot be mapped to gene coordinates directly (as `create_gene_indexed_qtl_file.py` does).
The gene link comes from the Open4Gene peak-to-gene table: each credible set row is joined
to the genes its peak is linked to, and the row is emitted once per linked gene with that
gene's coordinates in `trait_chr`/`trait_start`/`trait_end`.

The join is on BOTH peak and cell type, so a signal is only attributed to a gene in the cell
type where that peak-gene link was actually significant.

Output columns are exactly those of the other `*.qtl.tsv.gz` files (19 credible set columns +
trait_chr/trait_start/trait_end): the API takes the merged response header from whichever qtl
file it reads first, so a differing column set would misalign the cross-resource response.
Within those columns `trait` becomes the linked gene symbol and `trait_original` keeps the
peak id (`cs_id` also stays peak-based, so the signal identity is preserved).

Gene coordinates MUST come from the same GENCODE file the API resolves query genes with
(v32 for finngen_caqtl - see `gencode_version` in the API's credible_sets profile config):
the API filters returned rows by an exact string match on trait_start/trait_end, so
coordinates from any other version silently match nothing.
"""

import argparse

import polars as pl

CS_SCHEMA = {
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
}

OPEN4GENE_SCHEMA = {
    "peak_id": pl.Utf8,
    "gene_id": pl.Utf8,
    "symbol": pl.Utf8,
    "cell_type": pl.Utf8,
}

OUTPUT_COLUMNS = list(CS_SCHEMA.keys()) + ["trait_chr", "trait_start", "trait_end"]

# Open4Gene cell types are prefixed ("predicted.celltype.l1.B"), credible set ones are not
# ("l1.B"); the two vocabularies are otherwise identical
CELL_TYPE_PREFIX = "predicted.celltype."


def gene_coordinates(gene_mapping_file: str) -> pl.LazyFrame:
    """Read GENCODE gene coordinates in the API's chromosome convention (X=23, no Y/MT)."""
    return (
        pl.scan_csv(
            gene_mapping_file,
            separator="\t",
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
        .rename({"gene_start": "trait_start", "gene_end": "trait_end"})
        .with_columns(pl.col("gene_id").str.split(".").list.first().alias("ensg"))
        .select("ensg", "gene_name", "trait_chr", "trait_start", "trait_end")
    )


def peak_to_gene_links(
    open4gene_file: str, gene_coords: pl.LazyFrame
) -> pl.LazyFrame:
    """One row per (peak, cell type, linked gene) with the gene's coordinates."""
    return (
        pl.scan_csv(
            open4gene_file,
            separator="\t",
            schema_overrides=OPEN4GENE_SCHEMA,
            null_values=["NA"],
        )
        .select(
            "peak_id",
            "gene_id",
            pl.col("cell_type").str.strip_prefix(CELL_TYPE_PREFIX).alias("cell_type"),
        )
        # the input is already one row per (peak, gene, cell type); deduplicating guards
        # against a future input that repeats links, which would multiply credible set rows
        .unique()
        .join(gene_coords, left_on="gene_id", right_on="ensg", how="inner")
        .select("peak_id", "cell_type", "gene_name", "trait_chr", "trait_start", "trait_end")
    )


def gene_indexed_credible_sets(
    credible_set_file: str, links: pl.LazyFrame
) -> pl.LazyFrame:
    return (
        pl.scan_csv(
            credible_set_file,
            separator="\t",
            schema_overrides=CS_SCHEMA,
            infer_schema_length=100000,
            null_values=["NA"],
        )
        .join(
            links,
            left_on=["trait", "cell_type"],
            right_on=["peak_id", "cell_type"],
            how="inner",
        )
        # both expressions see the pre-update frame, so the peak id lands in trait_original
        # even though trait is replaced in the same call
        .with_columns(
            pl.col("gene_name").alias("trait"),
            pl.col("trait").alias("trait_original"),
        )
        .select(OUTPUT_COLUMNS)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--credible-sets",
        required=True,
        help="caQTL credible set TSV (uncompressed; the 19-column munged format)",
    )
    parser.add_argument(
        "--open4gene", required=True, help="Open4Gene peak-to-gene results TSV(.gz)"
    )
    parser.add_argument(
        "--gene-mapping",
        required=True,
        help="GENCODE genes TSV of the version the API uses for this dataset (v32)",
    )
    parser.add_argument(
        "--output", required=True, help="output TSV (unsorted, sort+bgzip+tabix after)"
    )
    args = parser.parse_args()

    gene_coords = gene_coordinates(args.gene_mapping)
    links = peak_to_gene_links(args.open4gene, gene_coords)

    link_stats = (
        links.select(
            pl.len().alias("links"),
            pl.col("peak_id").n_unique().alias("peaks"),
            pl.col("gene_name").n_unique().alias("genes"),
            pl.col("cell_type").n_unique().alias("cell_types"),
        )
        .collect()
        .row(0)
    )
    print(
        f"peak-to-gene links: {link_stats[0]} over {link_stats[1]} peaks, "
        f"{link_stats[2]} genes, {link_stats[3]} cell types"
    )

    unmapped_genes = (
        pl.scan_csv(
            args.open4gene,
            separator="\t",
            schema_overrides=OPEN4GENE_SCHEMA,
            null_values=["NA"],
        )
        .select("gene_id")
        .unique()
        .join(gene_coords, left_on="gene_id", right_on="ensg", how="anti")
        .select(pl.len())
        .collect()
        .item()
    )
    print(f"{unmapped_genes} Open4Gene genes without coordinates in the GENCODE file")

    gene_indexed_credible_sets(args.credible_sets, links).sink_csv(
        args.output, separator="\t", null_value="NA"
    )


if __name__ == "__main__":
    main()
