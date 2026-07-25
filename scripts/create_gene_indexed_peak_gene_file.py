"""Build the gene-indexed copy of the Open4Gene peak-to-gene table that backs the API's
`/gene_to_peaks/{gene}` endpoint.

The published table is tabix-indexed on PEAK coordinates, which answers "which genes is this
peak linked to" but has no inverse: a gene's peaks are scattered across the file. This script
appends the linked gene's coordinates as three trailing columns and sorts on them, so the same
rows can be queried by gene locus.

The three appended columns exist only for the index - the API drops them and re-derives gene
coordinates from the GENCODE version the request asks for, so `/gene_to_peaks` and
`/peak_to_genes` return the identical column set. `gene_chrom` therefore uses the API's
chromosome convention (numeric, X=23) rather than the `chr1` form of the peak columns.
"""

import argparse

import polars as pl

from create_caqtl_gene_indexed_qtl_file import gene_coordinates

OUTPUT_COLUMNS = ["gene_chrom", "gene_start", "gene_end"]


def gene_indexed_peaks(open4gene_file: str, gene_coords: pl.LazyFrame) -> pl.LazyFrame:
    return (
        pl.scan_csv(
            open4gene_file,
            separator="\t",
            schema_overrides={"gene_id": pl.Utf8},
            null_values=["NA"],
        )
        .join(gene_coords, left_on="gene_id", right_on="ensg", how="inner")
        .rename({"trait_chr": "gene_chrom", "trait_start": "gene_start", "trait_end": "gene_end"})
        .drop("gene_name")
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--open4gene", required=True, help="Open4Gene peak-to-gene results TSV(.gz)"
    )
    parser.add_argument(
        "--gene-mapping",
        required=True,
        help="GENCODE genes TSV of the version the Open4Gene gene ids come from (v32)",
    )
    parser.add_argument(
        "--output", required=True, help="output TSV (unsorted, sort+bgzip+tabix after)"
    )
    args = parser.parse_args()

    gene_coords = gene_coordinates(args.gene_mapping)
    out = gene_indexed_peaks(args.open4gene, gene_coords)

    columns = out.collect_schema().names()
    if columns[-3:] != OUTPUT_COLUMNS:
        raise SystemExit(f"unexpected output columns: {columns}")

    dropped = (
        pl.scan_csv(
            args.open4gene,
            separator="\t",
            schema_overrides={"gene_id": pl.Utf8},
            null_values=["NA"],
        )
        .join(gene_coords, left_on="gene_id", right_on="ensg", how="anti")
        .select(pl.len())
        .collect()
        .item()
    )
    print(f"{len(columns)} output columns, {dropped} rows dropped without gene coordinates")

    out.sink_csv(args.output, separator="\t", null_value="NA")


if __name__ == "__main__":
    main()
