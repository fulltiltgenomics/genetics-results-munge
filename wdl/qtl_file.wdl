version 1.0

task qtl_file {

    input {
        String docker
        File data_file
        File gene_mapping_file
        String map_from_column
        String map_to_column
        String output_file
        String output_unmapped_file
    }

    meta {
        description: "Create a QTL file (gene chromosome, start and end added to data, indexed by gene start)"
    }

    parameter_meta {
        docker: {
            help: "Docker image to use"
        }
        data_file: {
            help: "Location of data file"
        }
        gene_mapping_file: {
            help: "Location of gene mapping file (tsv with at least these columns: gene_id, chrom, gene_start, gene_end, gene_name)"
        }
        map_from_column: {
            help: "Column name to use for mapping in data file"
        }
        map_to_column: {
            help: "Column name to use for mapping in gene mapping file"
        }
        output_file: {
            help: "Output filename"
        }
        output_unmapped_file: {
            help: "Output filename for unmapped genes (tsv)"
        }
    }

    output {
        File qtl_file = output_file
        File qtl_file_tbi = output_file + ".tbi"
        File unmapped_file = output_unmapped_file
    }
    
    runtime {
        docker: docker
        memory: "16 GB"
        cpu: 2
        disks: "local-disk 200 HDD"
        preemptible: 2
        noAddress: true
    }

    command <<<

        set -euxo pipefail

        echo $(date) "create qtl file"

        python <<EOF > "~{output_file}.log"

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
                    maintain_order="right", # TODO need to figure sorting out here, now sorting afterwards; streaming=True for collect() is deprecated but not updated yet in polars 1.30.0
                )
                .drop(["ensg", "gene_name"], strict=False)
                .collect()
            )
            unmapped = df.filter(pl.col("trait_chr").is_null())
            return df, unmapped

        data_file = "~{data_file}"
        gene_mapping_file = "~{gene_mapping_file}"
        map_from_column = "~{map_from_column}"
        map_to_column = "~{map_to_column}"

        gene_mapping = (
            pl.scan_csv(
                gene_mapping_file,
                separator="\t",
                schema_overrides={
                    "gene_id": pl.Utf8,
                    "chrom": pl.Utf8,
                    "gene_start": pl.Int32,
                    "gene_end": pl.Int32,
                    "gene_name": pl.Utf8,
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
                "~{output_unmapped_file}",
                separator="\t",
                null_value="NA",
            )
        print(f"{len(set(unmapped.select('trait').to_series().to_list()))} unmapped genes")

        merged_df.filter(~pl.col("trait_chr").is_null()).write_csv(
            "~{output_file}.temp",
            separator="\t",
            null_value="NA",
        )
        EOF

        echo $(date) "sort, bgzip and tabix qtl file"
        # sort key is gene chr, gene start pos, gene end pos, data chr, data pos, data ref, data alt, trait
        sort -T . -k20,20g -k21,21g -k22,22g -k6,6g -k7,7g -k8,8 -k9,9 -k3,3 ~{output_file}.temp \
        | bgzip -@4 > ~{output_file} && tabix -@4 -s20 -b21 -e21 ~{output_file}

        echo $(date) "done"

    >>>

}
