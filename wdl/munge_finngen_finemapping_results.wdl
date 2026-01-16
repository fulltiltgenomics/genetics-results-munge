version 1.0

import "qtl_file.wdl" as qtl_file

workflow munge_finngen_finemapping_results {

    input {
        String docker
        String dataset
        File metadata_file
        Boolean trait_is_cell_type
        File? trait_mapping_file
        Boolean? require_trait_mapping
        File? variant_annotation_file
        String? cell_type
        String output_file
        Boolean create_qtl_file
    }

    meta {
        description: "Join and merge finngen format finemapping results (.snp.filter.tsv and .cred.summary.tsv files) into a single tabixed file"
    }

    parameter_meta {
        docker: {
            help: "Docker image to use"
        }
        dataset: {
            help: "Dataset name (e.g. FinnGen_R13)"
        }
        metadata_file: {
            help: "Tab-separated headerless file one trait per line with the following columns: data_type (e.g. GWAS or eQTL), trait (phenocode), snp_file, cred_file"
        }
        trait_is_cell_type: {
            help: "Whether the trait is not a trait but a cell type (e.g. predicted.celltype.CD4_T_cell)"
        }
        trait_mapping_file: {
            help: "Optional location of headerless trait mapping file (tsv with at least two columns: trait and mapped_trait)"
        }
        require_trait_mapping: {
            help: "Optional, whether to require trait mapping (if not given or false, traits that are not in the trait mapping file will be set to the original trait name)"
        }
        variant_annotation_file: {
            help: "Optional location of variant annotation file (tsv file with at least these columns: #variant, most_severe, gene_most_severe)"
        }
        cell_type: {
            help: "Cell type (optional, e.g. plasma)"
        }
        output_file: {
            help: "Output filename"
        }
        create_qtl_file: {
            help: "Whether to create a QTL file (gene chromosome, start and end added to data, indexed by gene start)"
        }
    }

    Array[Array[String]] metadata = transpose(read_tsv(metadata_file))

    call join_and_merge {
        input:
        docker = docker,
        dataset = dataset,
        trait_is_cell_type = trait_is_cell_type,
        cell_type = cell_type,
        trait_mapping_file = trait_mapping_file,
        require_trait_mapping = require_trait_mapping,
        variant_annotation_file = variant_annotation_file,
        data_types = metadata[0],
        traits = metadata[1],
        snp_files = metadata[2],
        cred_files = metadata[3],
        output_file = output_file
    }

    if (create_qtl_file) {
        call qtl_file.qtl_file as qtl_file {
            input:
            docker = docker,
            data_file = join_and_merge.munged_file,
            output_file = sub(output_file, ".tsv.gz", ".qtl.tsv.gz"),
            output_unmapped_file = sub(output_file, ".tsv.gz", ".qtl.unmapped.tsv")
        }
    }

    output {
        File munged_file = join_and_merge.munged_file
        File munged_file_tbi = join_and_merge.munged_file_tbi
        File log = join_and_merge.log
        Array[File] per_trait_munged_files = join_and_merge.per_trait_munged_files
        Array[File] per_trait_stats_files = join_and_merge.per_trait_stats_files
        File aggregate_stats_file = join_and_merge.aggregate_stats_file
        File? qtl_file_out = qtl_file.qtl_file
        File? qtl_file_out_tbi = qtl_file.qtl_file_tbi
        File? qtl_file_unmapped_out = qtl_file.unmapped_file
    }
}

task join_and_merge {

    input {
        String docker
        String dataset
        Array[String] data_types
        Array[String] traits
        Array[File] snp_files
        Array[File] cred_files
        Boolean trait_is_cell_type
        File? trait_mapping_file
        Boolean? require_trait_mapping
        File? variant_annotation_file
        String? cell_type
        String output_file
    }

    meta {
        description: "Join and merge finngen format finemapping results (.snp.filter.tsv and .cred.summary.tsv files) into a single tabixed file"
    }

    parameter_meta {
        docker: {
            help: "Docker image to use"
        }
        dataset: {
            help: "Dataset name (e.g. FinnGen_R13)"
        }
        data_types: {
            help: "Data types (e.g. GWAS or eQTL)"
        }
        traits: {
            help: "Traits (e.g. trait1, trait2, trait3)"
        }
        snp_files: {
            help: "Locations of snp.filter.tsv files"
        }
        cred_files: {
            help: "Locations of cred.summary.tsv files"
        }
        trait_is_cell_type: {
            help: "Whether the trait is not a trait but a cell type (e.g. predicted.celltype.CD4_T_cell)"
        }
        cell_type: {
            help: "Cell type (optional, e.g. plasma)"
        }
        require_trait_mapping: {
            help: "Optional, whether to require trait mapping (if not given or false, traits that are not in the trait mapping file will be set to the original trait name)"
        }
        trait_mapping_file: {
            help: "Optional location of headerless trait mapping file (tsv with at least two columns: trait and mapped_trait)"
        }
        variant_annotation_file: {
            help: "Optional location of variant annotation file (tsv file with at least these columns: #variant, most_severe, gene_most_severe)"
        }
        output_file: {
            help: "Output filename"
        }
    }

    output {
        File munged_file = output_file
        File munged_file_tbi = output_file + ".tbi"
        File log = "~{dataset}_join_and_merge.log"
        Array[File] per_trait_munged_files = glob("individual/*.SUSIE.munged.tsv")
        # Array[File] per_trait_munged_files = read_lines("list_of_munged_files")
        Array[File] per_trait_stats_files = glob("individual/*.SUSIE.munged.stats.json")
        File aggregate_stats_file = "credible_set_stats.tsv"
    }

    runtime {
        docker: docker
        memory: "16 GB"
        cpu: 2
        disks: "local-disk 100 HDD"
        preemptible: 2
        noAddress: true
    }

    command <<<

        set -euxo pipefail

        echo $(date) "join"
        # join .snp.filter.tsv and .cred.summary.tsv files
        # write selected columns to a new file per trait/study that is sorted by chr, pos, ref, alt, trait
        mkdir -p "individual"

        cat <<'CREDIBLE_SET_STATS_EOF' > credible_set_stats.py
        import json
        import polars as pl
        from dataclasses import dataclass, asdict

        CODING_VARIANTS = {
            "missense_variant",
            "synonymous_variant",
            "stop_gained",
            "stop_lost",
            "start_lost",
            "frameshift_variant",
            "inframe_insertion",
            "inframe_deletion",
            "splice_acceptor_variant",
            "splice_donor_variant",
        }

        LOF_VARIANTS = {
            "stop_gained",
            "stop_lost",
            "start_lost",
            "frameshift_variant",
            "splice_acceptor_variant",
            "splice_donor_variant",
        }

        @dataclass
        class CredibleSetStats:
            trait: str
            trait_original: str
            dataset: str
            data_type: str
            n_risk_cs: int
            n_risk_cs_with_coding: int
            n_risk_cs_with_coding_pip_gt_0_05: int
            n_risk_cs_with_lof: int
            n_risk_cs_with_lof_pip_gt_0_05: int
            n_protective_cs: int
            n_protective_cs_with_coding: int
            n_protective_cs_with_coding_pip_gt_0_05: int
            n_protective_cs_with_lof: int
            n_protective_cs_with_lof_pip_gt_0_05: int

        def calculate_stats(df: pl.DataFrame) -> CredibleSetStats:
            trait = df["trait"][0] if len(df) > 0 else ""
            trait_original = (
                df["trait_original"][0]
                if len(df) > 0 and "trait_original" in df.columns
                else trait
            )
            dataset = df["dataset"][0] if len(df) > 0 else ""
            data_type = df["data_type"][0] if len(df) > 0 else ""

            if len(df) == 0:
                return CredibleSetStats(
                    trait=trait, trait_original=trait_original, dataset=dataset, data_type=data_type,
                    n_risk_cs=0, n_risk_cs_with_coding=0, n_risk_cs_with_coding_pip_gt_0_05=0,
                    n_risk_cs_with_lof=0, n_risk_cs_with_lof_pip_gt_0_05=0,
                    n_protective_cs=0, n_protective_cs_with_coding=0, n_protective_cs_with_coding_pip_gt_0_05=0,
                    n_protective_cs_with_lof=0, n_protective_cs_with_lof_pip_gt_0_05=0,
                )

            if df["pip"].dtype == pl.Utf8:
                df = df.with_columns(pl.col("pip").cast(pl.Float64, strict=False))
            if df["beta"].dtype == pl.Utf8:
                df = df.with_columns(pl.col("beta").cast(pl.Float64, strict=False))
            if "aaf" in df.columns and df["aaf"].dtype == pl.Utf8:
                df = df.with_columns(pl.col("aaf").cast(pl.Float64, strict=False))

            lead_variants = df.sort("pip", descending=True).group_by("cs_id").first()

            if "aaf" in lead_variants.columns:
                lead_variants = lead_variants.with_columns(
                    pl.when(pl.col("aaf").is_not_null())
                    .then(
                        pl.when(
                            ((pl.col("beta") > 0) & (pl.col("aaf") < 0.5))
                            | ((pl.col("beta") < 0) & (pl.col("aaf") > 0.5))
                        ).then(pl.lit("risk")).otherwise(pl.lit("protective"))
                    )
                    .otherwise( # TODO fix data so that AF is always available - this assumes alternate allele is minor allele and is bad
                        pl.when(pl.col("beta") > 0).then(pl.lit("risk")).otherwise(pl.lit("protective"))
                    )
                    .alias("direction")
                )
            else:
                lead_variants = lead_variants.with_columns(
                    pl.when(pl.col("beta") > 0).then(pl.lit("risk")).otherwise(pl.lit("protective")).alias("direction")
                )

            risk_cs_ids = set(lead_variants.filter(pl.col("direction") == "risk")["cs_id"].to_list())
            protective_cs_ids = set(lead_variants.filter(pl.col("direction") == "protective")["cs_id"].to_list())

            risk_variants = df.filter(pl.col("cs_id").is_in(risk_cs_ids))
            protective_variants = df.filter(pl.col("cs_id").is_in(protective_cs_ids))

            def count_cs_with_variants(variants_df, variant_set, pip_threshold=None):
                if "most_severe" not in variants_df.columns:
                    return 0
                filtered = variants_df.filter(pl.col("most_severe").is_in(variant_set))
                if pip_threshold is not None:
                    filtered = filtered.filter(pl.col("pip") > pip_threshold)
                return filtered["cs_id"].n_unique()

            return CredibleSetStats(
                trait=trait, trait_original=trait_original, dataset=dataset, data_type=data_type,
                n_risk_cs=len(risk_cs_ids),
                n_risk_cs_with_coding=count_cs_with_variants(risk_variants, CODING_VARIANTS),
                n_risk_cs_with_coding_pip_gt_0_05=count_cs_with_variants(risk_variants, CODING_VARIANTS, 0.05),
                n_risk_cs_with_lof=count_cs_with_variants(risk_variants, LOF_VARIANTS),
                n_risk_cs_with_lof_pip_gt_0_05=count_cs_with_variants(risk_variants, LOF_VARIANTS, 0.05),
                n_protective_cs=len(protective_cs_ids),
                n_protective_cs_with_coding=count_cs_with_variants(protective_variants, CODING_VARIANTS),
                n_protective_cs_with_coding_pip_gt_0_05=count_cs_with_variants(protective_variants, CODING_VARIANTS, 0.05),
                n_protective_cs_with_lof=count_cs_with_variants(protective_variants, LOF_VARIANTS),
                n_protective_cs_with_lof_pip_gt_0_05=count_cs_with_variants(protective_variants, LOF_VARIANTS, 0.05),
            )

        def write_stats_json(stats, output_path):
            with open(output_path, "w") as f:
                json.dump(asdict(stats), f, indent=2)

        def stats_to_tsv_row(stats):
            d = asdict(stats)
            return "\t".join(str(v) for v in d.values())

        def get_tsv_header():
            return "\t".join([
                "trait", "trait_original", "dataset", "data_type",
                "n_risk_cs", "n_risk_cs_with_coding", "n_risk_cs_with_coding_pip_gt_0_05",
                "n_risk_cs_with_lof", "n_risk_cs_with_lof_pip_gt_0_05",
                "n_protective_cs", "n_protective_cs_with_coding", "n_protective_cs_with_coding_pip_gt_0_05",
                "n_protective_cs_with_lof", "n_protective_cs_with_lof_pip_gt_0_05",
            ])
        CREDIBLE_SET_STATS_EOF

        python <<EOF > "~{dataset}_join_and_merge.log"

        import warnings
        from scipy.special import log_ndtr
        from typing import Literal
        import numpy as np
        import polars as pl
        from credible_set_stats import calculate_stats, write_stats_json, get_tsv_header, stats_to_tsv_row

        Datatype = Literal["GWAS", "eQTL", "sQTL", "pQTL", "edQTL", "metaboQTL"]

        class NoDataException(Exception):
            pass

        def merge(
            dataset: str,
            data_type: Datatype,
            snp_file: str,
            cred_file: str,
            cell_type: str | None = None,
            trait_is_cell_type: bool = False,
            trait_mapping: pl.DataFrame | None = None,
            require_trait_mapping: bool = False,
            variant_annotation: pl.DataFrame | None = None,
        ) -> tuple[pl.DataFrame, pl.DataFrame]:
            snp_df = pl.read_csv(
                snp_file,
                separator="\t",
                schema_overrides={
                    "trait": pl.Utf8,
                    "region": pl.Utf8,
                    "chromosome": pl.Utf8,
                    "position": pl.Int32,
                    "allele1": pl.Utf8,
                    "allele2": pl.Utf8,
                    "maf": pl.Float64,
                    "cs": pl.Int8,
                    "cs_specific_prob": pl.Float64,
                    "p": pl.Float64,
                    "beta": pl.Float64,
                    "se": pl.Float64,
                    "most_severe": pl.Utf8,
                    "gene_most_severe": pl.Utf8,
                },
                infer_schema_length=100000,  # needs to be large because some float fields may not otherwise be detected correctly
                null_values=["NA"],
            ).with_columns(
                pl.concat_str(
                    [pl.col("region"), pl.col("cs").cast(pl.Utf8)], separator="_"
                ).alias("cs_id")
            )
            cred_df = (
                pl.read_csv(
                    cred_file,
                    separator="\t",
                    schema_overrides={
                        "trait": pl.Utf8,
                        "low_purity": pl.Boolean,
                        "cs_size": pl.Int32,
                        "cs_min_r2": pl.Float64,
                    },
                    null_values=["NA"],
                )
                .filter(~pl.col("low_purity"))
                .with_columns(
                    pl.concat_str(
                        [pl.col("region"), pl.col("cs").cast(pl.Utf8)], separator="_"
                    ).alias("cs_id")
                )
            )

            if len(cred_df) == 0:
                raise NoDataException(f"no good credible sets found")

            with warnings.catch_warnings():
                # zero division warnings are given although we handle them
                warnings.filterwarnings("ignore", message="divide by zero encountered in log10")
                merged_df = (
                    snp_df.join(cred_df, on=["trait", "cs_id"], how="inner")
                    .rename(
                        {
                            "chromosome": "chr",
                            "position": "pos",
                            "allele1": "ref",
                            "allele2": "alt",
                            "cs_specific_prob": "pip",
                        }
                    )
                    .with_columns(
                        # use numeric chromosome column
                        pl.col("chr").str.replace("chr", "").replace("X", "23").cast(pl.Int8),
                        # remove predicted.celltype prefix and chr suffix from cell type name
                        pl.when(trait_is_cell_type)
                        .then(
                            pl.col("trait")
                            .str.replace(r"^.*predicted\.celltype\.", "")
                            .str.replace(r"\.chr\d+$|\.chrX$", "")
                            .str.replace(r"\.mean\.inv\.SAIGE", "")
                            .str.replace(r"\.sum\.inv", "")
                            .cast(pl.Utf8)
                        )
                        .otherwise(pl.lit(cell_type))
                        .alias("cell_type"),
                        pl.when(trait_is_cell_type)
                        .then(pl.col("region"))
                        .otherwise(pl.col("trait"))
                        .alias("trait"),
                        pl.lit(dataset).alias("dataset"),
                        pl.lit(data_type).alias("data_type"),
                        # compute -log10(p) from beta and se when p-value underflows
                        # TODO this can be a bit off from the original p-value, e.g. 1.93e-357 can become 1.94e-357
                        # would be better if the original -log10(p) would be in the data file directly to avoid confusion
                        pl.when(pl.col("p").eq(0))
                        .then(
                            (
                                (-log_ndtr(-(pl.col("beta") / pl.col("se")).abs()) - np.log(2))
                                / np.log(10)
                            ).round(4)
                        )
                        .otherwise((-np.log10(pl.col("p")).round(4)))
                        .alias("mlog10p"),
                        # use scientific notation strings for af, beta and se to avoid rounding them to zero
                        pl.col("maf")
                        .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
                        .alias("aaf"), # the column is actually MAF from fine-mapping pipeline so we've lost allele information and need to join AAF from variant annotation file
                        pl.col("beta")
                        .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
                        .alias("beta"),
                        pl.col("se")
                        .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
                        .alias("se"),
                        pl.col("pip").round(4),
                        pl.col("cs_min_r2").round(4),
                    )
                    .select(
                        pl.col("dataset"),
                        pl.col("data_type"),
                        pl.col("trait"),
                        pl.col("cell_type"),
                        pl.col("chr"),
                        pl.col("pos"),
                        pl.col("ref"),
                        pl.col("alt"),
                        pl.col("mlog10p"),
                        pl.col("beta"),
                        pl.col("se"),
                        pl.col("pip"),
                        pl.col("cs_id"),
                        pl.col("cs_size"),
                        pl.col("cs_min_r2"),
                        pl.col("aaf"),
                        *([pl.col("most_severe")] if "most_severe" in snp_df.columns else []),
                        *(
                            [pl.col("gene_most_severe")]
                            if "gene_most_severe" in snp_df.columns
                            else []
                        ),
                    )
                )

            if trait_mapping is not None:
                merged_df = merged_df.join(trait_mapping, on="trait", how="left").rename(
                    {"trait": "trait_original", "mapped_trait": "trait"}
                )
                if not require_trait_mapping:
                    merged_df = merged_df.with_columns(
                        pl.when(pl.col("trait").is_null())
                        .then(pl.col("trait_original"))
                        .otherwise(pl.col("trait"))
                        .alias("trait")
                    )
                null_traits = merged_df.filter(pl.col("trait").is_null())
            else:
                merged_df = merged_df.with_columns(pl.col("trait").alias("trait_original"))
                null_traits = merged_df.filter(False)

            if variant_annotation is not None:
                # drop existing annotation columns if any
                merged_df = merged_df.drop(
                    [col for col in merged_df.columns if col.endswith("most_severe") or col.startswith("aaf")]
                )
                merged_df = (
                    merged_df.with_columns(
                        pl.concat_str(
                            [
                                pl.col("chr"),
                                pl.col("pos").cast(pl.Utf8),
                                pl.col("ref"),
                                pl.col("alt"),
                            ],
                            separator=":",
                        ).alias("variant_id")
                    )
                    .join(
                        variant_annotation,
                        on="variant_id",
                        how="left",
                    )
                    .drop("variant_id")
                )

            return (
                # redefine column order after adding trait_original
                merged_df.select(
                    "dataset",
                    "data_type",
                    "trait",
                    "trait_original",
                    *[
                        col
                        for col in merged_df.columns
                        if col not in ["dataset", "data_type", "trait", "trait_original"]
                    ],
                ).sort("chr", "pos", "ref", "alt", "trait"),
                null_traits.select(
                    "dataset",
                    "data_type",
                    "trait",
                    "trait_original",
                    *[
                        col
                        for col in null_traits.columns
                        if col not in ["dataset", "data_type", "trait", "trait_original"]
                    ],
                ).sort("chr", "pos", "ref", "alt", "trait"),
            )

        dataset = "~{dataset}"
        traits = "~{sep=',' traits}".split(",")
        data_types = "~{sep=',' data_types}".split(",")
        snp_files = "~{sep=',' snp_files}".split(",")
        cred_files = "~{sep=',' cred_files}".split(",")
        trait_is_cell_type = True if "~{trait_is_cell_type}" == "true" else False
        cell_type = "~{cell_type}" if "~{defined(cell_type)}" == "true" else None
        trait_mapping_file = "~{trait_mapping_file}" if "~{defined(trait_mapping_file)}" == "true" else None
        require_trait_mapping = True if "~{require_trait_mapping}" == "true" else False
        variant_annotation_file = "~{variant_annotation_file}" if "~{defined(variant_annotation_file)}" == "true" else None

        if trait_mapping_file is not None:
            trait_mapping = (
                pl.read_csv(trait_mapping_file, separator="\t", has_header=False)
                .rename({"column_1": "trait", "column_2": "mapped_trait"})
                .select("trait", "mapped_trait")
            )
        else:
            trait_mapping = None

        if variant_annotation_file is not None:
            variant_annotation = (
                pl.scan_csv(variant_annotation_file, separator="\t")
                .rename({"#variant": "variant_id"})
                .with_columns(
                    pl.col("AF")
                    .map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8)
                    .alias("aaf"),
                )
                .select("variant_id", "aaf", "most_severe", "gene_most_severe")
                .collect()
            )
        else:
            variant_annotation = None

        all_stats = []

        for i, trait in enumerate(traits):
            try:
                merged_df, null_traits = merge(
                    dataset,
                    data_types[i],
                    snp_files[i],
                    cred_files[i],
                    cell_type,
                    trait_is_cell_type,
                    trait_mapping,
                    require_trait_mapping,
                    variant_annotation,
                )
                if len(null_traits) > 0:
                    null_traits.write_csv(
                        f"individual/{trait}.SUSIE.munged.null_traits.tsv",
                        separator="\t",
                        null_value="NA",
                    )
                    print(
                        f"{i+1}/{len(traits)}: {trait}: {len(set(null_traits.select('trait').to_series().to_list()))} unmapped traits"
                    )

                merged_df.write_csv(
                    f"individual/{trait}.SUSIE.munged.tsv",
                    separator="\t",
                    null_value="NA",
                )

                # calculate stats for each unique trait and write JSON
                for unique_trait in merged_df.filter(pl.col("trait").is_not_null())["trait"].unique():
                    trait_df = merged_df.filter(pl.col("trait") == unique_trait)
                    stats = calculate_stats(trait_df)
                    unique_trait_clean = unique_trait.replace("/", "_") # gene names can contain / which are not great in path names
                    write_stats_json(stats, f"individual/{trait}.{unique_trait_clean}.SUSIE.munged.stats.json")
                    all_stats.append(stats)

                print(f"{i+1}/{len(traits)}: {trait}: OK ({len(merged_df)} rows)")
            except NoDataException as e:
                print(f"{i+1}/{len(traits)}: {trait}: {e}")
            except pl.exceptions.NoDataError as e:
                print(f"{i+1}/{len(traits)}: {trait}: {e}")

        # write aggregate statistics TSV
        with open("credible_set_stats.tsv", "w") as f:
            f.write(get_tsv_header() + "\n")
            for s in all_stats:
                f.write(stats_to_tsv_row(s) + "\n")
        print(f"Wrote aggregate stats for {len(all_stats)} traits")
        EOF

        # merge per-trait files to one file and tabix it
        # sort key is chr, pos, ref, alt, trait
        echo $(date) "merge"
        ls -1 individual/*.SUSIE.munged.tsv | tr '\n' '\0' > merge_these
        # use uniq to remove duplicate header rows
        cat \
        <(echo -n "#") \
        <(sort -m -T . --files0-from=merge_these -k6,6g -k7,7g -k8,8 -k9,9 -k3,3 | uniq) \
        | bgzip > ~{output_file}
        tabix -s6 -b7 -e7 ~{output_file}

        # sanity check that number of rows match the input data
        # because the merging will not work correctly e.g. if column order changes
        n_lines_orig=$(cat individual/*.SUSIE.munged.tsv | grep -Ev "^dataset" | sort -u | wc -l)
        n_lines_output=$(zcat ~{output_file} | tail -n+2 | wc -l)
        if [ "$n_lines_orig" -ne "$n_lines_output" ]; then
            echo "ERROR: number of rows in the output file ($n_lines_output) does not match the input data ($n_lines_orig)" >&2
            exit 1
        fi

        # TODO this doesn't work, Cromwell can't seem to find the files using read_lines, instead using a glob now but that messes up directory structure
        # ls -1 individual/*.SUSIE.munged.tsv > list_of_munged_files

        echo $(date) "done"

    >>>
}
