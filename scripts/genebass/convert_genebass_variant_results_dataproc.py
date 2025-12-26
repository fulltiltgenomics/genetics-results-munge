#!/usr/bin/env python3
"""
Convert Hail MatrixTable variant-level results to TSV format.
Filters by p-value threshold.
Reads from GCS bucket (requester pays).
"""

import hail as hl
import argparse


# to start a cluster:
# (only use a high number of workers if sure that this works)
# better use a US region/zone to avoid egress costs but this requires permission I didn't have
# hailctl dataproc start genebass --region europe-west1 --zone europe-west1-b --num-workers 10 --max-idle 30m --subnet projects/phewas-development/regions/europe-west1/subnetworks/default --requester-pays-allow-all
#
# to run the job:
# hailctl dataproc submit genebass --region europe-west1 convert_genebass_variant_results_dataproc.py --input gs://ukbb-exome-public/500k/results/variant_results.mt --output gs://finngen-commons/genebass/variant_results_mlog10p3.tsv.bgz --pval-threshold 1e-3
#
# to stop the cluster:
# (don't rely on the automatic idle shutdown)
# hailctl dataproc stop genebass --region europe-west1


def main():
    parser = argparse.ArgumentParser(
        description="Convert Hail variant-level results to TSV"
    )
    parser.add_argument(
        "--input",
        default="gs://ukbb-exome-public/500k/results/variant_results.mt",
        help="Input MatrixTable path (default: GCS bucket path)",
    )
    parser.add_argument(
        "--output",
        default="gs://finngen-commons/genebass/variant_results_mlog10p3.tsv.bgz",
        help="Output TSV file path (default: gs://finngen-commons/genebass/variant_results_mlog10p3.tsv.bgz)",
    )
    parser.add_argument(
        "--pval-threshold",
        type=float,
        default=1e-3,
        help="P-value threshold (default: 1e-3)",
    )
    parser.add_argument(
        "--preview",
        action="store_true",
        help="Preview schema and sample data without writing output",
    )
    args = parser.parse_args()

    hl.init(
        spark_conf={
            "spark.hadoop.fs.gs.requester.pays.mode": "AUTO",
            "spark.hadoop.fs.gs.requester.pays.buckets": "ukbb-exome-public",
            "spark.hadoop.fs.gs.requester.pays.project.id": "phewas-development",
        },
    )

    print(f"Reading MatrixTable from {args.input}...")
    mt = hl.read_matrix_table(args.input)

    if args.preview:
        print("\n=== MatrixTable Schema ===")
        mt.describe()
        print("\n=== Sample Rows ===")
        mt.rows().show(5)
        print("\n=== Sample Columns ===")
        mt.cols().show(5)
        print("\n=== Sample Entries ===")
        mt.entries().show(5)
        return

    print("Flattening MatrixTable to entries table...")
    entries = mt.entries()

    print("Adding computed columns...")
    entries = entries.annotate(
        # Extract chromosome from locus and convert to integer
        # Remove "chr" prefix if present, convert X->23, Y->24, MT->26
        chr=hl.bind(
            lambda chrom: hl.switch(chrom)
            .when("X", 23)
            .when("Y", 24)
            .when("MT", 26)
            .when("M", 26)
            .default(hl.int32(chrom)),
            hl.if_else(
                entries.locus.contig.startswith("chr"),
                entries.locus.contig[3:],
                entries.locus.contig,
            ),
        ),
        pos=entries.locus.position,
        ref=entries.alleles[0],
        alt=entries.alleles[1],
        rsid=entries.markerID,
        # compute -log10(p-value) with proper underflow handling
        mlog10p=hl.if_else(
            hl.is_defined(entries.Pvalue)
            & hl.is_defined(entries.BETA)
            & hl.is_defined(entries.SE),
            hl.if_else(
                entries.Pvalue == 0,
                # underflow case: compute from z-score
                hl.bind(
                    lambda z: hl.bind(
                        lambda log_phi: -(hl.log(2.0) + log_phi) / hl.log(10.0),
                        hl.log(hl.pnorm(-hl.abs(z))),
                    ),
                    entries.BETA / entries.SE,
                ),
                hl.if_else(
                    entries.Pvalue > 0,
                    -hl.log10(entries.Pvalue),
                    hl.missing(hl.tfloat64),
                ),
            ),
            hl.missing(hl.tfloat64),
        ),
        af_overall=entries.AF,
        af_cases=entries["AF.Cases"],
        af_controls=entries["AF.Controls"],
        ac=entries.AC,
        an=entries.call_stats.AN,
        beta=entries.BETA,
        se=entries.SE,
        trait=hl.delimit(
            [
                entries.trait_type,
                entries.phenocode,
                entries.pheno_sex,
                entries.coding,
                entries.modifier,
            ],
            "_",
        ),
    )

    print(f"Filtering by p-value < {args.pval_threshold}...")
    entries = entries.filter(
        hl.is_defined(entries.Pvalue) & (entries.Pvalue < args.pval_threshold)
    )

    entries = entries.key_by()

    #   ROW FIELDS (Variant-level information)
    #
    #   - locus - Genomic locus (chr:pos)
    #   - alleles - Array of alleles [ref, alt]
    #   - markerID - Variant ID (rsID or chr:pos_ref/alt)
    #   - gene - Gene symbol
    #   - annotation - Variant annotation type (LC, synonymous, missense, etc.)
    #   - call_stats - Population-level statistics (AC, AF, AN, homozygote_count)

    #   COLUMN FIELDS (Trait/Phenotype information)
    #
    #   - n_cases - Number of cases
    #   - n_controls - Number of controls
    #   - heritability - Estimated heritability
    #   - saige_version - SAIGE software version
    #   - inv_normalized - Whether trait was inverse-normalized
    #   - trait_type - Type of trait (categorical, continuous, etc.)
    #   - phenocode - Phenotype code
    #   - pheno_sex - Sex group (both_sexes, males, females)
    #   - coding - Coding/subcategory for the trait
    #   - modifier - Additional modifier
    #   - n_cases_defined - Number of defined cases
    #   - n_cases_both_sexes, n_cases_females, n_cases_males - Sex-stratified case counts
    #   - description - Trait description
    #   - description_more - Extended description
    #   - coding_description - Description of the coding
    #   - category - Phenotype category

    #   ENTRY FIELDS (Variant-Trait association results)
    #
    #   - AC - Allele count in tested samples
    #   - AF - Allele frequency in tested samples
    #   - BETA - Effect size
    #   - SE - Standard error
    #   - AF.Cases - Allele frequency in cases
    #   - AF.Controls - Allele frequency in controls
    #   - Pvalue - Association p-value

    entries = entries.select(
        chr=entries.chr,
        pos=entries.pos,
        ref=entries.ref,
        alt=entries.alt,
        rsid=entries.rsid,
        gene=entries.gene,
        annotation=entries.annotation,
        mlog10p=entries.mlog10p,
        beta=entries.beta,
        se=entries.se,
        af_overall=entries.af_overall,
        af_cases=entries.af_cases,
        af_controls=entries.af_controls,
        ac=entries.ac,
        an=entries.an,
        n_cases=entries.n_cases,
        n_controls=entries.n_controls,
        heritability=entries.heritability,
        trait=entries.trait,
        description=entries.description,
        coding_description=entries.coding_description,
        category=entries.category,
    )

    print("Sorting results...")
    entries = entries.order_by(
        entries.chr,
        entries.pos,
        entries.ref,
        entries.alt,
    )

    # round mlog10p to 4 decimal places
    entries = entries.annotate(
        mlog10p=hl.if_else(
            hl.is_defined(entries.mlog10p),
            hl.float64(hl.format("%.4f", entries.mlog10p)),
            hl.missing(hl.tfloat64),
        )
    )

    # convert missing values to "NA" strings for export
    entries = entries.select(
        chr=hl.or_else(hl.str(entries.chr), "NA"),
        pos=hl.or_else(hl.str(entries.pos), "NA"),
        ref=hl.or_else(entries.ref, "NA"),
        alt=hl.or_else(entries.alt, "NA"),
        rsid=hl.or_else(entries.rsid, "NA"),
        gene=hl.or_else(entries.gene, "NA"),
        annotation=hl.or_else(entries.annotation, "NA"),
        mlog10p=hl.or_else(hl.str(entries.mlog10p), "NA"),
        beta=hl.or_else(hl.str(entries.beta), "NA"),
        se=hl.or_else(hl.str(entries.se), "NA"),
        af_overall=hl.or_else(hl.str(entries.af_overall), "NA"),
        af_cases=hl.or_else(hl.str(entries.af_cases), "NA"),
        af_controls=hl.or_else(hl.str(entries.af_controls), "NA"),
        ac=hl.or_else(hl.str(entries.ac), "NA"),
        an=hl.or_else(hl.str(entries.an), "NA"),
        n_cases=hl.or_else(hl.str(entries.n_cases), "NA"),
        n_controls=hl.or_else(hl.str(entries.n_controls), "NA"),
        heritability=hl.or_else(hl.str(entries.heritability), "NA"),
        trait=hl.or_else(entries.trait, "NA"),
        description=hl.or_else(entries.description, "NA"),
        coding_description=hl.or_else(entries.coding_description, "NA"),
        category=hl.or_else(entries.category, "NA"),
    )

    print(f"Exporting to {args.output}...")
    entries.export(args.output)

    n_results = entries.count()
    print(f"\nDone! Exported {n_results} significant results to {args.output}")
    print(f"Filtered by p-value < {args.pval_threshold}")


if __name__ == "__main__":
    main()
