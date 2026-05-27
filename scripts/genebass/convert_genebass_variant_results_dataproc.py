#!/usr/bin/env python3
"""
Convert Hail MatrixTable variant-level results to TSV format.

Reads from GCS bucket (requester pays), applies quality filters
(AF*n_cases >= 5, annotation merge, phenotype blacklist), and exports TSV.

To start a cluster:
  hailctl dataproc start genebass --region europe-west1 --zone europe-west1-b \
    --num-workers 10 --max-idle 30m \
    --subnet projects/phewas-development/regions/europe-west1/subnetworks/default \
    --requester-pays-allow-all

To run the job:
  hailctl dataproc submit genebass --region europe-west1 \
    convert_genebass_variant_results_dataproc.py \
    --input gs://ukbb-exome-public/500k/results/variant_results.mt \
    --output gs://finngen-commons/genebass/variant_results.tsv.bgz \
    --pheno-blacklist gs://finngen-commons/genebass/pheno_blacklist.txt

To stop the cluster:
  hailctl dataproc stop genebass --region europe-west1
"""

import hail as hl
import argparse


def load_blacklist(path):
    """Load phenotype blacklist from a text file (one trait per line)."""
    if path.startswith("gs://"):
        import subprocess
        result = subprocess.run(
            ["gsutil", "cat", path], capture_output=True, text=True, check=True
        )
        lines = result.stdout.strip().split("\n")
    else:
        with open(path) as f:
            lines = f.read().strip().split("\n")
    return set(line.strip() for line in lines if line.strip())


def main():
    parser = argparse.ArgumentParser(
        description="Convert Hail variant-level results to TSV"
    )
    parser.add_argument(
        "--input",
        default="gs://ukbb-exome-public/500k/results/variant_results.mt",
        help="Input MatrixTable path",
    )
    parser.add_argument(
        "--output",
        default="gs://finngen-commons/genebass/variant_results.tsv.bgz",
        help="Output TSV file path",
    )
    parser.add_argument(
        "--pval-threshold",
        type=float,
        default=1e-3,
        help="P-value threshold for pre-filtering (default: 1e-3)",
    )
    parser.add_argument(
        "--pheno-blacklist",
        help="Path to phenotype blacklist file (one trait per line, local or gs://)",
    )
    parser.add_argument(
        "--preview",
        action="store_true",
        help="Preview schema and sample data without writing output",
    )
    parser.add_argument(
        "--test-region",
        help="Restrict to a genomic region for testing, e.g. 'chr22:20000000-21000000'",
    )
    parser.add_argument(
        "--test-n-traits",
        type=int,
        help="Restrict to first N traits for testing",
    )
    args = parser.parse_args()

    blacklist = set()
    if args.pheno_blacklist:
        blacklist = load_blacklist(args.pheno_blacklist)
        print(f"Loaded {len(blacklist)} blacklisted phenotypes")

    hl.init(
        spark_conf={
            "spark.hadoop.fs.gs.requester.pays.mode": "AUTO",
            "spark.hadoop.fs.gs.requester.pays.buckets": "ukbb-exome-public",
            "spark.hadoop.fs.gs.requester.pays.project.id": "phewas-development",
        },
    )

    print(f"Reading MatrixTable from {args.input}...")
    mt = hl.read_matrix_table(args.input)

    if args.test_region:
        interval = hl.parse_locus_interval(args.test_region, reference_genome="GRCh38")
        mt = mt.filter_rows(interval.contains(mt.locus))
        print(f"Test mode: restricted to {args.test_region}")

    if args.test_n_traits:
        mt = mt.add_col_index("_col_idx")
        mt = mt.filter_cols(mt._col_idx < args.test_n_traits)
        mt = mt.drop("_col_idx")
        print(f"Test mode: restricted to first {args.test_n_traits} traits")

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

    # drop synonymous and undefined annotations, keep original pLoF/missense/LC
    mt = mt.filter_rows(
        hl.is_defined(mt.annotation) & (mt.annotation != "synonymous")
    )

    # filter out blacklisted phenotypes
    if blacklist:
        # build trait key from column fields to match blacklist format
        mt = mt.annotate_cols(
            _trait_key=hl.delimit(
                [mt.trait_type, mt.phenocode, mt.pheno_sex, mt.coding, mt.modifier],
                "_",
            )
        )
        mt = mt.filter_cols(~hl.literal(blacklist).contains(mt._trait_key))

    print("Flattening MatrixTable to entries table...")
    entries = mt.entries()

    # filter: AF * n_cases >= 5 (remove false positives)
    entries = entries.filter(
        hl.is_defined(entries.AF)
        & hl.is_defined(entries.n_cases)
        & (entries.AF * entries.n_cases >= 5)
    )

    print(f"Filtering by p-value < {args.pval_threshold}...")
    entries = entries.filter(
        hl.is_defined(entries.Pvalue) & (entries.Pvalue < args.pval_threshold)
    )

    print("Adding computed columns...")
    entries = entries.annotate(
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
        mlog10p=hl.if_else(
            hl.is_defined(entries.Pvalue)
            & hl.is_defined(entries.BETA)
            & hl.is_defined(entries.SE),
            hl.if_else(
                entries.Pvalue == 0,
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
        trait_original=hl.delimit(
            [
                entries.trait_type,
                entries.phenocode,
                entries.pheno_sex,
                entries.coding,
                entries.modifier,
            ],
            "_",
        ),
        # readable trait name from description and coding_description
        trait=hl.bind(
            lambda desc, coding: hl.if_else(
                hl.is_defined(desc) & (desc != "") & hl.is_defined(coding) & (coding != ""),
                desc + " | " + coding,
                hl.if_else(
                    hl.is_defined(desc) & (desc != ""),
                    desc,
                    hl.if_else(
                        hl.is_defined(coding) & (coding != ""),
                        coding,
                        hl.missing(hl.tstr),
                    ),
                ),
            ),
            entries.description,
            entries.coding_description,
        ),
    )

    entries = entries.key_by()

    # round mlog10p to 4 decimal places
    entries = entries.annotate(
        mlog10p=hl.if_else(
            hl.is_defined(entries.mlog10p),
            hl.float64(hl.format("%.4f", entries.mlog10p)),
            hl.missing(hl.tfloat64),
        )
    )

    # format and select output columns
    entries = entries.select(
        **{
            "#dataset": hl.literal("genebass"),
        },
        chr=hl.or_else(hl.str(entries.chr), "NA"),
        pos=hl.or_else(hl.str(entries.pos), "NA"),
        ref=hl.or_else(entries.ref, "NA"),
        alt=hl.or_else(entries.alt, "NA"),
        gene=hl.or_else(entries.gene, "NA"),
        annotation=hl.or_else(entries.annotation, "NA"),
        mlog10p=hl.or_else(hl.str(entries.mlog10p), "NA"),
        beta=hl.or_else(hl.format("%.3e", entries.beta), "NA"),
        se=hl.or_else(hl.format("%.3e", entries.se), "NA"),
        af_overall=hl.or_else(hl.format("%.3e", entries.af_overall), "NA"),
        af_cases=hl.or_else(hl.format("%.3e", entries.af_cases), "NA"),
        af_controls=hl.or_else(hl.format("%.3e", entries.af_controls), "NA"),
        ac=hl.or_else(hl.str(entries.ac), "NA"),
        an=hl.or_else(hl.str(entries.an), "NA"),
        n_cases=hl.or_else(hl.str(entries.n_cases), "NA"),
        n_controls=hl.or_else(hl.str(entries.n_controls), "NA"),
        trait=hl.or_else(entries.trait, "NA"),
        trait_original=hl.or_else(entries.trait_original, "NA"),
    )

    print("Sorting results...")
    entries = entries.order_by(
        entries.chr,
        entries.pos,
        entries.ref,
        entries.alt,
    )

    print(f"Exporting to {args.output}...")
    entries.export(args.output)

    n_results = entries.count()
    print(f"\nDone! Exported {n_results} results to {args.output}")


if __name__ == "__main__":
    main()
