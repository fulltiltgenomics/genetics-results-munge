#!/usr/bin/env python3
"""
Convert Hail MatrixTable gene burden results to TSV format.
Filters by p-value threshold.
Joins and uses gene positions from gencode.

Output is NOT sorted: the unfiltered export is one row per gene x annotation x
phenotype (343M rows for the 2022 release), and ordering that in Hail is a full
shuffle of ~65 GB. convert_genebass_gene_results.sh sorts downstream instead,
where the work is per trait (~76k rows) or on the significant-hits subset.
Give --output a .bgz extension to have Hail block-gzip the export directly, so
the plain TSV never lands on disk.
"""

import hail as hl
import argparse
import gzip
import io
import urllib.request

TMP_DIR = "/mnt/disks/data/hail_tmp" # need some space for sorting
GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gff3.gz"


def load_gencode_gene_positions():
    """Load gene positions from gencode v35 gff3 file, returns a Hail table keyed by gene_name."""
    print(f"Loading gencode data from {GENCODE_URL}...")
    genes = []
    response = urllib.request.urlopen(GENCODE_URL)
    f = io.TextIOWrapper(gzip.GzipFile(fileobj=response), encoding="utf-8")

    for line in f:
        if line.startswith("#"):
            continue
        s = line.strip().split("\t")
        if s[2] != "gene":
            continue
        chrom = s[0].replace("chr", "").replace("X", "23").replace("Y", "24").replace("M", "26")
        chrom = int(chrom)
        gene_start = int(s[3])
        gene_end = int(s[4])
        fields = s[8].split(";")
        gene_name = [field for field in fields if field.startswith("gene_name=")][0].split("=")[1]
        genes.append({"gene_name": gene_name, "gencode_chr": chrom, "gencode_start": gene_start, "gencode_end": gene_end})

    f.close()
    response.close()

    print(f"Loaded {len(genes)} genes from gencode")
    gencode_ht = hl.Table.parallelize(genes)
    gencode_ht = gencode_ht.key_by("gene_name")
    return gencode_ht


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
        description="Convert Hail gene burden results to TSV"
    )
    parser.add_argument(
        "--input", default="/mnt/disks/data/results.mt", help="Input MatrixTable path"
    )
    parser.add_argument(
        "--output", default="gene_burden_results.tsv", help="Output TSV file path"
    )
    parser.add_argument(
        "--pval-threshold",
        type=float,
        default=1,
        help="P-value threshold, keeps Pvalue_Burden <= threshold (default: 1, no filtering)",
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
    args = parser.parse_args()

    blacklist = set()
    if args.pheno_blacklist:
        blacklist = load_blacklist(args.pheno_blacklist)
        print(f"Loaded {len(blacklist)} blacklisted phenotypes")

    hl.init(log="/tmp/hail.log", tmp_dir=TMP_DIR, local_tmpdir=TMP_DIR, spark_conf={"spark.local.dir": TMP_DIR}, quiet=True)

    gencode_ht = load_gencode_gene_positions()

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

    # use gene positions from gencode as the positions in genebass data are something else than regular gene start and end.
    # joined on the MT rows (~76k) rather than on the flattened entries (~343M), which would shuffle the whole export
    print("Joining with gencode data for gene positions...")
    gencode_row = gencode_ht[mt.gene_symbol]
    mt = mt.annotate_rows(
        gene_chr=hl.int32(gencode_row.gencode_chr),
        gene_start_pos=hl.int32(gencode_row.gencode_start),
        gene_end_pos=hl.int32(gencode_row.gencode_end),
    )

    print("Flattening MatrixTable to entries table...")
    entries = mt.entries()

    # Add computed columns
    print("Adding computed columns...")
    entries = entries.annotate(
        mlog10p_burden=hl.if_else(
            hl.is_defined(entries.Pvalue_Burden),
            hl.if_else(
                entries.Pvalue_Burden == 0,
                324.0,
                hl.if_else(
                    entries.Pvalue_Burden > 0,
                    -hl.log10(entries.Pvalue_Burden),
                    hl.missing(hl.tfloat64),
                ),
            ),
            hl.missing(hl.tfloat64),
        ),
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
        beta=entries.BETA_Burden,
        se=entries.SE_Burden,
    )

    print(f"Filtering by burden p-value <= {args.pval_threshold}...")
    entries = entries.filter(
        hl.is_defined(entries.Pvalue_Burden)
        & (entries.Pvalue_Burden <= args.pval_threshold)
    )

    # drop keys so we can freely select/rename columns
    entries = entries.key_by()

    #   ROW FIELDS (Gene-level information)

    #   - gene_id - Ensembl gene ID
    #   - gene_symbol - Gene symbol/name
    #   - annotation - Variant annotation type (pLoF, missense|LC, synonymous, etc.)
    #   - interval - Genomic interval (chromosome and position range)
    #   - markerIDs - List of variant IDs
    #   - markerAFs - List of allele frequencies
    #   - total_variants - Total number of variants in gene
    #   - Nmarker_MACCate_1 through Nmarker_MACCate_8 - Number of markers in each Minor Allele Count (MAC) category

    #   COLUMN FIELDS (Trait/Phenotype information)

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

    #   ENTRY FIELDS (Gene-Trait association results)

    #   - Pvalue - SKAT-O p-value
    #   - Pvalue_Burden - Burden test p-value
    #   - Pvalue_SKAT - SKAT test p-value
    #   - BETA_Burden - Effect size from burden test
    #   - SE_Burden - Standard error from burden test
    #   - Pvalue.NA, Pvalue_Burden.NA, Pvalue_SKAT.NA, BETA_Burden.NA, SE_Burden.NA - NA-adjusted versions of these statistics
    #   - total_variants_pheno - Number of variants tested for this specific phenotype

    # select and rename columns
    entries = entries.select(
        trait=entries.trait,
        gene=entries.gene_symbol,
        gene_id=entries.gene_id,
        gene_chr=entries.gene_chr,
        gene_start_pos=entries.gene_start_pos,
        gene_end_pos=entries.gene_end_pos,
        annotation=entries.annotation,
        mlog10p_burden=entries.mlog10p_burden,
        beta=entries.beta,
        se=entries.se,
        total_variants=entries.total_variants,
        total_variants_pheno=entries.total_variants_pheno,
        n_cases=entries.n_cases,
        n_controls=entries.n_controls,
        trait_original=entries.trait_original,
    )

    # convert missing values and empty strings to "NA" strings for export
    entries = entries.select(
        dataset=hl.literal("genebass"),
        # fall back to trait_original when no readable description is available
        trait=hl.if_else(
            hl.is_missing(entries.trait) | (entries.trait == ""),
            hl.if_else(
                hl.is_missing(entries.trait_original) | (entries.trait_original == ""),
                "NA",
                entries.trait_original,
            ),
            entries.trait,
        ),
        gene=hl.if_else(hl.is_missing(entries.gene) | (entries.gene == ""), "NA", entries.gene),
        gene_id=hl.if_else(hl.is_missing(entries.gene_id) | (entries.gene_id == ""), "NA", entries.gene_id),
        gene_chr=hl.bind(
            lambda s: hl.if_else(s == "", "NA", s),
            hl.or_else(hl.str(entries.gene_chr), "")
        ),
        gene_start_pos=hl.bind(
            lambda s: hl.if_else(s == "", "NA", s),
            hl.or_else(hl.str(entries.gene_start_pos), "")
        ),
        gene_end_pos=hl.bind(
            lambda s: hl.if_else(s == "", "NA", s),
            hl.or_else(hl.str(entries.gene_end_pos), "")
        ),
        annotation=hl.if_else(hl.is_missing(entries.annotation) | (entries.annotation == ""), "NA", entries.annotation),
        mlog10p_burden=hl.bind(
            lambda s: hl.if_else(s == "", "NA", s),
            hl.or_else(hl.str(entries.mlog10p_burden), "")
        ),
        beta=hl.bind(
            lambda s: hl.if_else(s == "", "NA", s),
            hl.or_else(hl.str(entries.beta), "")
        ),
        se=hl.bind(
            lambda s: hl.if_else(s == "", "NA", s),
            hl.or_else(hl.str(entries.se), "")
        ),
        total_variants=hl.bind(
            lambda s: hl.if_else(s == "", "NA", s),
            hl.or_else(hl.str(entries.total_variants), "")
        ),
        total_variants_pheno=hl.bind(
            lambda s: hl.if_else(s == "", "NA", s),
            hl.or_else(hl.str(entries.total_variants_pheno), "")
        ),
        n_cases=hl.bind(
            lambda s: hl.if_else(s == "", "NA", s),
            hl.or_else(hl.str(entries.n_cases), "")
        ),
        n_controls=hl.bind(
            lambda s: hl.if_else(s == "", "NA", s),
            hl.or_else(hl.str(entries.n_controls), "")
        ),
        trait_original=hl.if_else(hl.is_missing(entries.trait_original) | (entries.trait_original == ""), "NA", entries.trait_original),
        flags=hl.literal("NA"),
    )
    
    entries = entries.rename({"dataset": "#dataset"})

    print(f"Exporting to {args.output}...")
    entries.export(args.output)

    # no count() here: it would re-run the whole pipeline for a number the
    # downstream split already reports
    print(f"\nDone! Exported to {args.output}")
    print(f"Filtered by p-value <= {args.pval_threshold}")


if __name__ == "__main__":
    main()
