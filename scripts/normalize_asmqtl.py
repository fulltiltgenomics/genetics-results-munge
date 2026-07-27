#!/usr/bin/env python3
"""Normalize ASM-QTL indel alleles using bcftools norm (left-align + trim).

Produces a TSV mapping file with columns:
    chrom, pos, ref, alt, norm_pos, norm_ref, norm_alt

Only indels that changed are included. The munge script reads this mapping
and applies it during munging.

Usage:
    python scripts/normalize_asmqtl.py \
        --input /mnt/disks/data/CPGMethylation/Publish/Data-S1.tab \
        --ref /mnt/disks/data/ref/GRCh38.fa
"""

import argparse
import subprocess
import tempfile

import polars as pl


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Path to Data-S1.tab or Data-S3.tab")
    parser.add_argument("--ref", required=True, help="Path to GRCh38 reference FASTA (indexed)")
    parser.add_argument("--output", help="Output mapping TSV (default: <input>.norm_mapping.tsv)")
    return parser.parse_args()


def _detect_chr_prefix(ref_path: str) -> bool:
    fai = ref_path + ".fai"
    with open(fai) as f:
        first_chr = f.readline().split("\t")[0]
    return first_chr.startswith("chr")


def main():
    args = parse_args()

    input_path = args.input
    output_path = args.output or input_path.replace(".tab", ".norm_mapping.tsv")
    ref_has_chr = _detect_chr_prefix(args.ref)

    print(f"Reading {input_path}...")
    df = pl.read_csv(input_path, separator="\t", null_values=["NA"],
                     infer_schema_length=0)
    print(f"  {df.height} rows")

    variants = df.select("Chrom", "SeqVariant_start", "SeqVariant_ref", "SeqVariant_alt", "SeqVariant_vartype").unique()
    print(f"  {variants.height} unique variants")

    # only normalize standard indels (skip SNVs, SVs, microsatellites)
    is_indel = pl.col("SeqVariant_ref").str.len_bytes() != pl.col("SeqVariant_alt").str.len_bytes()
    acgtn = pl.col("SeqVariant_ref").str.contains(r"^[ACGTNacgtn]+$") & pl.col("SeqVariant_alt").str.contains(r"^[ACGTNacgtn]+$")
    is_standard_indel = is_indel & acgtn & (pl.col("SeqVariant_vartype") == "Indel")
    indels = variants.filter(is_standard_indel)
    n_sv = variants.filter(is_indel & ~acgtn).height
    n_micro = variants.filter(is_indel & acgtn & (pl.col("SeqVariant_vartype") != "Indel")).height
    print(f"  {indels.height} indels to normalize ({n_sv} structural variants, {n_micro} microsatellites skipped)")

    if indels.height == 0:
        print("No indels to normalize.")
        # write empty mapping file with header
        with open(output_path, "w") as f:
            f.write("chrom\tpos\tref\talt\tnorm_pos\tnorm_ref\tnorm_alt\n")
        return

    # unique ID per indel for tracking through bcftools
    indels = indels.with_columns(
        (pl.col("Chrom") + ":" + pl.col("SeqVariant_start") + ":"
         + pl.col("SeqVariant_ref") + ":" + pl.col("SeqVariant_alt")).alias("_vid"),
    )

    tmpdir = tempfile.mkdtemp()
    vcf_in = f"{tmpdir}/indels.vcf"
    vcf_out = f"{tmpdir}/indels.norm.vcf"

    contig_lines = []
    with open(args.ref + ".fai") as fai:
        for line in fai:
            name, length = line.split("\t")[:2]
            contig_lines.append(f"##contig=<ID={name},length={length}>")

    with open(vcf_in, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for cl in contig_lines:
            f.write(cl + "\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in indels.iter_rows(named=True):
            chrom = row["Chrom"]
            if ref_has_chr and not chrom.startswith("chr"):
                chrom = "chr" + chrom
            elif not ref_has_chr and chrom.startswith("chr"):
                chrom = chrom[3:]
            f.write(f"{chrom}\t{row['SeqVariant_start']}\t{row['_vid']}\t{row['SeqVariant_ref']}\t{row['SeqVariant_alt']}\t.\t.\t.\n")

    print("Running bcftools norm...")
    result = subprocess.run(
        ["bcftools", "norm", "-f", args.ref, "-o", vcf_out, vcf_in],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"bcftools norm stderr:\n{result.stderr}")
        raise RuntimeError("bcftools norm failed")
    for line in result.stderr.splitlines():
        print(f"  {line}")

    # read normalized VCF, keyed by variant ID
    norm_map = {}
    with open(vcf_out) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            norm_map[fields[2]] = (fields[1], fields[3], fields[4])

    # write mapping: only rows that actually changed
    n_changed = 0
    with open(output_path, "w") as out:
        out.write("chrom\tpos\tref\talt\tnorm_pos\tnorm_ref\tnorm_alt\n")
        for row in indels.iter_rows(named=True):
            vid = row["_vid"]
            norm = norm_map.get(vid)
            if norm is None:
                continue
            norm_pos, norm_ref, norm_alt = norm
            orig_pos = row["SeqVariant_start"]
            orig_ref = row["SeqVariant_ref"]
            orig_alt = row["SeqVariant_alt"]
            if orig_pos != norm_pos or orig_ref != norm_ref or orig_alt != norm_alt:
                out.write(f"{row['Chrom']}\t{orig_pos}\t{orig_ref}\t{orig_alt}\t{norm_pos}\t{norm_ref}\t{norm_alt}\n")
                n_changed += 1

    print(f"  {n_changed} indels changed by normalization")
    print(f"  mapping written to {output_path}")


if __name__ == "__main__":
    main()
