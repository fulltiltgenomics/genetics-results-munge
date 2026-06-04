#!/usr/bin/env python3

### converts gencode gff file to tsv file with gene id, chromosome, gene start, gene end, gene strand, gene name, and gene type
### usage, e.g. for gencode v49: python3 gencode_to_gene_pos_tsv.py https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gff3.gz gencode.v49.annotation.genes.tsv

import polars as pl
import sys
import gzip
import io
import urllib.request

gencode_gff_file = sys.argv[1]
output_tsv_file = sys.argv[2]

outfile_handle = open(output_tsv_file, "wt")
outfile_handle.write(
    "gene_id\tchrom\tgene_start\tgene_end\tgene_strand\tgene_name\tgene_type\n"
)

response = None
if gencode_gff_file.startswith("http://") or gencode_gff_file.startswith("https://"):
    response = urllib.request.urlopen(gencode_gff_file)
    f = io.TextIOWrapper(gzip.GzipFile(fileobj=response), encoding="utf-8")
else:
    f = gzip.open(gencode_gff_file, "rt")

for line in f:
    if line.startswith("#"):
        continue
    s = line.strip().split("\t")
    if s[2] != "gene":
        continue
    chrom = (
        s[0].replace("chr", "").replace("X", "23").replace("Y", "24").replace("M", "26")
    )
    gene_start = s[3]
    gene_end = s[4]
    gene_strand = s[6]
    fields = s[8].split(";")
    gene_id = [f for f in fields if f.startswith("gene_id=")][0].split("=")[1]
    gene_name = [f for f in fields if f.startswith("gene_name=")][0].split("=")[1]
    gene_type = [f for f in fields if f.startswith("gene_type=")][0].split("=")[1]
    outfile_handle.write(
        f"{gene_id}\t{chrom}\t{gene_start}\t{gene_end}\t{gene_strand}\t{gene_name}\t{gene_type}\n"
    )

f.close()
if response is not None:
    response.close()

outfile_handle.close()
