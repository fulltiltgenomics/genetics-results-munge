#!/usr/bin/env python3

### converts a gencode gff3 file to a tsv file with one row per exon: gene id, transcript id,
### chromosome, exon start/end, strand, exon number/id, gene name/type, transcript type/name,
### the coding portion of that exon (NA where the exon is entirely untranslated), and the
### Ensembl-canonical / MANE Select flags
###
### the chromosome encoding is the one gencode_to_gene_pos_tsv_gff.py already writes
### (X=23, Y=24, M=26), so this file joins to gencode.vNN.annotation.genes.tsv on chrom
###
### the coding portion is carried per exon rather than as separate CDS rows: a gene track draws
### the exon as a box and the translated part of it thicker, so both live on the same row and
### the file stays one row per exon instead of two rows per coding exon
###
### usage, e.g. for gencode v49: python3 gencode_to_exon_tsv_gff.py https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gff3.gz gencode.v49.annotation.exons.tsv

import sys
import gzip
import io
import urllib.request

gencode_gff_file = sys.argv[1]
output_tsv_file = sys.argv[2]

COLUMNS = [
    "gene_id",
    "transcript_id",
    "chrom",
    "exon_start",
    "exon_end",
    "exon_strand",
    "exon_number",
    "exon_id",
    "gene_name",
    "gene_type",
    "transcript_type",
    "transcript_name",
    "cds_start",
    "cds_end",
    "is_canonical",
    "is_mane_select",
]


def attributes(field: str) -> dict[str, str]:
    out = {}
    for pair in field.strip().rstrip(";").split(";"):
        if "=" not in pair:
            continue
        key, _, value = pair.partition("=")
        # `tag` repeats within one attribute field; keep every value
        if key == "tag" and key in out:
            out[key] += "," + value
        else:
            out[key] = value
    return out


def to_chrom(seqid: str) -> str:
    c = seqid.replace("chr", "")
    return {"X": "23", "Y": "24", "M": "26"}.get(c, c)


def open_gff(path: str):
    if path.startswith("http://") or path.startswith("https://"):
        response = urllib.request.urlopen(path)
        return io.TextIOWrapper(gzip.GzipFile(fileobj=response), encoding="utf-8"), response
    return gzip.open(path, "rt"), None


# pass 1: which transcripts are Ensembl canonical / MANE Select, and where their CDS lies.
# the tags are on the transcript line and the CDS on its own lines, but both are wanted on the
# exon rows, so neither can be resolved in a single streaming pass
canonical: set[str] = set()
mane: set[str] = set()
cds_by_transcript: dict[str, list[tuple[int, int]]] = {}

f, response = open_gff(gencode_gff_file)
for line in f:
    if line.startswith("#"):
        continue
    s = line.rstrip("\n").split("\t")
    if s[2] == "transcript":
        attrs = attributes(s[8])
        tags = attrs.get("tag", "").split(",")
        tx = attrs["transcript_id"]
        if "Ensembl_canonical" in tags:
            canonical.add(tx)
        if "MANE_Select" in tags:
            mane.add(tx)
    elif s[2] == "CDS":
        attrs = attributes(s[8])
        cds_by_transcript.setdefault(attrs["transcript_id"], []).append(
            (int(s[3]), int(s[4]))
        )
f.close()
if response is not None:
    response.close()

# pass 2: one row per exon
outfile_handle = open(output_tsv_file, "wt")
outfile_handle.write("\t".join(COLUMNS) + "\n")

f, response = open_gff(gencode_gff_file)
for line in f:
    if line.startswith("#"):
        continue
    s = line.rstrip("\n").split("\t")
    if s[2] != "exon":
        continue
    attrs = attributes(s[8])
    tx = attrs["transcript_id"]
    exon_start, exon_end = int(s[3]), int(s[4])

    # the translated part of THIS exon: a CDS block never spans an intron, so at most one of
    # the transcript's blocks overlaps, and the intersection is the coding sub-interval to draw
    cds_start, cds_end = "NA", "NA"
    for block_start, block_end in cds_by_transcript.get(tx, ()):
        lo, hi = max(exon_start, block_start), min(exon_end, block_end)
        if lo <= hi:
            cds_start, cds_end = str(lo), str(hi)
            break

    outfile_handle.write(
        "\t".join(
            [
                attrs["gene_id"],
                tx,
                to_chrom(s[0]),
                str(exon_start),
                str(exon_end),
                s[6],
                attrs.get("exon_number", "NA"),
                attrs.get("exon_id", "NA"),
                attrs.get("gene_name", "NA"),
                attrs.get("gene_type", "NA"),
                attrs.get("transcript_type", "NA"),
                attrs.get("transcript_name", "NA"),
                cds_start,
                cds_end,
                "1" if tx in canonical else "0",
                "1" if tx in mane else "0",
            ]
        )
        + "\n"
    )

f.close()
if response is not None:
    response.close()
outfile_handle.close()
