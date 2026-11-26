#!/usr/bin/env python3
"""Lift GRCh37 variants (and plain intervals) to GRCh38 with UCSC liftOver, then hold every
lifted variant to the GRCh38 reference sequence.

liftOver moves coordinates only: it knows nothing of alleles. So a variant is lifted as the
interval its REF allele covers and then checked here:
  - the interval must land on the same chromosome with its length unchanged
  - where the chain maps it to the minus strand, REF and ALT are reverse-complemented
  - REF must then equal the GRCh38 sequence at the lifted interval. A variant whose REF does not
    match is dropped, including the ones whose ALT is the GRCh38 reference: recovering those
    means swapping REF and ALT, and with them every per-allele count the caller carries
  - indels are left-aligned in GRCh38 and anchored on the base before them. On the minus strand
    the reverse complement moves the shared base to the other end, so this is not optional there
  - variants that land on one GRCh38 position with the same alleles are all dropped: nothing says
    which of them the GRCh38 variant is
The mitochondrion is not lifted: the hg19 chain's chrM is not the sequence GRCh37 (b37) calls MT.

Every drop is counted by reason in the `counts` dict the caller passes.

CROSS-CHECK AGAINST PICARD (2026-11-26, the 8,865,268 variants of the SCHEMA 2021 release, same
chain and the same hg38.fa for both): Picard 3.5.0 LiftoverVcf with default settings agreed on
position and alleles for every variant both methods kept, except 32 minus-strand indels that
Picard wrote one repeat unit to the right of the left-aligned form written here (same
haplotype, checked against the reference). Picard also kept 366 variants that land on another
contig, which are dropped here. Rerun the comparison with
  java -Xmx4g -jar picard.jar LiftoverVcf -I sites.hg19.vcf.gz -O picard.hg38.vcf.gz \
      --CHAIN hg19ToHg38.over.chain.gz --REJECT picard.reject.vcf.gz -R hg38.fa
on a sites-only VCF with chr-prefixed contigs; Picard holds the whole target genome in memory,
so it needs about 4 GB of heap where this needs the chain and a memory map.

Chain and binary (UCSC, not vendored here):
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
  https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
The fasta is an uncompressed, .fai-indexed GRCh38 with chr-prefixed names (UCSC hg38.fa).
"""

import mmap
import subprocess
import tempfile
from pathlib import Path

import polars as pl

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


class Fasta:
    """Upper-cased sequence from an uncompressed fasta through its .fai index."""

    def __init__(self, path: str | Path):
        self._index = {}
        with open(f"{path}.fai") as fai:
            for line in fai:
                name, length, offset, bases, width = line.split("\t")[:5]
                self._index[name] = (int(length), int(offset), int(bases), int(width))
        self._file = open(path, "rb")
        self._map = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """0-based half-open; empty when the interval is not inside the contig."""
        length, offset, bases, width = self._index[chrom]
        if start < 0 or end > length:
            return ""
        a = offset + start // bases * width + start % bases
        b = offset + end // bases * width + end % bases
        return self._map[a:b].replace(b"\n", b"").decode().upper()


def left_align(fasta: Fasta, chrom: str, pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Left-align an indel (1-based pos) and leave exactly one anchor base before it."""
    while ref and alt and ref[-1] == alt[-1]:
        if len(ref) == 1 or len(alt) == 1:
            base = fasta.fetch(chrom, pos - 2, pos - 1)
            if not base:
                break
            ref, alt, pos = base + ref, base + alt, pos - 1
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def run_liftover(liftover_bin: str, chain: str, bed_rows: list[str], work: Path) -> tuple[list[list[str]], int]:
    """liftOver over BED6 lines; the lifted rows and the number liftOver could not map."""
    bed_in, bed_out, unmapped = work / "in.bed", work / "out.bed", work / "unmapped.bed"
    bed_in.write_text("".join(bed_rows))
    subprocess.run([liftover_bin, str(bed_in), chain, str(bed_out), str(unmapped)],
                   check=True, stderr=subprocess.DEVNULL)
    lifted = [line.rstrip("\n").split("\t") for line in bed_out.open()]
    return lifted, len(bed_rows) - len(lifted)


def _count(counts: dict, key: str, n: int = 1) -> None:
    counts[key] = counts.get(key, 0) + n


def lift_variants(
    df: pl.DataFrame, *, chain: str, liftover_bin: str, fasta: Fasta, counts: dict,
) -> pl.DataFrame:
    """`df` has a seqname column `contig` (with or without `chr`), and `pos`, `ref`, `alt` in
    GRCh37, unique on the four. Returns the rows that lift, with `pos`, `ref`, `alt` in GRCh38;
    `contig` is unchanged because a variant that changes chromosome is dropped.
    """
    df = df.with_row_index("_row")
    is_mt = pl.col("contig").str.replace("(?i)^chr", "").str.to_uppercase().is_in(["M", "MT"])
    _count(counts, "liftover_mitochondrial", df.filter(is_mt).height)
    todo = df.filter(~is_mt)

    bed_rows = [
        f"chr{c.removeprefix('chr')}\t{pos - 1}\t{pos - 1 + len(ref)}\t{row}\t0\t+\n"
        for row, c, pos, ref in todo.select("_row", "contig", "pos", "ref").iter_rows()
    ]
    with tempfile.TemporaryDirectory() as tmp:
        lifted, n_unmapped = run_liftover(liftover_bin, chain, bed_rows, Path(tmp))
    _count(counts, "liftover_unmapped", n_unmapped)

    alleles = {row: (c, ref, alt) for row, c, ref, alt in todo.select("_row", "contig", "ref", "alt").iter_rows()}
    rows, new_pos, new_ref, new_alt = [], [], [], []
    for chrom, start, end, row, _, strand in lifted:
        row, start, end = int(row), int(start), int(end)
        contig, ref, alt = alleles[row]
        if chrom != f"chr{contig.removeprefix('chr')}":
            _count(counts, "liftover_different_chromosome")
            continue
        if end - start != len(ref):
            _count(counts, "liftover_span_length_changed")
            continue
        if strand == "-":
            ref, alt = reverse_complement(ref), reverse_complement(alt)
        target = fasta.fetch(chrom, start, end)
        if target != ref:
            swapped = len(ref) == 1 and len(alt) == 1 and target == alt
            _count(counts, "liftover_alt_is_grch38_ref" if swapped else "liftover_ref_mismatch")
            continue
        pos = start + 1
        if len(ref) != len(alt):
            pos, ref, alt = left_align(fasta, chrom, pos, ref, alt)
        if strand == "-":
            _count(counts, "liftover_minus_strand_kept")
        rows.append(row)
        new_pos.append(pos)
        new_ref.append(ref)
        new_alt.append(alt)

    lifted_df = pl.DataFrame(
        {"_row": rows, "pos": new_pos, "ref": new_ref, "alt": new_alt},
        schema={"_row": pl.UInt32, "pos": pl.Int64, "ref": pl.Utf8, "alt": pl.Utf8},
    )
    out = df.drop("pos", "ref", "alt").join(lifted_df, on="_row", how="inner").drop("_row")
    collides = pl.struct("contig", "pos", "ref", "alt").is_duplicated()
    _count(counts, "liftover_collision", out.filter(collides).height)
    return out.filter(~collides)


def lift_intervals(
    df: pl.DataFrame, *, chain: str, liftover_bin: str, counts: dict,
) -> pl.DataFrame:
    """`df` has `contig`, 1-based inclusive `start` and `end` in GRCh37. Returns the rows whose
    interval lifts in one piece (liftOver's default minMatch) onto the same chromosome, with
    `start` and `end` in GRCh38."""
    df = df.with_row_index("_row")
    bed_rows = [
        f"chr{c.removeprefix('chr')}\t{start - 1}\t{end}\t{row}\t0\t+\n"
        for row, c, start, end in df.select("_row", "contig", "start", "end").iter_rows()
    ]
    with tempfile.TemporaryDirectory() as tmp:
        lifted, n_unmapped = run_liftover(liftover_bin, chain, bed_rows, Path(tmp))
    _count(counts, "liftover_unmapped", n_unmapped)
    contigs = dict(df.select("_row", "contig").iter_rows())
    rows, starts, ends = [], [], []
    for chrom, start, end, row, _, _ in lifted:
        if chrom != f"chr{contigs[int(row)].removeprefix('chr')}":
            _count(counts, "liftover_different_chromosome")
            continue
        rows.append(int(row))
        starts.append(int(start) + 1)
        ends.append(int(end))
    lifted_df = pl.DataFrame(
        {"_row": rows, "start": starts, "end": ends},
        schema={"_row": pl.UInt32, "start": pl.Int64, "end": pl.Int64},
    )
    return df.drop("start", "end").join(lifted_df, on="_row", how="inner").drop("_row")
