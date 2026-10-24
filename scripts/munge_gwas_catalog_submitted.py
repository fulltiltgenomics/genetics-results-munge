#!/usr/bin/env python3
"""Munge a GWAS Catalog study from the SUBMITTED file, orienting every row against gnomAD and the reference.

Why not the harmonised file: the Catalog's harmoniser drops every variant it cannot match
to its reference VCF (hm_code 15/16) and, for GCST90444202, all of chromosome X -- 9.8M + 6M
of 146.2M submitted rows, and the dropped ones skew rare, which is what that study is for.
The submitted file has them all but no rsids, its effect allele is the reference allele on
about half the rows, and some alleles carry shared padding bases (`TA>AA` for a `T>A`). So
this script re-derives the orientation itself, for all rows, from the reference fasta first
and gnomAD where the fasta cannot tell, and `munge_gwas_catalog.py`'s output from the
harmonised file becomes a validation set rather than the product (`compare`).

Steps, each a subcommand, each streaming or an external sort, so 146M rows fit a 15 GB machine:

  normalize  read the submitted GWAS-SSF file; trim the shared suffix then the shared prefix of
             each allele pair (one base always kept, the position advancing with the prefix),
             map X to 23, drop other contigs, write and externally sort by chr, pos
  lookup     one BigQuery query: the normalized rows (a GCS external table) joined to the gnomAD
             table on chr, pos and the allele pair in either orientation, exported back to GCS
             and sorted locally
  orient     decide ref/alt per row: the allele that equals the fasta sequence at pos is ref,
             which settles every SNP and every indel whose alleles differ in their first base.
             When both alleles match (the shorter is a prefix of the longer: an insertion and
             a deletion sharing a position, e.g. A>AT and AT>A, two different variants) the
             orientation whose gnomAD AF is closer to the submitted EAF on the logit scale
             wins, an orientation gnomAD does not carry counting as AF 1e-5; a tie, or neither
             allele matching the fasta, drops the row. beta is negated and af becomes 1 - EAF
             where the effect allele turned out to be ref; rsid comes from the gnomAD row in
             the chosen orientation, NA when there is none. Duplicate (chr, pos, ref, alt)
             after orientation keeps the first. Writes the munged file and its mlog10p > 4
             companion, bgzipped and tabix-indexed (`--stage` uploads them), plus a
             `.decisions.tsv.gz` with one line per row saying how it was decided or why dropped
  compare    lockstep the munged file against the harmonised-file munge: for every variant in
             both, count agreement of ref/alt, beta sign, af and rsid, plus the rows only one
             side has -- the orientation check against an independent harmonisation

Run from `scripts/` so the sibling modules import. Example (GCST90444202):

    python3 munge_gwas_catalog_submitted.py normalize --submitted portal.tsv.bgz --out norm.tsv.gz
    python3 munge_gwas_catalog_submitted.py lookup --rows norm.tsv.gz --work gs://bucket/dir/work/ --out gnomad.tsv.gz
    python3 munge_gwas_catalog_submitted.py orient --rows norm.tsv.gz --gnomad gnomad.tsv.gz \\
        --fasta /mnt/disks/data/ref/hg38.fa --out GCST90444202.munged.tsv.gz --stage gs://bucket/dir/
    python3 munge_gwas_catalog_submitted.py compare --munged GCST90444202.munged.tsv.gz \\
        --harmonised GCST90444202.h.munged.tsv.gz
"""

import argparse
import glob
import gzip
import math
import os
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

from munge_gwas_catalog import OUTPUT_HEADER, BgzipWriter, chromosome_number, mlog10p_from_z
from sumstat_utils import upload_to_gcs

GNOMAD_TABLE = "daly-finngenie.genetics_results.gnomad_variant_annotation"
ROW_HEADER = "chr\tpos\ta1\ta2\tbeta\tse\teaf\tp\n"
SORT = "LC_ALL=C sort -u -t $'\\t' -k1,1n -k2,2n -S 3G --parallel=4"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="step", required=True)

    p = sub.add_parser("normalize")
    p.add_argument("--submitted", required=True, help="the submitted (portal) file, gzipped, GWAS-SSF columns")
    p.add_argument("--out", required=True, help="normalized, sorted rows .tsv.gz")
    p.add_argument("--tmp", default="/mnt/disks/data/tmp", help="directory for sort's temporary files")
    p.add_argument("--max-rows", type=int, help="stop after this many submitted rows (smoke tests)")

    p = sub.add_parser("lookup")
    p.add_argument("--rows", required=True, help="output of normalize")
    p.add_argument("--work", required=True, help="gs://bucket/prefix/ for the external table and the export")
    p.add_argument("--out", required=True, help="gnomAD rows for the input, .tsv.gz, sorted by chr and pos")
    p.add_argument("--tmp", default="/mnt/disks/data/tmp")
    p.add_argument("--gnomad-table", default=GNOMAD_TABLE)

    p = sub.add_parser("orient")
    p.add_argument("--rows", required=True, help="output of normalize")
    p.add_argument("--gnomad", required=True, help="output of lookup")
    p.add_argument("--fasta", required=True, help="GRCh38 reference fasta (UCSC hg38.fa, chr-prefixed)")
    p.add_argument("--out", required=True, help="munged .tsv.gz; the mlog10p_gt4 companion is written next to it")
    p.add_argument("--stage", help="gs://bucket/dir/ to upload both files and their indexes to")

    p = sub.add_parser("compare")
    p.add_argument("--munged", required=True)
    p.add_argument("--harmonised", required=True, help="munge_gwas_catalog.py output from the harmonised file")
    return parser.parse_args()


def open_text(path: str):
    return subprocess.Popen(["zcat", path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20).stdout


def rows_by_position(lines, chr_col: int = 0, pos_col: int = 1):
    """Group consecutive rows sharing (chr, pos); yields ((chr, pos), rows) and asserts order."""
    group, key, last = [], None, (0, -1)
    for line in lines:
        f = line.rstrip("\n").split("\t")
        k = (int(f[chr_col]), int(f[pos_col]))
        if k != key:
            if group:
                yield key, group
            if k < last:
                raise SystemExit(f"not sorted: {k} after {last}")
            group, key, last = [], k, k
        group.append(f)
    if group:
        yield key, group


def trim_alleles(pos: int, a: str, b: str) -> tuple[int, str, str]:
    """Drop the bases two alleles share: suffix first, then prefix (which moves pos), one base always kept."""
    while len(a) > 1 and len(b) > 1 and a[-1] == b[-1]:
        a, b = a[:-1], b[:-1]
    while len(a) > 1 and len(b) > 1 and a[0] == b[0]:
        a, b, pos = a[1:], b[1:], pos + 1
    return pos, a, b


# --- normalize ---------------------------------------------------------------------------

def normalize(args: argparse.Namespace) -> None:
    sub = open_text(args.submitted)
    header = sub.readline().rstrip("\n").split("\t")
    col = {c: header.index(c) for c in (
        "chromosome", "base_pair_location", "effect_allele", "other_allele",
        "beta", "standard_error", "effect_allele_frequency", "p_value",
    )}
    i_chr, i_pos, i_a1, i_a2 = col["chromosome"], col["base_pair_location"], col["effect_allele"], col["other_allele"]
    i_beta, i_se, i_eaf, i_p = col["beta"], col["standard_error"], col["effect_allele_frequency"], col["p_value"]

    os.makedirs(args.tmp, exist_ok=True)
    sort = subprocess.Popen(
        f"{SORT} -T {args.tmp} | gzip > {args.out}", shell=True, executable="/bin/bash",
        stdin=subprocess.PIPE, text=True, bufsize=1 << 20,
    )
    counts = Counter()
    n = 0
    for line in sub:
        n += 1
        if args.max_rows and n > args.max_rows:
            break
        if n % 10_000_000 == 0:
            print(f"  {n / 1e6:.0f}M submitted rows", flush=True)
        f = line.rstrip("\n").split("\t")
        chrom = chromosome_number(f[i_chr])
        if chrom is None:
            counts[f"dropped: contig {f[i_chr]}"] += 1
            continue
        a1, a2 = f[i_a1].upper(), f[i_a2].upper()
        if a1 == a2 or not a1 or not a2:
            counts["dropped: identical or empty alleles"] += 1
            continue
        pos, t1, t2 = trim_alleles(int(f[i_pos]), a1, a2)
        if (t1, t2) != (a1, a2):
            counts["alleles trimmed"] += 1
        counts["kept"] += 1
        sort.stdin.write(f"{chrom}\t{pos}\t{t1}\t{t2}\t{f[i_beta]}\t{f[i_se]}\t{f[i_eaf]}\t{f[i_p]}\n")
    sort.stdin.close()
    if sort.wait() != 0:
        raise SystemExit("sort failed")
    print(f"  submitted rows: {n}")
    for key, v in counts.most_common():
        print(f"  {key}: {v}")


# --- lookup ------------------------------------------------------------------------------

def lookup(args: argparse.Namespace) -> None:
    work = args.work.rstrip("/") + "/"
    # the inline external-table definition takes comma-separated input with no header
    csv_local = tempfile.mktemp(suffix=".csv.gz", dir=args.tmp)
    subprocess.run(
        f"zcat {args.rows} | cut -f1-4 | tr '\\t' ',' | gzip > {csv_local}",
        shell=True, check=True, executable="/bin/bash",
    )
    csv_gcs = work + "rows.csv.gz"
    subprocess.run(["gcloud", "storage", "cp", csv_local, csv_gcs], check=True)
    os.unlink(csv_local)

    export = work + "lookup-*.tsv.gz"
    # overwrite=true replaces shards of the same name; a previous export with more shards
    # leaves the rest behind and every row of them is read twice
    subprocess.run(["gcloud", "storage", "rm", export], capture_output=True)
    query = f"""
    EXPORT DATA OPTIONS (uri='{export}', format='CSV', field_delimiter='\\t', compression='GZIP', overwrite=true) AS
    SELECT c.chr, c.pos, c.a1, c.a2, g.ref, g.alt, g.rsids, g.AF, g.AF_nfe, g.genome_or_exome
    FROM submitted c
    JOIN `{args.gnomad_table}` g
      ON g.chr = c.chr AND g.pos = c.pos
     AND ((g.ref = c.a2 AND g.alt = c.a1) OR (g.ref = c.a1 AND g.alt = c.a2))
    """
    subprocess.run([
        "bq", "query", "--use_legacy_sql=false", "--nouse_cache",
        f"--external_table_definition=submitted::chr:INTEGER,pos:INTEGER,a1:STRING,a2:STRING@CSV={csv_gcs}",
        query,
    ], check=True)

    tmpdir = tempfile.mkdtemp(dir=args.tmp)
    subprocess.run(["gcloud", "storage", "cp", export, tmpdir + "/"], check=True)
    # the export is thousands of small files (too many for one argument list) in no
    # particular order; concatenated gzip members decompress as one stream, and the sort
    # is what lets orient merge by position
    subprocess.run(
        f"find {tmpdir} -name 'lookup-*.tsv.gz' -print0 | xargs -0 cat | zcat | {SORT} -T {args.tmp} | gzip > {args.out}",
        shell=True, check=True, executable="/bin/bash",
    )
    for p in glob.glob(tmpdir + "/lookup-*.tsv.gz"):
        os.unlink(p)
    os.rmdir(tmpdir)
    subprocess.run(["gcloud", "storage", "rm", csv_gcs, export], check=True)
    print(f"  gnomAD rows written to {args.out}")


# --- orient ------------------------------------------------------------------------------

def orient(args: argparse.Namespace) -> None:
    from pyfaidx import Fasta

    fasta = Fasta(args.fasta, as_raw=True, sequence_always_upper=True)
    sequence = {"name": None, "seq": ""}

    def reference_matches(chrom: int, pos: int, allele: str) -> bool:
        # the whole chromosome is held as one string: 146M rows of pyfaidx seeks would take hours
        name = "chrX" if chrom == 23 else f"chr{chrom}"
        if sequence["name"] != name:
            sequence["name"], sequence["seq"] = name, fasta[name][:]
        return sequence["seq"][pos - 1: pos - 1 + len(allele)] == allele

    def logit_af(value: float | None) -> float:
        # frequencies are compared on the logit scale, so 0.18 against an absent variant is
        # further apart than 0.82 against 0.59; absent, 0 and 1 are clamped to the floor
        value = 1e-5 if value is None else min(max(value, 1e-5), 1 - 1e-5)
        return math.log(value / (1 - value))

    def gnomad_af(rows: list) -> float | None:
        try:
            return float(rows[0][7]) if rows else None
        except ValueError:
            return None

    rows = open_text(args.rows)
    gnomad = rows_by_position(open_text(args.gnomad))
    g_key, g_rows = next(gnomad, ((99, 0), []))

    filtered_path = args.out.replace(".tsv.gz", ".mlog10p_gt4.tsv.gz")
    full = BgzipWriter(args.out)
    filt = BgzipWriter(filtered_path)
    full.write(OUTPUT_HEADER)
    filt.write(OUTPUT_HEADER)
    decisions = gzip.open(args.out.replace(".tsv.gz", "") + ".decisions.tsv.gz", "wt")
    decisions.write("chr\tpos\ta1\ta2\tdecision\n")
    counts = Counter()
    n = 0
    for key, group in rows_by_position(rows):
        n += len(group)
        if n // 10_000_000 != (n - len(group)) // 10_000_000:
            print(f"  {n / 1e6:.0f}M rows, {full.n - 1} written", flush=True)
        while g_key < key:
            g_key, g_rows = next(gnomad, ((99, 0), []))
        at_pos = g_rows if g_key == key else []
        seen = set()
        for f in group:
            chrom, pos, a1, a2 = key[0], key[1], f[2], f[3]
            try:
                beta, se, eaf = float(f[4]), float(f[5]), float(f[6])
            except ValueError:
                counts["dropped: beta/se/eaf unparseable"] += 1
                decisions.write(f"{chrom}\t{pos}\t{a1}\t{a2}\tdropped_unparseable\n")
                continue
            if not se > 0 or math.isnan(beta):
                counts["dropped: se <= 0 or beta NA"] += 1
                decisions.write(f"{chrom}\t{pos}\t{a1}\t{a2}\tdropped_se\n")
                continue

            # gnomAD rows whose ref/alt is this pair in either orientation; genomes preferred
            as_alt = sorted((g for g in at_pos if g[4] == a2 and g[5] == a1), key=lambda g: g[9] != "g")
            as_ref = sorted((g for g in at_pos if g[4] == a1 and g[5] == a2), key=lambda g: g[9] != "g")
            m1, m2 = reference_matches(chrom, pos, a1), reference_matches(chrom, pos, a2)
            if m1 != m2:
                # one allele is the reference sequence at pos: a SNP, or an indel whose
                # alleles differ in their first base. Nothing else can decide it differently
                effect_is_alt = m2
                decision = "reference"
            elif not m1:
                counts["dropped: neither allele matches the reference"] += 1
                decisions.write(f"{chrom}\t{pos}\t{a1}\t{a2}\tdropped_reference_neither\n")
                continue
            else:
                # both match: the shorter allele is a prefix of the longer and the genome reads
                # the longer, so the pair is a deletion (ref = longer) or an insertion (ref =
                # shorter), two different variants. The one whose alt frequency agrees with
                # gnomAD wins; an orientation gnomAD does not carry counts as the floor, so
                # where gnomAD has neither the minor allele becomes alt
                d_alt = abs(logit_af(eaf) - logit_af(gnomad_af(as_alt)))
                d_ref = abs(logit_af(1.0 - eaf) - logit_af(gnomad_af(as_ref)))
                if d_alt == d_ref:
                    counts["dropped: ambiguous indel, AF undecidable"] += 1
                    decisions.write(f"{chrom}\t{pos}\t{a1}\t{a2}\tdropped_ambiguous_indel\n")
                    continue
                effect_is_alt = d_alt < d_ref
                decision = "ambiguous_indel_by_af"
            chosen = (as_alt if effect_is_alt else as_ref) or [None]
            chosen = chosen[0]
            rsid = chosen[6].split(",")[0] if chosen is not None and chosen[6] not in ("", "NA") else "NA"
            ref, alt = (a2, a1) if effect_is_alt else (a1, a2)
            if (ref, alt) in seen:
                counts["dropped: duplicate variant"] += 1
                decisions.write(f"{chrom}\t{pos}\t{a1}\t{a2}\tdropped_duplicate\n")
                continue
            seen.add((ref, alt))
            orientation = "alt" if effect_is_alt else "ref"
            counts[f"{decision}, effect is {orientation}"] += 1
            decisions.write(f"{chrom}\t{pos}\t{a1}\t{a2}\t{decision}_{orientation}\n")

            if not effect_is_alt:
                beta, eaf = -beta, 1.0 - eaf
            beta += 0.0  # a negated zero would print as -0.000e+00
            try:
                p = float(f[7])
            except ValueError:
                p = 0.0
            if p > 0:
                mlog10p = round(-math.log10(p), 4) + 0.0
            else:
                counts["p underflow (mlog10p from beta/se)"] += 1
                mlog10p = round(mlog10p_from_z(beta, se), 4)
            line = f"{chrom}\t{pos}\t{ref}\t{alt}\t{rsid}\t{beta:.3e}\t{se:.3e}\t{mlog10p}\t{round(eaf, 6)}\n"
            full.write(line)
            if mlog10p > 4:
                filt.write(line)
    full.close()
    filt.close()
    decisions.close()
    print(f"  rows in: {n}, written: {full.n - 1}, mlog10p > 4: {filt.n - 1}")
    for key, v in counts.most_common():
        print(f"  {key}: {v}")
    if args.stage:
        prefix = args.stage.rstrip("/") + "/"
        for local in (args.out, filtered_path):
            upload_to_gcs(local, prefix + Path(local).name)


# --- compare -----------------------------------------------------------------------------

def compare(args: argparse.Namespace) -> None:
    munged = rows_by_position(_skip_header(open_text(args.munged)))
    harmonised = rows_by_position(_skip_header(open_text(args.harmonised)))
    counts = Counter()
    m_key, m_rows = next(munged, ((99, 0), []))
    for h_key, h_rows in harmonised:
        while m_key < h_key:
            counts["munged only"] += len(m_rows)
            m_key, m_rows = next(munged, ((99, 0), []))
        # a shared position is consumed here, or the next advance would count it as munged only
        at_pos = {(r[2], r[3]): r for r in m_rows} if m_key == h_key else {}
        if m_key == h_key:
            m_key, m_rows = next(munged, ((99, 0), []))
        matched = set()
        for h in h_rows:
            m = at_pos.get((h[2], h[3]))
            if m is not None:
                matched.add((h[2], h[3]))
                counts["same ref/alt"] += 1
                bm, bh = float(m[5]), float(h[5])
                # the harmoniser negates a zero beta into -2.2e-308; a sign is only a sign away from zero
                counts["beta sign agrees" if (bm >= 0) == (bh >= 0) or abs(bh) < 1e-300 or abs(bm) < 1e-300 else "BETA SIGN DIFFERS"] += 1
                counts["af agrees" if abs(float(m[8]) - float(h[8])) < 1e-5 else "AF DIFFERS"] += 1
                counts["rsid agrees" if m[4] == h[4] else "rsid differs"] += 1
            elif (h[3], h[2]) in at_pos:
                matched.add((h[3], h[2]))
                counts["REF/ALT SWAPPED"] += 1
            else:
                counts["harmonised only"] += 1
        counts["munged only"] += len(set(at_pos) - matched)
    while m_key != (99, 0):
        counts["munged only"] += len(m_rows)
        m_key, m_rows = next(munged, ((99, 0), []))
    for key, v in sorted(counts.items()):
        print(f"  {key}: {v}")


def _skip_header(lines):
    next(lines)
    return lines


def main() -> None:
    args = parse_args()
    {"normalize": normalize, "lookup": lookup, "orient": orient, "compare": compare}[args.step](args)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
