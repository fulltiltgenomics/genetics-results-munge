#!/usr/bin/env python3
"""Munge a GWAS Catalog harmonised summary statistics file into the sumstat schema.

Input: the Catalog's harmonised GWAS-SSF file (`<GCST id>.h.tsv.gz`, under `harmonised/`
in the study's FTP directory) with its `-meta.yaml`. Columns read: `chromosome`,
`base_pair_location` (GRCh38 after harmonisation), `effect_allele`, `other_allele`,
`beta`, `standard_error`, `effect_allele_frequency`, `p_value`, `rsid`, `hm_code`.

Allele decision: for every row it could resolve, the harmoniser has already put the
reference-genome allele in `other_allele` and the alternative in `effect_allele`,
negating `beta` and taking 1 - EAF where the submitter had them the other way round
(hm_code 2, 4, 6, 8, 11, 13). So `ref` = `other_allele`, `alt` = `effect_allele`, and
`beta` and `af` pass through unchanged -- nothing is flipped here. Rows whose hm_code
reads "Not harmonised" (9 and 14-18) still carry the submitter's orientation, which this
script cannot verify, so they are dropped and counted per code. Indel representation is
the submitter's; `check_sumstat_vs_gnomad.py` measures it against gnomAD afterwards
instead of this script re-deriving it.

Chromosomes: `chr` prefix stripped case-insensitively, X -> 23; any other non-numeric
contig (Y, MT, scaffolds) is dropped and counted. mlog10p = -log10(p) rounded to 4
decimals, recomputed from beta/se where p underflows to 0.

The file is read line by line and written through bgzip, because a Catalog file can hold
well over a hundred million rows (GCST90444202: 146M) and this is run on a 15 GB machine.
The harmoniser sorts its output (`is_sorted: true` in the meta yaml); that order is
asserted while streaming rather than re-created, so an unsorted input fails instead of
indexing wrong.

Output: `<stem>.munged.tsv.gz` plus the `.mlog10p_gt4.tsv.gz` companion, both bgzipped and
tabix-indexed (`-s1 -b2 -e2`), columns `#chr pos ref alt rsid beta se mlog10p af`.
`--stage gs://bucket/dir/` uploads both with their indexes.

A harmonised file can be a subset of the submission (chromosome X and every variant the
harmoniser's reference VCF lacks are missing from GCST90444202's): read the `running.log`
next to it first. Where that matters, `munge_gwas_catalog_submitted.py` munges the submitted
file and uses this script's output as the independent check of its orientation.

Example (GCST90444202, Huerta-Chagoya et al. 2024 T2D):

    python3 scripts/munge_gwas_catalog.py \\
        --input /mnt/disks/data/gwascat/GCST90444202/GCST90444202.h.tsv.gz
"""

import argparse
import math
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

from scipy.special import log_ndtr

from sumstat_utils import upload_to_gcs

REQUIRED_COLUMNS = [
    "chromosome", "base_pair_location", "effect_allele", "other_allele",
    "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid", "hm_code",
]
OUTPUT_HEADER = "#chr\tpos\tref\talt\trsid\tbeta\tse\tmlog10p\taf\n"

# hm_code values whose description does not end in "Not harmonised"; the others are 9
# (palindromic, dropped by the harmoniser's own policy) and 14-18 (fields missing, no or
# several reference matches, palindromic strand not inferable)
HARMONISED_CODES = frozenset(str(c) for c in (1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13))
ALLELE_RE = re.compile(r"^[ACGT]+$")
BATCH_LINES = 20_000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="harmonised .h.tsv.gz file")
    parser.add_argument("--output", help="munged output path (default: <input stem>.munged.tsv.gz next to the input)")
    parser.add_argument("--stage", help="gs://bucket/dir/ to upload the munged files and their indexes to")
    parser.add_argument("--max-rows", type=int, help="stop after this many input rows (smoke tests)")
    return parser.parse_args()


def mlog10p_from_z(beta: float, se: float) -> float:
    # two-sided p from |z|, in log space so it does not underflow where p did
    return (-log_ndtr(-abs(beta / se)) - math.log(2)) / math.log(10)


class BgzipWriter:
    """Lines go to `bgzip -c` in batches; `close` indexes the file."""

    def __init__(self, path: str):
        self.path = path
        self.proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=open(path, "wb"))
        self.batch: list[str] = []
        self.n = 0

    def write(self, line: str) -> None:
        self.batch.append(line)
        self.n += 1
        if len(self.batch) >= BATCH_LINES:
            self.flush()

    def flush(self) -> None:
        if self.batch:
            self.proc.stdin.write("".join(self.batch).encode())
            self.batch = []

    def close(self) -> None:
        self.flush()
        self.proc.stdin.close()
        if self.proc.wait() != 0:
            raise SystemExit(f"bgzip failed for {self.path}")
        subprocess.run(["tabix", "-f", "-s1", "-b2", "-e2", self.path], check=True)


def chromosome_number(value: str) -> int | None:
    name = value[3:] if value[:3].lower() == "chr" else value
    if name == "X":
        return 23
    return int(name) if name.isdigit() else None


def munge(input_path: str, output_path: str, max_rows: int | None) -> dict:
    filtered_path = output_path.replace(".tsv.gz", ".mlog10p_gt4.tsv.gz")
    zcat = subprocess.Popen(["zcat", input_path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    header = zcat.stdout.readline().rstrip("\n").split("\t")
    missing = [c for c in REQUIRED_COLUMNS if c not in header]
    if missing:
        raise SystemExit(f"{input_path} lacks required columns: {missing}")
    col = {c: header.index(c) for c in REQUIRED_COLUMNS}
    i_chr, i_pos, i_eff, i_oth = col["chromosome"], col["base_pair_location"], col["effect_allele"], col["other_allele"]
    i_beta, i_se, i_af, i_p, i_rs, i_hm = (
        col["beta"], col["standard_error"], col["effect_allele_frequency"], col["p_value"], col["rsid"], col["hm_code"],
    )
    width = len(header)

    full = BgzipWriter(output_path)
    filt = BgzipWriter(filtered_path)
    full.write(OUTPUT_HEADER)
    filt.write(OUTPUT_HEADER)

    dropped = Counter()
    seen_chromosomes: set[int] = set()
    current_chr, current_pos = None, -1
    n_in = n_p_zero = 0
    for line in zcat.stdout:
        n_in += 1
        if max_rows and n_in > max_rows:
            zcat.kill()
            break
        if n_in % 10_000_000 == 0:
            print(f"  {n_in / 1e6:.0f}M rows read, {full.n - 1} written", flush=True)
        f = line.rstrip("\n").split("\t")
        if len(f) != width:
            dropped["malformed row"] += 1
            continue
        if f[i_hm] not in HARMONISED_CODES:
            dropped[f"hm_code {f[i_hm]} (not harmonised)"] += 1
            continue
        chrom = chromosome_number(f[i_chr])
        if chrom is None:
            dropped[f"contig {f[i_chr]}"] += 1
            continue
        ref, alt = f[i_oth], f[i_eff]
        if not (ALLELE_RE.match(ref) and ALLELE_RE.match(alt)) or ref == alt:
            dropped["malformed alleles"] += 1
            continue
        try:
            pos = int(f[i_pos])
            beta = float(f[i_beta])
            se = float(f[i_se])
        except ValueError:
            dropped["unparseable pos/beta/se"] += 1
            continue
        if not (se > 0) or math.isnan(beta):
            dropped["se <= 0 or beta NA"] += 1
            continue

        if chrom != current_chr:
            if chrom in seen_chromosomes:
                raise SystemExit(f"input is not sorted: chromosome {chrom} reappears at row {n_in}")
            seen_chromosomes.add(chrom)
            current_chr, current_pos = chrom, -1
        elif pos < current_pos:
            raise SystemExit(f"input is not sorted: {chrom}:{pos} after {chrom}:{current_pos} at row {n_in}")
        current_pos = pos

        try:
            p = float(f[i_p])
        except ValueError:
            p = 0.0
        if p > 0:
            mlog10p = round(-math.log10(p), 4) + 0.0  # round(-log10(1.0)) is -0.0
        else:
            n_p_zero += 1
            mlog10p = round(mlog10p_from_z(beta, se), 4)

        try:
            af = f"{round(float(f[i_af]), 6)}"
        except ValueError:
            af = "NA"
        rsid = f[i_rs] if f[i_rs] else "NA"

        beta += 0.0  # the harmoniser negates a zero beta into -2.2e-308; -0.0 prints with a sign
        out = f"{chrom}\t{pos}\t{ref}\t{alt}\t{rsid}\t{beta:.3e}\t{se:.3e}\t{mlog10p}\t{af}\n"
        full.write(out)
        if mlog10p > 4:
            filt.write(out)

    if not max_rows and zcat.wait() != 0:
        raise SystemExit(f"zcat failed on {input_path}")
    full.close()
    filt.close()
    return {
        "rows_in": n_in, "rows_out": full.n - 1, "rows_mlog10p_gt4": filt.n - 1,
        "p_underflow": n_p_zero, "dropped": dropped,
        "filtered_path": filtered_path,
    }


def main() -> None:
    args = parse_args()
    output = args.output or str(Path(args.input).with_name(Path(args.input).name.replace(".tsv.gz", "") + ".munged.tsv.gz"))

    print(f"Munging {args.input} -> {output}")
    result = munge(args.input, output, args.max_rows)
    print(f"  rows in: {result['rows_in']}, written: {result['rows_out']}, mlog10p > 4: {result['rows_mlog10p_gt4']}")
    print(f"  p underflow (mlog10p from beta/se): {result['p_underflow']}")
    for reason, n in sorted(result["dropped"].items(), key=lambda x: -x[1]):
        print(f"  dropped {n}: {reason}")

    if args.stage:
        prefix = args.stage.rstrip("/") + "/"
        for local in (output, result["filtered_path"]):
            upload_to_gcs(local, prefix + Path(local).name)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
