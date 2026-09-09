#!/usr/bin/env python3
"""One-off measurement: lift the Collins et al. 2022 rCNV sliding windows GRCh37 -> GRCh38.

Answers the descope trigger on the rCNV epic ("more than ~2% of windows fail to lift cleanly
-> the sliding-window product is dropped"). It creates no table and writes nothing outside
--work-dir; the table-building munge is a separate script.

INPUT
  --sumstats-dir: the unpacked Collins_rCNV_2022.sliding_window_sumstats directory (108
  bgzipped BED files, 54 phenotypes x DEL/DUP). Every file carries the SAME 267,237 windows
  (200 kb, 10 kb step, GRCh37 autosomes) in the same order; --verify-files N re-reads N files
  and asserts the (chr,start,end) triples are byte-identical before using the first one.

PROCEDURE (follows scripts/munge_calderon.py, the repo's only liftOver munge)
  1. write the unique windows as BED4 with name = "<chr>-<start>-<end>" (GRCh37 key), seqnames
     chr-prefixed because the hg19ToHg38 chain is keyed on them;
  2. `liftOver in.bed chain out.bed unmapped.bed` with UCSC defaults (minMatch 0.95, no
     -multiple) — the same invocation munge_calderon.py makes;
  3. drop keys that map more than once (multi-mapped / split);
  4. drop keys that land on a different chromosome than their GRCh37 chromosome.
  munge_calderon.py stops there because its peaks have no expected width. A 200 kb window does,
  so this measurement adds the task's fifth filter:
  5. drop keys whose lifted length falls outside [--min-len, --max-len] (default 180-220 kb).

  liftOver reports a whole interval as unmapped without saying which end failed, so a SECOND
  liftOver pass lifts the two 1 bp endpoints of every window independently. That pass is used
  only to attribute an interval-level failure to "unmapped start" / "unmapped end" / "both
  ends unmapped"; it never rescues a window. A window whose two endpoints both lift but whose
  interval liftOver still refused is classed "interior_deleted" — liftOver's "Partially deleted
  in new" at the default minMatch 0.95, i.e. more than 5% of the interval has no hg38 image.

OUTPUT
  stdout: the counts table, the lifted-length distribution, the neighbour-geometry counts and
  the DECISION line. docs/rcnv-sliding-windows.md pastes that output verbatim.
  --work-dir also keeps dropped_windows.tsv (every dropped window with its class) and the raw
  liftOver in/out/unmapped BEDs.

Re-run:
  python3 scripts/rcnv_liftover_windows.py --sumstats-dir <dir> \
      --liftover-bin <path/liftOver> --chain <path/hg19ToHg38.over.chain.gz> --work-dir <dir>
Chain and binary source (UCSC, not vendored here):
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
  https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
"""

import argparse
import gzip
import hashlib
import shutil
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
BIN_URL = "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver"

# ordered so the report reads as a funnel: every window is in exactly one bucket
FAILURE_ORDER = [
    "unmapped_start", "unmapped_end", "unmapped_both", "split",
    "chromosome_changed", "length_out_of_tolerance", "interior_deleted",
]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sumstats-dir", required=True,
                   help="unpacked Collins_rCNV_2022.sliding_window_sumstats directory")
    p.add_argument("--liftover-bin", default="liftOver", help="UCSC liftOver binary")
    p.add_argument("--chain", required=True, help=f"hg19ToHg38 chain(.gz); source: {CHAIN_URL}")
    p.add_argument("--work-dir", required=True, help="scratch dir for BEDs and dropped_windows.tsv")
    p.add_argument("--min-len", type=int, default=180_000, help="lifted-length floor (default 180000)")
    p.add_argument("--max-len", type=int, default=220_000, help="lifted-length ceiling (default 220000)")
    p.add_argument("--step", type=int, default=10_000, help="GRCh37 window step (default 10000)")
    p.add_argument("--verify-files", type=int, default=2,
                   help="how many sumstats files to read and compare window-for-window (default 2)")
    return p.parse_args()


def read_windows(path: Path):
    """(chrom, start, end) triples of one sumstats BED, header dropped."""
    out = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c, s, e = line.split("\t", 3)[:3]
            out.append((c, int(s), int(e)))
    return out


def load_window_set(sumstats_dir: Path, verify_files: int):
    files = sorted(p for p in sumstats_dir.iterdir() if p.name.endswith(".bed.gz"))
    if not files:
        raise SystemExit(f"no *.bed.gz under {sumstats_dir}")
    checked = files[:max(1, verify_files)]
    windows, digests = None, []
    for p in checked:
        w = read_windows(p)
        digests.append(hashlib.sha256(
            "\n".join(f"{c}\t{s}\t{e}" for c, s, e in w).encode()).hexdigest())
        if windows is None:
            windows = w
    print(f"sumstats files:            {len(files)}")
    def display_label(p):
        parts = p.name.split(".")
        return ".".join(parts[:1] + parts[2:3]) if len(parts) >= 3 else p.name

    print(f"window set verified on:    {len(checked)} files "
          f"({', '.join(display_label(p) for p in checked)})")
    for p, d in zip(checked, digests):
        print(f"  sha256(chr,start,end)    {d[:16]}  {p.name}")
    if len(set(digests)) != 1:
        raise SystemExit("window sets differ between files — the measurement's premise is false")
    print(f"window set identical:      yes ({len(windows)} windows per file)")
    return windows


def run_liftover(binary, chain, bed_in, bed_out, unmapped):
    try:
        subprocess.run([binary, str(bed_in), str(chain), str(bed_out), str(unmapped)],
                       check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.stderr.decode() if isinstance(e.stderr, bytes) else str(e.stderr))
        raise


def write_bed(path, rows):
    with open(path, "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def read_mapped(path):
    """name -> list of (chrom, start, end); a list longer than 1 is a multi-map."""
    out = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            c, s, e, name = line.rstrip("\n").split("\t")[:4]
            out[name].append((c, int(s), int(e)))
    return out


def read_unmapped(path):
    """name -> liftOver's own reason string ('Deleted in new', 'Split in new', ...)."""
    out, reason = {}, "unknown"
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                reason = line.lstrip("#").strip()
                continue
            out[line.rstrip("\n").split("\t")[3]] = reason
    return out


def main():
    args = parse_args()
    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    if shutil.which(args.liftover_bin) is None and not Path(args.liftover_bin).exists():
        raise SystemExit(f"liftOver binary not found: {args.liftover_bin} (source: {BIN_URL})")
    if not Path(args.chain).exists():
        raise SystemExit(f"chain not found: {args.chain} (source: {CHAIN_URL})")

    md5 = hashlib.md5()
    with open(args.chain, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            md5.update(chunk)
    chain_md5 = md5.hexdigest()
    print("=" * 78)
    print("rCNV sliding windows — UCSC liftOver hg19 -> hg38")
    print("=" * 78)
    print(f"liftOver binary:           {args.liftover_bin}")
    print(f"chain:                     {Path(args.chain).name}  md5 {chain_md5}")
    print(f"chain source:              {CHAIN_URL}")
    print(f"tolerance:                 lifted length in [{args.min_len:,}, {args.max_len:,}] bp, "
          f"same chromosome")
    print()

    windows = load_window_set(Path(args.sumstats_dir), args.verify_files)
    total = len(windows)
    keys = [f"{c}-{s}-{e}" for c, s, e in windows]
    src = {k: w for k, w in zip(keys, windows)}
    if len(src) != total:
        raise SystemExit(f"{total - len(src)} duplicate (chr,start,end) triples in the window set")
    print()

    # pass 1 — the interval, exactly as munge_calderon.py lifts its peaks
    iv_in, iv_out, iv_un = work / "windows.hg19.bed", work / "windows.hg38.bed", work / "windows.unmapped.bed"
    write_bed(iv_in, [(f"chr{c}", s, e, k) for (c, s, e), k in zip(windows, keys)])
    run_liftover(args.liftover_bin, args.chain, iv_in, iv_out, iv_un)
    mapped = read_mapped(iv_out)
    unmapped_reason = read_unmapped(iv_un)

    # pass 2 — the two endpoints alone, ONLY to attribute an interval failure to an end
    ep_in, ep_out, ep_un = work / "ends.hg19.bed", work / "ends.hg38.bed", work / "ends.unmapped.bed"
    ends = []
    for (c, s, e), k in zip(windows, keys):
        ends.append((f"chr{c}", s, s + 1, k + "|s"))
        ends.append((f"chr{c}", e - 1, e, k + "|e"))
    write_bed(ep_in, ends)
    run_liftover(args.liftover_bin, args.chain, ep_in, ep_out, ep_un)
    ends_mapped = set(read_mapped(ep_out))

    clean, failures, lifted = {}, {}, []
    for k in keys:
        c37, s37, e37 = src[k]
        hits = mapped.get(k, [])
        if not hits:
            ms, me = (k + "|s") in ends_mapped, (k + "|e") in ends_mapped
            reason = unmapped_reason.get(k, "unknown")
            if "Split" in reason:
                failures[k] = "split"
            elif not ms and not me:
                failures[k] = "unmapped_both"
            elif not ms:
                failures[k] = "unmapped_start"
            elif not me:
                failures[k] = "unmapped_end"
            else:
                failures[k] = "interior_deleted"
            continue
        if len(hits) > 1:
            failures[k] = "split"
            continue
        c38, s38, e38 = hits[0]
        if c38 != f"chr{c37}":
            failures[k] = "chromosome_changed"
            continue
        if not (args.min_len <= e38 - s38 <= args.max_len):
            failures[k] = "length_out_of_tolerance"
            continue
        clean[k] = (c38, s38, e38)
        lifted.append(e38 - s38)

    n_clean, n_fail = len(clean), len(failures)
    fail_pct = 100.0 * n_fail / total
    counts = Counter(failures.values())

    print("liftOver reasons reported for intervals not mapped (raw, pass 1):")
    for reason, n in sorted(Counter(unmapped_reason.values()).items(), key=lambda kv: -kv[1]):
        print(f"  {reason:<34} {n:>8,}")
    if not unmapped_reason:
        print("  (none)")
    print()

    print("WINDOW COUNTS")
    print(f"  {'total windows':<34} {total:>8,}  100.000%")
    print(f"  {'mapped cleanly':<34} {n_clean:>8,}  {100.0 * n_clean / total:7.3f}%")
    for cls in FAILURE_ORDER:
        n = counts.get(cls, 0)
        print(f"  {'  failed: ' + cls:<34} {n:>8,}  {100.0 * n / total:7.3f}%")
    print(f"  {'failed (all classes)':<34} {n_fail:>8,}  {fail_pct:7.3f}%")
    print()

    print("LIFTED LENGTH (clean windows, bp)")
    if lifted:
        srt = sorted(lifted)
        print(f"  {'min':<34} {srt[0]:>8,}")
        print(f"  {'median':<34} {int(statistics.median(srt)):>8,}")
        print(f"  {'max':<34} {srt[-1]:>8,}")
        print(f"  {'exactly 200,000':<34} {sum(1 for x in srt if x == 200_000):>8,}  "
              f"{100.0 * sum(1 for x in srt if x == 200_000) / len(srt):7.3f}%")
    else:
        print("  (no clean windows)")
    print()

    # neighbour geometry: consecutive GRCh37 windows on a chromosome are `step` apart and
    # overlap by (200kb - step). hpa1.12 stores start_grch37/end_grch37 beside the lifted
    # coords, so a pair whose lifted spacing no longer matches is what a reader would trip on.
    by_chrom = defaultdict(list)
    for k, (c38, s38, e38) in clean.items():
        c37, s37, e37 = src[k]
        by_chrom[int(c37)].append((s37, e37, s38, e38, k))
    reordered, spacing_off, involved = 0, 0, set()
    pairs = 0
    for c in by_chrom:
        rows = sorted(by_chrom[c])
        for (s37a, e37a, s38a, e38a, ka), (s37b, e37b, s38b, e38b, kb) in zip(rows, rows[1:]):
            pairs += 1
            if s38b <= s38a:
                reordered += 1
                involved.update((ka, kb))
                continue
            ov37, ov38 = e37a - s37b, e38a - s38b
            if abs(ov38 - ov37) > args.step:
                spacing_off += 1
                involved.update((ka, kb))
    print(f"NEIGHBOUR GEOMETRY (adjacent clean windows on the same chromosome, {pairs:,} pairs)")
    print(f"  {'pairs reordered after lifting':<34} {reordered:>8,}  {100.0 * reordered / pairs:7.3f}%")
    print(f"  {'pairs whose overlap moved > step':<34} {spacing_off:>8,}  {100.0 * spacing_off / pairs:7.3f}%")
    print(f"  {'clean windows in such a pair':<34} {len(involved):>8,}  "
          f"{100.0 * len(involved) / n_clean:7.3f}%")
    print()

    # dropped windows: full list to disk, merged runs to stdout
    dropped_path = work / "dropped_windows.tsv"
    with open(dropped_path, "w") as fh:
        fh.write("#chr\tstart\tend\tfailure_class\tliftover_reason\n")
        for k in keys:
            if k in failures:
                c, s, e = src[k]
                fh.write(f"{c}\t{s}\t{e}\t{failures[k]}\t{unmapped_reason.get(k, '')}\n")
    per_chrom = defaultdict(Counter)
    per_chrom_total = Counter()
    for k, cls in failures.items():
        c = int(src[k][0])
        per_chrom[c][cls] += 1
    for c, _, _ in windows:
        per_chrom_total[int(c)] += 1
    short = {"unmapped_start": "unm_s", "unmapped_end": "unm_e", "unmapped_both": "unm_se",
             "split": "split", "chromosome_changed": "chr_chg",
             "length_out_of_tolerance": "len_oot", "interior_deleted": "int_del"}
    print("DROPPED WINDOWS BY CHROMOSOME")
    print("  " + f"{'chr':>3} {'windows':>9} {'dropped':>8} {'%':>7}  " +
          " ".join(f"{short[c]:>7}" for c in FAILURE_ORDER))
    for c in sorted(per_chrom_total):
        d = sum(per_chrom[c].values())
        print("  " + f"{c:>3} {per_chrom_total[c]:>9,} {d:>8,} {100.0 * d / per_chrom_total[c]:6.3f}%  " +
              " ".join(f"{per_chrom[c].get(cls, 0):>7,}" for cls in FAILURE_ORDER))
    print()

    runs = defaultdict(list)
    for k in keys:
        if k in failures:
            c, s, e = src[k]
            cls = failures[k]
            cur = runs[int(c)]
            if cur and cur[-1][2] == cls and s <= cur[-1][1]:
                cur[-1] = (cur[-1][0], max(cur[-1][1], e), cls, cur[-1][3] + 1)
            else:
                cur.append((s, e, cls, 1))
    print(f"DROPPED WINDOWS BY REGION (contiguous runs of one class; "
          f"full list: {dropped_path.relative_to(work)})")
    print(f"  {'chr':>3}  {'GRCh37 region':<26} {'class':<24} {'windows':>7}")
    for c in sorted(runs):
        for s, e, cls, n in runs[c]:
            print(f"  {c:>3}  {f'{s:,}-{e:,}':<26} {cls:<24} {n:>7,}")
    if not runs:
        print("  (none)")
    print()

    ships = fail_pct <= 2.0
    print("=" * 78)
    print(f"DECISION: {fail_pct:.3f}% of windows fail to lift cleanly vs the ~2% descope trigger "
          f"-> windows product {'SHIPS (hpa1.12 proceeds)' if ships else 'is DROPPED'}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
