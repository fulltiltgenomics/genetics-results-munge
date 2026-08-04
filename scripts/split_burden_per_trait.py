#!/usr/bin/env python3
"""Split a combined gene burden / exome result file into one tabixed file per trait.

Input is a bgzipped TSV whose first line is the `#`-prefixed header, e.g. the
unfiltered genebass gene burden export or any `.munged.tsv.gz` written by
`write_exome_output()`. Output is `<output-dir>/<trait>.tsv.gz` plus a `.tbi`
for every distinct value of the trait column, each carrying the header and
sorted so tabix can index it.

The input does NOT need to be sorted. The genebass export is row-major (all
traits of one gene x annotation, then the next), which is exactly the order
that makes a single streaming split cheap, so each per-trait file is sorted
here instead — 76k rows per trait rather than 343M rows globally.

Two passes: the first writes one plain-gzip temp file per trait (all handles
stay open, so each input row is written once), the second sorts each temp file
on the tabix columns, bgzips it, indexes it and uploads it if --output-dir is
a gs:// path.

Usage:
  scripts/split_burden_per_trait.py \
    --input /mnt/disks/data/genebass_gene_burden/gene_burden_results.tsv.bgz \
    --output-dir gs://finngen-commons/results_api_data/exome_results/genebass/gene_burden_per_trait

Gene burden files are indexed on the gene locus (`-s5 -b6 -e6`, the default);
pass `--tabix-args "-s2 -b3 -e3"` and `--trait-col 19` for exome variant files.
"""

import argparse
import gzip
import re
import resource
import shutil
import subprocess
import tempfile
from pathlib import Path

# a trait name becomes a file name and, in the API, a path component appended to a
# configured prefix; anything outside this set could escape the output directory
SAFE_TRAIT = re.compile(r"^[A-Za-z0-9._|+-]+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Combined bgzipped TSV to split")
    parser.add_argument("--output-dir", required=True, help="Output directory (local or gs://)")
    parser.add_argument("--trait-col", type=int, default=16,
                        help="1-based column holding the trait key (default: 16, trait_original in the gene burden layout)")
    parser.add_argument("--tabix-args", default="-s5 -b6 -e6",
                        help="tabix -s/-b/-e args; also define the sort key (default: gene locus)")
    parser.add_argument("--suffix", default=".tsv.gz", help="Output file suffix (default: .tsv.gz)")
    parser.add_argument("--tmp-dir", help="Where to keep the intermediate per-trait files "
                                          "(default: next to a local --output-dir, else system temp). "
                                          "The genebass split needs ~11 GB here, which /tmp usually does not have")
    return parser.parse_args()


def _sort_columns(tabix_args: list[str]) -> tuple[int, int, int]:
    """Return 0-based (seq, begin, end) column indices from tabix -s/-b/-e args."""
    flags = {}
    for arg in tabix_args:
        m = re.fullmatch(r"-([sbe])(\d+)", arg)
        if m:
            flags[m.group(1)] = int(m.group(2)) - 1
    missing = {"s", "b", "e"} - set(flags)
    if missing:
        raise ValueError(f"tabix args {tabix_args} do not define {sorted(missing)}")
    return flags["s"], flags["b"], flags["e"]


def split_by_trait(input_path: str, tmp_dir: Path, trait_idx: int) -> tuple[dict[str, Path], int, bytes]:
    """Stream the input into one gzipped temp file per trait. Returns paths, row count and header."""
    handles: dict[str, gzip.GzipFile] = {}
    paths: dict[str, Path] = {}
    n_rows = 0

    # gzip, not bgzip -d: BGZF is gzip-compatible, and this also accepts a plain gzip input
    with subprocess.Popen(["gzip", "-cd", input_path], stdout=subprocess.PIPE) as proc:
        stream = proc.stdout
        header = stream.readline()
        if not header.startswith(b"#"):
            raise ValueError(f"{input_path} does not start with a #-prefixed header line")
        for line in stream:
            trait = line.rstrip(b"\n").split(b"\t")[trait_idx].decode()
            handle = handles.get(trait)
            if handle is None:
                if not SAFE_TRAIT.match(trait):
                    raise ValueError(f"trait {trait!r} is not safe to use as a file name")
                path = tmp_dir / f"{trait}.part.gz"
                # level 1: these are read back within minutes, so speed beats ratio
                handle = gzip.open(path, "wb", compresslevel=1)
                handles[trait] = handle
                paths[trait] = path
            handle.write(line)
            n_rows += 1
            if n_rows % 20_000_000 == 0:
                print(f"  {n_rows:,} rows into {len(handles):,} traits")
        if proc.wait() != 0:
            raise RuntimeError(f"gzip -cd failed on {input_path}")

    for handle in handles.values():
        handle.close()
    return paths, n_rows, header


def write_trait_file(part: Path, out_path: Path, header: bytes,
                     sort_cols: tuple[int, int, int], tabix_args: list[str]) -> tuple[int, int]:
    """Sort one trait's rows on the tabix columns, bgzip and index. Returns (written, dropped)."""
    seq_i, beg_i, end_i = sort_cols
    rows = []
    dropped = 0
    with gzip.open(part, "rb") as f:
        for line in f:
            fields = line.rstrip(b"\n").split(b"\t")
            try:
                key = (int(fields[seq_i]), int(fields[beg_i]), int(fields[end_i]))
            except ValueError:
                # tabix cannot index a row without a locus (gene missing from the gencode join)
                dropped += 1
                continue
            rows.append((key, line))
    rows.sort(key=lambda r: r[0])

    with open(out_path, "wb") as out, subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=out) as proc:
        proc.stdin.write(header)
        for _, line in rows:
            proc.stdin.write(line)
        proc.stdin.close()
        if proc.wait() != 0:
            raise RuntimeError(f"bgzip failed for {out_path}")
    subprocess.run(["tabix", "-f"] + tabix_args + [str(out_path)], check=True)
    return len(rows), dropped


def main() -> None:
    args = parse_args()
    tabix_args = args.tabix_args.split()
    sort_cols = _sort_columns(tabix_args)
    trait_idx = args.trait_col - 1
    is_gcs = args.output_dir.startswith("gs://")

    # one open handle per trait for the whole first pass (~4.5k for genebass)
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (min(hard, max(soft, 16384)), hard))

    # the intermediates are as big as the output, so default them next to it rather
    # than onto whatever /tmp happens to be mounted on
    tmp_dir = args.tmp_dir
    if tmp_dir is None and not is_gcs:
        tmp_dir = str(Path(args.output_dir).resolve().parent)
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    tmp_root = Path(tempfile.mkdtemp(dir=tmp_dir, prefix="per_trait_"))
    out_dir = tmp_root / "out" if is_gcs else Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        print(f"Splitting {args.input} by column {args.trait_col}...")
        paths, n_rows, header = split_by_trait(args.input, tmp_root, trait_idx)
        print(f"  {n_rows:,} rows in {len(paths):,} traits")

        print(f"Sorting, bgzipping and indexing {len(paths):,} per-trait files...")
        total_written = total_dropped = 0
        for i, (trait, part) in enumerate(sorted(paths.items()), start=1):
            out_path = out_dir / f"{trait}{args.suffix}"
            written, dropped = write_trait_file(part, out_path, header, sort_cols, tabix_args)
            total_written += written
            total_dropped += dropped
            part.unlink()
            if i % 500 == 0:
                print(f"  {i:,}/{len(paths):,} traits")

        print(f"Done! {total_written:,} rows in {len(paths):,} per-trait files")
        if total_dropped:
            print(f"  dropped {total_dropped:,} rows with a non-numeric locus (not indexable)")

        if is_gcs:
            # one rsync rather than two gcloud invocations per trait; it does not
            # delete anything at the destination unless asked to
            dest = args.output_dir.rstrip("/")
            print(f"Uploading {len(paths) * 2:,} files to {dest}/...")
            subprocess.run(["gcloud", "storage", "rsync", str(out_dir), dest], check=True)
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)


if __name__ == "__main__":
    main()
