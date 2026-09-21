#!/usr/bin/env python3
"""Split the deCODE 2021 plasma pQTL summary statistics into one munged sumstat per aptamer.

    python3 scripts/munge_decode_pqtl.py --output-dir /mnt/disks/data/decode/munged \
        --stage gs://finngen-commons/results_api_data/sumstats/deCODE_pQTL_2021/ \
        --input-array external_sumstats_input.decode.tsv

The per-aptamer files exist so that the FinnGen autoreporting tool, which takes one
tabix-indexed sumstat per phenotype, can LD-clump each aptamer's associations; those
reports then feed `wdl/create_pseudo_credible_sets.wdl`. See
`docs/pseudo-credible-sets.md` for the deCODE run.

Source file
-----------
`gs://finngen-commons/decode/deCODE_pQTLs_NatGen2021_aligned_p0.005.tsv.gz` (+ `.tbi`,
`.probes`): Ferkingstad et al. 2021 (Nat Genet 53:1712), 35,559 Icelanders, 4,907 SomaScan
v4 aptamers, restricted to p < 0.005. Tab separated with header

    #resource dataset data_type trait chr pos ref alt mlog10p beta se

`trait` is the aptamer id (`seq.10000.28`). The file was derived in 2023 from the sibling
`deCODE_pQTLs_NatGen2021_p0.005.tsv.gz` (`#pheno chr pos ref alt maf beta sebeta mlogp n`)
by normalising alleles to the gnomAD `ref`/`alt` representation; the sibling is not used
here because the aligned file is the cleaned one, but it is what the format assumptions
below were checked against.

Format assumptions (checked against the delivered bytes)
--------------------------------------------------------
- positions are GRCh38 (the paper's release build); contigs are `1`..`22` and `X`, no `chr`
  prefix, sorted by contig then position (the file is tabix indexed on `chr`/`pos`);
- `ref`/`alt` already follow gnomAD and `beta` refers to `alt`: rows whose orientation the
  alignment did not change carry the source's beta unchanged, and rows it re-represented
  (deCODE writes some SNPs with padded alleles, `AAGT`/`GAGT`, and indels in its own
  left-alignment) carry the same beta under the trimmed alleles. The 7 % of source rows
  the alignment could not place (chr21: 9,564,833 → 8,863,739) are not in this file;
- there is no allele frequency: the source carries only `maf`, which cannot be oriented
  to `alt`, so the output has no `af` column and downstream `aaf` is NA. There is no
  sample size either (the source's `n` is per row, ~35,3xx);
- `mlog10p` is finite and never 0 or NA in the p < 0.005 subset (chr21 minimum 2.3008,
  maximum 21,686 — cis-pQTL signals far past any p-value's float range, which is why the
  input is -log10(p) and stays so). `beta` and `se` are never 0 or NA;
- **duplicates**: 1,831 of the 783,081,834 rows share (trait, chr, pos, ref, alt) with a
  sibling row carrying slightly different statistics — two source representations of one
  indel that normalise to the same gnomAD variant. One row per variant is kept, the one
  with the larger `mlog10p`; `--summary` counts them per aptamer;
- per aptamer 31k–491k rows (median 157k); the strongest association per aptamer has
  `mlog10p` from 6.4 to 29,138 (median 47), and 4,844 of the 4,907 reach p ≤ 5e-8.

Output
------
`<output-dir>/<aptamer>.munged.tsv.gz` + `.tbi` (`tabix -s1 -b2 -e2`), columns

    #chr pos ref alt mlog10p beta se

with `X` mapped to `23`, `mlog10p` rounded to 4 decimals and `beta`/`se` formatted `.3e`.
No `mlog10p > 4` companion is written: autoreporting reads the full file, and the API
does not serve these.

Why this script writes its own output
-------------------------------------
The product is a directory of 4,907 files cut from one 12 GB stream, not one DataFrame.
Pass 1 reads the input contig by contig through `tabix` (so `--jobs` workers read
different contigs at once) and appends each row to one of `--buckets` plain-text files
chosen by a hash of the aptamer id, under `<split-dir>/<bucket>/<contig>.tsv`. Pass 2
takes one bucket at a time: reads its contig files in contig order into polars, splits by
aptamer, and for each aptamer deduplicates, formats and hands the frame to
`sumstat_utils.write_bgzip`, which owns the bgzip + tabix step. The bucket is deleted once
every aptamer in it is written.

Buckets exist because the first version wrote 4,907 per-aptamer plain files directly from
one stream and managed 60k rows/s against 290k in isolation: a write buffer per file kept
1.2 GB of half-filled buffers, and the page cache flushed them as interleaved small writes.
Sixty-four sequential streams per worker cost nothing either way; the whole run takes
about an hour on four cores. Disk: the plain text of
the whole input (~2× the compressed size) next to the compressed output, shrinking as
buckets finish; memory: one bucket's rows in polars per worker (~1 GB at 64 buckets).

Flags
-----
--input         aligned deCODE file, gs:// (fetched once into --cache-dir, its .tbi with
                it) or local; `tabix` needs the index either way
--cache-dir     where a gs:// input is cached (default /mnt/disks/data/decode)
--output-dir    local directory for the per-aptamer files (default <cache-dir>/munged)
--split-dir     pass-1 bucket directory (default <output-dir>.split)
--buckets       number of pass-1 buckets (default 64)
--skip-split    do not read the input again: finish pass 2 for whatever buckets are still
                in --split-dir (none, after a complete run) and reuse --summary for the
                rest, so --stage and --input-array can be run on a finished directory
--probes        comma-separated aptamer ids or a file of ids: munge only these (the
                whole input is still read); for a test run
--jobs          workers, for both passes (default 4)
--stage         gs:// prefix to rsync the finished directory to (opt-in, never default)
--input-array   write the autoreporting input array here: one line per aptamer that has
                at least one row at or above --lead-threshold, pointing at the staged
                path when --stage is given and at the local file otherwise
--lead-threshold  autoreporting's `sign_treshold`; an aptamer with no row that
                significant cannot form a locus, so it is left out of the array (default 5e-8)
--summary       per-aptamer row counts, duplicates dropped and max mlog10p
                (default <output-dir>.summary.tsv)
"""

import argparse
import math
import shutil
import subprocess
import sys
import zlib
from multiprocessing import Pool
from pathlib import Path

import polars as pl

from sumstat_utils import fetch_gs, write_bgzip

INPUT_DEFAULT = "gs://finngen-commons/decode/deCODE_pQTLs_NatGen2021_aligned_p0.005.tsv.gz"
CACHE_DIR_DEFAULT = "/mnt/disks/data/decode"
EMPTY_FILE = "gs://misc-analysis/EMPTY_FILE"

INPUT_HEADER = ["#resource", "dataset", "data_type", "trait", "chr", "pos", "ref", "alt", "mlog10p", "beta", "se"]
SPLIT_COLUMNS = ["trait", "chr", "pos", "ref", "alt", "mlog10p", "beta", "se"]
OUTPUT_COLUMNS = ["#chr", "pos", "ref", "alt", "mlog10p", "beta", "se"]
SUMMARY_SCHEMA = {
    "aptamer": pl.Utf8, "rows_in": pl.Int64, "rows_out": pl.Int64, "null_chr": pl.Int64,
    "duplicates": pl.Int64, "max_mlog10p": pl.Float64,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--input", default=INPUT_DEFAULT)
    parser.add_argument("--cache-dir", default=CACHE_DIR_DEFAULT)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--split-dir", default=None)
    parser.add_argument("--buckets", type=int, default=64)
    parser.add_argument("--skip-split", action="store_true")
    parser.add_argument("--probes", default=None)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--stage", default=None)
    parser.add_argument("--input-array", default=None)
    parser.add_argument("--lead-threshold", type=float, default=5e-8)
    parser.add_argument("--summary", default=None)
    args = parser.parse_args()
    args.output_dir = Path(args.output_dir or f"{args.cache_dir}/munged")
    args.split_dir = Path(args.split_dir or f"{args.output_dir}.split")
    args.summary = Path(args.summary or f"{args.output_dir}.summary.tsv")
    if args.stage and not args.stage.startswith("gs://"):
        parser.error("--stage takes a gs:// prefix")
    return args


def read_probe_filter(spec: str | None) -> set[str] | None:
    if spec is None:
        return None
    if Path(spec).exists():
        return {line.strip() for line in Path(spec).read_text().splitlines() if line.strip()}
    return {p.strip() for p in spec.split(",") if p.strip()}


def contigs(input_path: Path) -> list[str]:
    """The input's contigs in file order, which is also the output's sort order once X is 23."""
    out = subprocess.run(["tabix", "-l", str(input_path)], check=True, capture_output=True, text=True)
    return out.stdout.split()


def check_header(input_path: Path) -> int:
    proc = subprocess.Popen(["bgzip", "-d", "-c", str(input_path)], stdout=subprocess.PIPE)
    header = proc.stdout.readline().rstrip(b"\n").decode().split("\t")
    proc.kill()
    if header != INPUT_HEADER:
        raise SystemExit(f"unexpected input header {header}, expected {INPUT_HEADER}")
    return INPUT_HEADER.index("trait")


def split_contigs(job: tuple[Path, list[str], Path, int, set[str] | None, int]) -> int:
    """Pass 1 for one worker: its contigs, each into `<split_dir>/<bucket>/<contig>.tsv`."""
    input_path, my_contigs, split_dir, buckets, keep, trait_idx = job
    keep_bytes = {k.encode() for k in keep} if keep is not None else None
    n = 0
    for contig in my_contigs:
        handles = [open(split_dir / f"{b:03d}" / f"{contig}.tsv", "wb", buffering=1 << 20) for b in range(buckets)]
        proc = subprocess.Popen(["tabix", str(input_path), contig], stdout=subprocess.PIPE)
        for line in proc.stdout:
            # columns before `trait` are constant for the whole file and are not carried
            parts = line.split(b"\t", trait_idx + 1)
            trait = parts[trait_idx]
            if keep_bytes is not None and trait not in keep_bytes:
                continue
            handles[zlib.crc32(trait) % buckets].write(trait + b"\t" + parts[trait_idx + 1])
            n += 1
        for handle in handles:
            handle.close()
        if proc.wait() != 0:
            raise SystemExit(f"tabix failed on contig {contig}")
        print(f"  contig {contig} split", file=sys.stderr, flush=True)
    return n


def split_input(input_path: Path, split_dir: Path, buckets: int, jobs: int, keep: set[str] | None) -> None:
    """Pass 1: every row into a bucket file per contig, contigs spread over `jobs` workers."""
    if split_dir.exists():
        shutil.rmtree(split_dir)
    for b in range(buckets):
        (split_dir / f"{b:03d}").mkdir(parents=True)
    trait_idx = check_header(input_path)
    all_contigs = contigs(input_path)
    # contigs come largest first, so dealing them round-robin balances the workers
    per_worker = [all_contigs[i::jobs] for i in range(jobs)]
    with Pool(jobs) as pool:
        counts = pool.map(split_contigs, [(input_path, c, split_dir, buckets, keep, trait_idx) for c in per_worker])
    print(f"  pass 1: {sum(counts)} rows into {buckets} buckets", file=sys.stderr)


def munge_aptamer(df: pl.DataFrame, out: Path) -> dict:
    """Dedup, map X→23, format and write one aptamer's rows."""
    n_in = df.height
    df = df.with_columns(
        pl.col("chr").str.replace(r"(?i)^chr", "").str.replace(r"^X$", "23").cast(pl.Int32, strict=False).alias("#chr")
    )
    n_null_chr = df["#chr"].null_count()
    df = (
        df.drop_nulls("#chr")
        .sort(["#chr", "pos", "ref", "alt", "mlog10p"], descending=[False, False, False, False, True])
        .unique(subset=["#chr", "pos", "ref", "alt"], keep="first", maintain_order=True)
    )
    df = df.with_columns(
        pl.col("mlog10p").round(4),
        pl.col("beta").map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8),
        pl.col("se").map_elements(lambda x: f"{x:.3e}", return_dtype=pl.Utf8),
    ).select(OUTPUT_COLUMNS)
    write_bgzip(df, str(out))
    return {
        "aptamer": out.name.removesuffix(".munged.tsv.gz"),
        "rows_in": n_in,
        "rows_out": df.height,
        "null_chr": n_null_chr,
        "duplicates": n_in - n_null_chr - df.height,
        "max_mlog10p": df["mlog10p"].max(),
    }


def munge_bucket(job: tuple[Path, list[str], Path, set[str] | None]) -> list[dict]:
    """Pass 2 for one bucket: its contig files in contig order, one output per aptamer."""
    bucket_dir, all_contigs, output_dir, keep = job
    frames = []
    for contig in all_contigs:
        path = bucket_dir / f"{contig}.tsv"
        if path.exists() and path.stat().st_size > 0:
            frames.append(pl.read_csv(
                path, separator="\t", has_header=False, new_columns=SPLIT_COLUMNS,
                schema_overrides={"trait": pl.Utf8, "chr": pl.Utf8, "pos": pl.Int64, "ref": pl.Utf8, "alt": pl.Utf8,
                                  "mlog10p": pl.Float64, "beta": pl.Float64, "se": pl.Float64},
            ))
    results = []
    if frames:
        for (trait,), part in pl.concat(frames).partition_by("trait", as_dict=True).items():
            if keep is not None and trait not in keep:
                continue
            # the aptamer id becomes a file name and, in the API, a path component
            if "/" in trait or trait in (".", ".."):
                raise SystemExit(f"aptamer id {trait!r} is not safe to use as a file name")
            results.append(munge_aptamer(part.drop("trait"), output_dir / f"{trait}.munged.tsv.gz"))
    shutil.rmtree(bucket_dir)
    return results


def main() -> None:
    args = parse_args()
    keep = read_probe_filter(args.probes)
    lead_mlog10p = -math.log10(args.lead_threshold)

    input_path = fetch_gs(args.input, Path(args.cache_dir))
    fetch_gs(args.input + ".tbi", Path(args.cache_dir))
    if args.skip_split:
        if not args.split_dir.is_dir() and not args.summary.exists():
            raise SystemExit(f"--skip-split given but neither {args.split_dir} nor {args.summary} exists")
    else:
        split_input(input_path, args.split_dir, args.buckets, args.jobs, keep)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    all_contigs = contigs(input_path)
    buckets = sorted(p for p in args.split_dir.glob("*") if p.is_dir()) if args.split_dir.is_dir() else []
    print(f"  pass 2: {len(buckets)} buckets with {args.jobs} workers", file=sys.stderr)
    results = []
    with Pool(args.jobs) as pool:
        for i, bucket_results in enumerate(
            pool.imap_unordered(munge_bucket, [(b, all_contigs, args.output_dir, keep) for b in buckets]), 1
        ):
            results.extend(bucket_results)
            print(f"  {i}/{len(buckets)} buckets done, {len(results)} aptamers written", file=sys.stderr, flush=True)

    # the summary persists across runs: an aptamer written earlier keeps its row unless this
    # run wrote it again, so a --skip-split rerun sees the whole directory
    done = pl.DataFrame(results, schema=SUMMARY_SCHEMA)
    if args.summary.exists():
        earlier = pl.read_csv(args.summary, separator="\t", schema_overrides=SUMMARY_SCHEMA)
        done = pl.concat([earlier.filter(~pl.col("aptamer").is_in(done["aptamer"])), done])
    summary = done.sort("aptamer")
    summary.write_csv(args.summary, separator="\t")
    # the threshold is applied here rather than stored, so a rerun may change it
    with_lead = summary.filter(pl.col("max_mlog10p") >= lead_mlog10p)
    print(
        f"  {summary['rows_in'].sum()} rows in, {summary['rows_out'].sum()} out, "
        f"{summary['duplicates'].sum()} duplicates dropped, {summary['null_chr'].sum()} rows with unmappable chr; "
        f"{with_lead.height} of {summary.height} aptamers reach p <= {args.lead_threshold:g}",
        file=sys.stderr,
    )
    print(f"  summary: {args.summary}", file=sys.stderr)

    if args.stage:
        prefix = args.stage.rstrip("/") + "/"
        subprocess.run(["gcloud", "storage", "rsync", "--recursive", str(args.output_dir), prefix], check=True)
        print(f"  staged {args.output_dir} to {prefix}", file=sys.stderr)

    if args.input_array:
        prefix = args.stage.rstrip("/") + "/" if args.stage else str(args.output_dir.resolve()) + "/"
        with open(args.input_array, "w") as out:
            for aptamer in with_lead["aptamer"]:
                out.write(f"{aptamer}\t{prefix}{aptamer}.munged.tsv.gz\t{EMPTY_FILE}\t{EMPTY_FILE}\n")
        print(f"  autoreporting input array with {with_lead.height} aptamers: {args.input_array}", file=sys.stderr)


if __name__ == "__main__":
    main()
