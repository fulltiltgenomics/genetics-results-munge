#!/usr/bin/env python3
"""Build the rsID -> variant tabix file behind results-api's `/rsid/variants`.

Source file
-----------
The gnomAD v4.0 genomes+exomes sites table served by results-api
(`gnomad.genomes.exomes.v4.0.sites.v2.tsv.bgz`, bgzip + tabix), tab separated with header

    #chr pos ref alt rsids filters AN AF AF_afr ...

Only the first six columns are read.

Output
------
`<output>` (bgzip) and `<output>.csi`, tab separated with header

    #rs id rsid chr pos ref alt

one row per (rsID, allele), sorted by `id` (the rs number) and indexed on a pseudo-contig
`rs` with `id` as the position, so `tabix --csi <output> rs:N-N` answers "which variants
carry rsN". results-api's `RsidDB` parses columns 4-7 of every returned row, so the column
order is a contract. `chr` is coded the way the served files are: X is 23, Y is 24.

Format assumptions (checked against the delivered sites file)
-------------------------------------------------------------
- `rsids` is `NA` or a comma-separated list of rs ids; a few in ten thousand rows carry two;
- an allele present in both genomes and exomes appears as two rows, which may differ only
  in `filters`;
- `filters` is empty or a comma-separated subset of `AC0`, `AS_VQSR`, `InbreedingCoeff`;
- contigs are 1-22, X and Y; a contig that is still non-numeric after the X/Y mapping is
  dropped and counted, because results-api would reject it anyway.

One row per allele, and the AC0 rule
------------------------------------
A dbSNP id names a position, not an allele, so one rsID can carry several alt alleles. The
previous build kept one row per rsID — the alphabetically first alt — which for a
multi-allelic rsID was as likely as not an allele with a handful of carriers, or none:
rs16890065 resolved to 8:40888741:C:A (AF 2e-5) while the variant people mean is C:T
(AF 0.13). Every allele is emitted; `RsidDB` already returns a list per rsID.

Alleles whose every row is filtered `AC0` (no carriers passing QC in gnomAD) are dropped
when the same rsID has an allele that is not, so a ghost allele does not sit next to the
real one. An rsID whose alleles are all AC0 keeps them, so no rsID that resolved before
stops resolving. `--keep-ac0` disables the rule.

Two passes, because the rows must end up sorted by rs number and the input is far larger
than memory: pass 1 streams the sites file once and partitions exploded rows into gzip
buckets of `--bucket-width` consecutive rs numbers; pass 2 sorts each bucket with GNU sort
under `--sort-memory`, applies the rules in one streaming pass and appends it to a single
bgzip stream. Nothing holds more than one rsID's alleles in Python. Peak disk is roughly
the gzip size of the exploded rows plus the output, since each bucket is deleted once
written; peak memory is `--sort-memory` plus the bgzip and gzip buffers, a little over 1 GB
at the defaults. Run it under a memory cap anyway on a shared machine, e.g.
`systemd-run --user --scope -p MemoryMax=4G nice -n 10 python scripts/...`.

Usage
-----
    python scripts/build_gnomad_rsid_index.py \\
        --gnomad gs://daly-genetics-results/gnomad/gnomad.genomes.exomes.v4.0.sites.v2.tsv.bgz \\
        --output /mnt/disks/data/gnomad/gnomad.genomes.exomes.v4.0.rsid.v2.tsv.gz \\
        --upload gs://daly-genetics-results/gnomad/

`--gnomad` may be a local path or `gs://`, streamed in retried byte ranges with the
credentials of `gcloud auth`. `--upload` copies the output and its `.csi` to that
directory when the build and its self-check pass.
Upload under a new name rather than over the served file: results-api caches the parsed
index per path and treats the file as immutable, so a file replaced in place is read with
a stale index until every pod restarts.
"""

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

BUCKET_PREFIX = "rsid_bucket_"
HEADER = b"#rs\tid\trsid\tchr\tpos\tref\talt\n"

# the rs number is kept as the string `substr(id, 3)` and only used numerically for the
# bucket index: mawk prints a numeric value through CONVFMT (%.6g), which would turn
# rs1638893388 into 1.63889e+09
PASS1_AWK = r"""
NR > 1 && $5 != "NA" {
    chr = $1
    if (chr == "X") chr = 23
    else if (chr == "Y") chr = 24
    if (chr !~ /^[0-9]+$/) { skipped++; next }
    ac0 = (index($6, "AC0") > 0) ? 1 : 0
    n = split($5, ids, ",")
    for (i = 1; i <= n; i++) {
        num = substr(ids[i], 3)
        f = sprintf("%s/%s%05d.gz", D, P, int(num / W))
        print num, chr, $2, $3, $4, ac0 | ("gzip -1 > " f)
    }
}
END { printf "skipped %d rows on non-numeric contigs\n", skipped > "/dev/stderr" }
"""


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def gcs_token() -> str:
    return subprocess.run(
        ["gcloud", "auth", "print-access-token"], capture_output=True, text=True, check=True
    ).stdout.strip()


def gcs_object_url(gs_path: str) -> str:
    bucket, name = gs_path[len("gs://"):].split("/", 1)
    return f"https://storage.googleapis.com/storage/v1/b/{bucket}/o/{urllib.parse.quote(name, safe='')}"


def stream_gcs(gs_path: str, out, chunk_bytes: int = 256 << 20) -> None:
    """Copy a GCS object to `out` in retried byte ranges.

    A single `gcloud storage cat` of this 31 GB file ran at ~135 MB/s (`gsutil cat` at ~6)
    but got a connection reset a few minutes in, and a reset mid-stream cannot be resumed
    into the pipeline downstream. Ranges make every failure local to one chunk. The token
    is re-read on 401 rather than on a timer.
    """
    url = gcs_object_url(gs_path)
    token = gcs_token()
    meta = urllib.request.Request(url, headers={"Authorization": f"Bearer {token}"})
    size = int(json.load(urllib.request.urlopen(meta, timeout=60))["size"])
    pos = 0
    while pos < size:
        end = min(pos + chunk_bytes, size) - 1
        # a retry resumes at the first byte not yet written: bytes already handed to `out`
        # cannot be taken back, so re-requesting the whole range would duplicate them
        for attempt in range(8):
            req = urllib.request.Request(
                url + "?alt=media",
                headers={"Authorization": f"Bearer {token}", "Range": f"bytes={pos}-{end}"},
            )
            try:
                with urllib.request.urlopen(req, timeout=120) as resp:
                    while True:
                        buf = resp.read(8 << 20)
                        if not buf:
                            break
                        out.write(buf)
                        pos += len(buf)
                if pos != end + 1:
                    raise OSError(f"short read, stopped at byte {pos} of {end + 1}")
                break
            except (urllib.error.HTTPError, OSError) as e:
                if isinstance(e, urllib.error.HTTPError) and e.code == 401:
                    token = gcs_token()
                log(f"  range {pos}-{end} attempt {attempt + 1} failed: {e}")
                time.sleep(2 * (attempt + 1))
        else:
            raise OSError(f"giving up on {gs_path} at byte {pos}")


def partition(gnomad: str, tmp_dir: str, bucket_width: int, threads: int) -> None:
    os.makedirs(tmp_dir, exist_ok=True)
    stale = [f for f in os.listdir(tmp_dir) if f.startswith(BUCKET_PREFIX)]
    if stale:
        sys.exit(f"{tmp_dir} already holds {len(stale)} bucket files; remove them or pick another --tmp-dir")
    awk = (
        f"mawk -F '\\t' -v OFS='\\t' -v W={bucket_width} -v D={shlex.quote(tmp_dir)} "
        f"-v P={BUCKET_PREFIX} {shlex.quote(PASS1_AWK)}"
    )
    # bgzip's threaded decompression is ~4x zcat on this file, which was the next bottleneck
    cmd = f"set -o pipefail; bgzip -d -c -@{threads} | cut -f1-6 | {awk}"
    log("pass 1: partitioning into buckets")
    pipeline = subprocess.Popen(["bash", "-c", cmd], stdin=subprocess.PIPE)
    try:
        if gnomad.startswith("gs://"):
            stream_gcs(gnomad, pipeline.stdin)
        else:
            with open(gnomad, "rb") as f:
                shutil.copyfileobj(f, pipeline.stdin, 8 << 20)
    finally:
        pipeline.stdin.close()
    if pipeline.wait() != 0:
        sys.exit("pass 1 pipeline failed")


def sorted_bucket(path: str, tmp_dir: str, sort_memory: str):
    """Yield the bucket's rows in (rs number, chr, pos, ref, alt) order.

    The sort is GNU sort with a fixed buffer, spilling to `tmp_dir`, not a Python list: the
    rs-number space is dense around the ids dbSNP assigned most recently, so one bucket can
    hold tens of millions of rows, and materialising that as tuples took the machine down.
    """
    cmd = (
        f"gzip -dc {shlex.quote(path)} | LC_ALL=C sort -t '\t' -k1,1n -k2,2n -k3,3n -k4,4 -k5,5 "
        f"-S {shlex.quote(sort_memory)} -T {shlex.quote(tmp_dir)} --compress-program=gzip"
    )
    proc = subprocess.Popen(["bash", "-c", f"set -o pipefail; {cmd}"], stdout=subprocess.PIPE)
    for line in proc.stdout:
        num, chr_, pos, ref, alt, ac0 = line.rstrip(b"\n").split(b"\t")
        yield int(num), (int(chr_), int(pos), ref.decode(), alt.decode()), ac0 == b"1"
    if proc.wait() != 0:
        sys.exit(f"sorting {path} failed")


def emit_rsid(out, num: int, alleles: dict, keep_ac0: bool) -> int:
    """`alleles` maps (chr, pos, ref, alt) -> True when every source row was AC0."""
    if not keep_ac0 and not all(alleles.values()):
        alleles = {k: v for k, v in alleles.items() if not v}
    for chr_, pos, ref, alt in sorted(alleles):
        out.write(f"rs\t{num}\trs{num}\t{chr_}\t{pos}\t{ref}\t{alt}\n".encode())
    return len(alleles)


def merge(tmp_dir: str, output: str, keep_ac0: bool, threads: int, sort_memory: str) -> tuple[int, int]:
    buckets = sorted(f for f in os.listdir(tmp_dir) if f.startswith(BUCKET_PREFIX))
    log(f"pass 2: merging {len(buckets)} buckets into {output}")
    n_rsids = n_rows = 0
    with open(output, "wb") as fh:
        bgzip = subprocess.Popen(["bgzip", "-c", f"-@{threads}"], stdin=subprocess.PIPE, stdout=fh)
        bgzip.stdin.write(HEADER)
        for i, name in enumerate(buckets, 1):
            path = os.path.join(tmp_dir, name)
            current = None
            alleles: dict = {}
            for num, key, ac0 in sorted_bucket(path, tmp_dir, sort_memory):
                if num != current:
                    if current is not None:
                        n_rsids += 1
                        n_rows += emit_rsid(bgzip.stdin, current, alleles, keep_ac0)
                    current, alleles = num, {}
                alleles[key] = alleles.get(key, True) and ac0
            if current is not None:
                n_rsids += 1
                n_rows += emit_rsid(bgzip.stdin, current, alleles, keep_ac0)
            os.remove(path)
            if i % 20 == 0 or i == len(buckets):
                log(f"  {i}/{len(buckets)} buckets, {n_rsids:,} rsids, {n_rows:,} rows")
        bgzip.stdin.close()
        if bgzip.wait() != 0:
            sys.exit("bgzip failed")
    return n_rsids, n_rows


def index(output: str) -> None:
    log("indexing")
    subprocess.run(["tabix", "-f", "--csi", "-s1", "-b2", "-e2", output], check=True)


def self_check(output: str) -> None:
    """The multi-allelic rsID that started this: both alleles must resolve, and the rs
    numbers above 2^31 that a .tbi cannot address must be reachable through the .csi."""
    for rsid, expect in (("16890065", {"8\t40888741\tC\tA", "8\t40888741\tC\tT"}), ("1638893388", None)):
        r = subprocess.run(
            ["tabix", "--csi", output, f"rs:{rsid}-{rsid}"], capture_output=True, text=True, check=True
        )
        got = {"\t".join(l.split("\t")[3:7]) for l in r.stdout.splitlines()}
        if not got or (expect is not None and got != expect):
            sys.exit(f"self-check failed for rs{rsid}: got {got or 'nothing'}")
    log("self-check passed")


def upload(output: str, dest: str) -> None:
    dest = dest.rstrip("/") + "/"
    log(f"uploading to {dest}")
    subprocess.run(["gsutil", "-q", "cp", output, output + ".csi", dest], check=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--gnomad", required=True, help="gnomAD sites TSV (bgzip), local path or gs://")
    p.add_argument("--output", required=True, help="output .tsv.gz path; the .csi is written next to it")
    p.add_argument("--tmp-dir", help="bucket directory (default: next to --output)")
    p.add_argument("--bucket-width", type=int, default=10_000_000, help="rs numbers per bucket")
    p.add_argument("--keep-ac0", action="store_true", help="keep AC0 alleles next to observed ones")
    p.add_argument("--threads", type=int, default=2, help="bgzip threads, decompressing in pass 1 and compressing in pass 2")
    p.add_argument("--sort-memory", default="512M", help="GNU sort buffer per bucket (its -S)")
    p.add_argument("--upload", help="gs:// directory to copy the output and its .csi into")
    p.add_argument("--skip-partition", action="store_true", help="reuse buckets already in --tmp-dir")
    args = p.parse_args()

    tmp_dir = args.tmp_dir or os.path.join(os.path.dirname(os.path.abspath(args.output)), "rsid_tmp")
    if not args.skip_partition:
        partition(args.gnomad, tmp_dir, args.bucket_width, args.threads)
    n_rsids, n_rows = merge(tmp_dir, args.output, args.keep_ac0, args.threads, args.sort_memory)
    log(f"wrote {n_rsids:,} rsids, {n_rows:,} rows")
    index(args.output)
    self_check(args.output)
    if args.upload:
        upload(args.output, args.upload)
    log("done")


if __name__ == "__main__":
    main()
