#!/usr/bin/env python3
"""Stamp gnomAD consequence onto a served credible-set file without re-munging it.

Overwrites `most_severe` and `gene_most_severe` from the consequence file that
`build_gnomad_annotation.py --consequence-output` writes. A variant the consequence file
does not hold gets `NA` in both. Every other byte of every row is written back as read, the
header included, and the output has the container of the input: bgzip in, bgzip out with a
tabix index built from the settings of the input's own index; plain text in, plain text out.
The container is read off the file's first bytes, not its name.

Inputs
------
- the served file: tab separated, first line a header naming `chr`, `pos`, `ref`, `alt`,
  `most_severe` and `gene_most_severe` (a leading `#` is ignored). A file missing one of
  them is refused. `chr` and `pos` are plain integers, X being 23;
- `--consequence`: bgzip + tabix, sorted by numeric chromosome then position, one row per
  variant, the same six columns located by name. Every row must have the header's
  column count.

Modes
-----
`--mode merge` is for files sorted by variant position. It walks the file and the
consequence file together, holding one position's alleles at a time.

`--mode lookup` is for the gene-indexed QTL copies, which are sorted on the trait gene's
coordinates. A first pass collects the file's distinct variants and stops once there are
more than `--max-lookup-rows`, before anything is read from the consequence file; the
consequence file is then cut down to those variants and a second pass stamps from that
table. Row order is untouched. Both modes write the same bytes for the same rows.

`--clear` writes `NA` to both columns of any file and reads no consequence file.

`--verify ORIGINAL STAMPED` compares two files row by row and exits non-zero unless the
headers are identical and every column other than the two annotation columns is too.

How the consequence file is read
--------------------------------
Always through tabix. `--consequence-access regions` asks for just the neighbourhoods of
the input's positions, which takes time proportional to the input once its positions are
known. `scan` reads each chromosome the input touches from start to end, which takes the
same time whatever the input holds. `auto` is `regions` in lookup mode, where the positions
are collected anyway; in merge mode it is `regions` for an input of at most
`--regions-max-bytes` on disk and `scan` above that, so that a large file is read once and
not twice.

In merge mode both streams must be in order and the run stops, naming the position, when one
is not.

Guards
------
Resident memory is checked against `--max-rss-gb` as rows go by and enforced again on
address space through RLIMIT_AS, which the bgzip and tabix children inherit. Rows written
must equal rows read. Outputs appear under their final names only when the run succeeds.

The last line on stdout is a JSON summary. In it a value counts as missing when it is `NA`
or empty, and `changed_non_na` counts values that were not missing and came out different.

Usage
-----
    python scripts/annotate_consequence.py --consequence gnomad.v4.0.consequence.tsv.bgz \\
        --input PGC_SCZ_2022_credible_sets.tsv.gz --output out/PGC_SCZ_2022_credible_sets.tsv.gz
    python scripts/annotate_consequence.py --consequence gnomad.v4.0.consequence.tsv.bgz \\
        --mode lookup --input X_credible_sets.qtl.tsv.gz --output out/X_credible_sets.qtl.tsv.gz
    python scripts/annotate_consequence.py --clear --input in.tsv --output out.tsv
    python scripts/annotate_consequence.py --verify in.tsv.gz out.tsv.gz
"""

import argparse
import gzip
import json
import math
import os
import resource
import struct
import subprocess
import sys
import time
from array import array
from collections import defaultdict
from itertools import zip_longest

TAB = b"\t"
NA = b"NA"
MISSING = (NA, b"")
COLUMNS = (b"chr", b"pos", b"ref", b"alt", b"most_severe", b"gene_most_severe")
ANNOTATION = [c.decode() for c in COLUMNS[4:]]
NOTHING: dict = {}
# positions closer than this share one tabix region: a seek re-inflates a block and re-reads
# up to an index window of rows, which costs more than reading a few kb of rows it skips
REGION_GAP = 5000
# lookup mode keeps a position as chr << CHROM_SHIFT | pos, eight bytes that sort correctly
CHROM_SHIFT = 32
POS_MASK = (1 << CHROM_SHIFT) - 1
CHECKPOINT_ROWS = 1 << 18
PAGE_SIZE = os.sysconf("SC_PAGE_SIZE")
# files that must not outlive a failed run
SCRATCH: list[str] = []


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def die(msg: str) -> None:
    sys.exit(f"annotate_consequence: {msg}")


def rss_guard(max_rss: int):
    def check() -> None:
        with open("/proc/self/statm") as f:
            rss = int(f.read().split()[1]) * PAGE_SIZE
        if rss > max_rss:
            die(f"resident memory {rss >> 20} MB is over --max-rss-gb")

    return check


def column_indexes(header: bytes, path: str) -> tuple[list[int], int]:
    """Positions of COLUMNS in a header line, and the header's column count."""
    cols = header.rstrip(b"\n").split(TAB)
    cols[0] = cols[0].removeprefix(b"#")
    bad = [c.decode() for c in COLUMNS if cols.count(c) != 1]
    if bad:
        die(f"{path}: header must name each of {', '.join(bad)} exactly once")
    return [cols.index(c) for c in COLUMNS], len(cols)


def is_bgzf(path: str) -> bool:
    with open(path, "rb") as f:
        head = f.read(16)
    if head[:2] != b"\x1f\x8b":
        return False
    if head[3:4] != b"\x04" or head[12:14] != b"BC":
        die(f"{path}: gzip but not bgzip, which tabix cannot index")
    return True


class Layout:
    """Where a served file's rows are split and which fields are read out of the split."""

    def __init__(self, header: bytes, path: str):
        if not header.endswith(b"\n"):
            die(f"{path}: no header line")
        self.indexes, ncols = column_indexes(header, path)
        top = max(self.indexes)
        # a row is split one field past the last column needed, so that a row longer or
        # shorter than the header shows up as a wrong field count instead of being
        # absorbed into the last field
        self.nsplit = top + 1
        self.nfields = min(ncols, top + 2)
        self.last = top if top == ncols - 1 else None
        if self.last is not None and self.last not in self.indexes[4:]:
            die(f"{path}: a key column is the last column, which is not supported")


class Reader:
    """Rows of a served file, bgzip or plain, with the header read off."""

    def __init__(self, path: str, threads: int):
        self.path = path
        self.bgzf = is_bgzf(path)
        self.proc = None
        if self.bgzf:
            self.proc = subprocess.Popen(
                ["bgzip", "-dc", f"-@{threads}", path], stdout=subprocess.PIPE, bufsize=1 << 20
            )
            self.lines = self.proc.stdout
        else:
            self.lines = open(path, "rb", buffering=1 << 20)
        self.header = self.lines.readline()
        self.layout = Layout(self.header, path)

    def close(self) -> None:
        self.lines.close()
        # a truncated file decompresses to fewer rows, which no row count would notice
        if self.proc and self.proc.wait() != 0:
            die(f"bgzip failed reading {self.path}")


class Writer:
    """The output under a scratch name, in the container of the input."""

    def __init__(self, path: str, bgzf: bool, threads: int):
        self.path, self.part = path, path + ".part"
        SCRATCH.append(self.part)
        self.file = open(self.part, "wb", buffering=1 << 20)
        self.proc = None
        self.write = self.file.write
        if bgzf:
            self.proc = subprocess.Popen(
                ["bgzip", "-c", f"-@{threads}"], stdin=subprocess.PIPE, stdout=self.file, bufsize=1 << 20
            )
            self.write = self.proc.stdin.write

    def close(self) -> None:
        if self.proc:
            self.proc.stdin.close()
            if self.proc.wait() != 0:
                die(f"bgzip failed writing {self.part}")
        self.file.close()

    def publish(self, index: "Index | None", threads: int) -> None:
        if index:
            subprocess.run(["tabix", "-f", f"-@{threads}", *index.args, self.part], check=True)
            SCRATCH.append(self.part + index.ext)
            if Index(self.part).settings != index.settings:
                die(f"the index built for {self.path} does not have the settings of the input's")
            os.replace(self.part + index.ext, self.path + index.ext)
        os.replace(self.part, self.path)


class Index:
    """The settings of the tabix index next to a bgzip file, as arguments that rebuild it."""

    def __init__(self, path: str):
        self.ext = next((e for e in (".tbi", ".csi") if os.path.exists(path + e)), None)
        if not self.ext:
            die(f"{path}: no .tbi or .csi index next to it")
        with gzip.open(path + self.ext, "rb") as f:
            raw = f.read(64)
        csi = []
        if raw[:4] == b"TBI\x01":
            conf = struct.unpack_from("<6i", raw, 8)
        elif raw[:4] == b"CSI\x01" and struct.unpack_from("<i", raw, 12)[0] >= 28:
            csi = ["-C", "-m", str(struct.unpack_from("<i", raw, 4)[0])]
            conf = struct.unpack_from("<6i", raw, 16)
        else:
            die(f"{path}{self.ext}: not a tabix index")
        fmt, seq, beg, end, meta, skip = conf
        if fmt & 0xFFFF:
            die(f"{path}{self.ext}: built with a tabix preset, expected explicit columns")
        self.settings = (self.ext, *csi, *conf)
        self.args = [*csi, f"-s{seq}", f"-b{beg}", f"-S{skip}", "-c", chr(meta)]
        if end:
            self.args.append(f"-e{end}")
        # the bit tabix sets for -0, half-open zero-based coordinates
        if fmt & 0x10000:
            self.args.append("-0")


class Stamper:
    """Rewrites the two annotation fields of a split row and counts what changed."""

    def __init__(self, reader: Reader, write):
        lay = reader.layout
        self.path, self.write = reader.path, write
        self.ms, self.gene = lay.indexes[4:]
        self.nfields, self.last = lay.nfields, lay.last
        self.rows = self.matched = 0
        # rows per (old most_severe, old gene, new most_severe, new gene): one count per
        # row instead of several, and the summary is derived from it at the end
        self.transitions = defaultdict(int)

    def stamp(self, f: list, hit) -> None:
        """`hit` is the (most_severe, gene_most_severe) to write, or None for a variant
        the consequence file does not hold."""
        self.rows += 1
        if len(f) != self.nfields:
            die(f"{self.path}: row {self.rows + 1} does not have the header's column count")
        eol = b""
        if self.last is not None and f[self.last].endswith(b"\n"):
            f[self.last] = f[self.last][:-1]
            eol = b"\n"
        if hit:
            self.matched += 1
            ms, gene = hit
        else:
            ms = gene = NA
        self.transitions[f[self.ms], f[self.gene], ms, gene] += 1
        f[self.ms], f[self.gene] = ms, gene
        if eol:
            f[self.last] += eol
        self.write(TAB.join(f))


def transition_counts(transitions: dict) -> dict:
    """Per-column summary counts from the {(old pair, new pair): rows} table."""
    counts = {}
    for i, name in enumerate(ANNOTATION):
        before = after = changed = 0
        for key, n in transitions.items():
            old, new = key[i], key[i + 2]
            before += n * (old in MISSING)
            after += n * (new in MISSING)
            changed += n * (old not in MISSING and old != new)
        counts[name] = {"na_before": before, "na_after": after, "changed_non_na": changed}
    return counts


class RegionFile:
    """A BED of the intervals covering every position added, for `tabix -R`. Positions
    arrive in ascending order."""

    def __init__(self, path: str):
        self.path = path
        SCRATCH.append(path)
        self.bed = open(path, "w")
        self.chrom = None

    def add(self, chrom: int, pos: int) -> None:
        if chrom != self.chrom or pos > self.end + REGION_GAP:
            self.flush()
            self.chrom, self.start = chrom, pos
        self.end = pos

    def flush(self) -> None:
        if self.chrom is not None:
            self.bed.write(f"{self.chrom}\t{self.start - 1}\t{self.end}\n")

    def close(self) -> None:
        self.flush()
        self.chrom = None
        self.bed.close()


class Consequence:
    """The consequence file as a position-ordered stream out of tabix: the regions of a
    RegionFile, or whole chromosomes one at a time."""

    def __init__(self, path: str, by_chrom: bool):
        self.path, self.by_chrom = path, by_chrom
        if not any(os.path.exists(path + e) for e in (".tbi", ".csi")):
            die(f"{path}: no .tbi or .csi index next to it")
        with gzip.open(path, "rb") as f:
            self.indexes, self.ncols = column_indexes(f.readline(), path)
        self.proc = None
        self.chrom = None

    def open(self, *query: str) -> None:
        self.close()
        self.proc = subprocess.Popen(["tabix", *query], stdout=subprocess.PIPE, bufsize=1 << 20)
        self.lines = self.proc.stdout
        self.name = None
        self.held = (None, 0, 0)

    def open_regions(self, regions: RegionFile) -> None:
        regions.close()
        self.open("-R", regions.path, self.path)

    def open_chrom(self, chrom: int) -> None:
        self.chrom = chrom
        self.open(self.path, str(chrom))

    def exhausted(self) -> None:
        # a tabix that died mid-stream would otherwise look like "the rest is absent"
        if self.proc.wait() != 0:
            die(f"tabix failed reading {self.path}")

    def close(self) -> None:
        if self.proc:
            self.lines.close()
            self.proc.kill()
            self.proc.wait()

    def fields(self, line: bytes) -> list:
        f = line[:-1].split(TAB)
        if len(f) != self.ncols:
            die(f"{self.path}: malformed row {line[:200].decode(errors='replace')!r}")
        return f

    def duplicate(self, line: bytes) -> None:
        die(f"{self.path}: more than one row for the variant of {line[:200].decode(errors='replace')!r}")

    def add(self, found: dict, line: bytes) -> None:
        _, _, ri, ai, mi, gi = self.indexes
        f = self.fields(line)
        key = f[ri] + TAB + f[ai]
        if key in found:
            self.duplicate(line)
        found[key] = (f[mi], f[gi])

    def variants_at(self, chrom: int, pos: int) -> dict:
        """{ref TAB alt: (most_severe, gene_most_severe)} of one position. Positions must be
        asked for in ascending order."""
        if self.by_chrom and chrom != self.chrom:
            self.open_chrom(chrom)
        found = {}
        line, cchr, cpos = self.held
        if (cchr, cpos) > (chrom, pos):
            return found
        if line is not None and (cchr, cpos) == (chrom, pos):
            self.add(found, line)
        name, ci, pi = self.name, self.indexes[0], self.indexes[1]
        nsplit = max(ci, pi) + 1
        try:
            for line in self.lines:
                f = line.split(TAB, nsplit)
                p = int(f[pi])
                if f[ci] != name:
                    name = f[ci]
                    c = int(name)
                    if c < cchr:
                        die(f"{self.path}: not sorted at {c}:{p}, after chromosome {cchr}")
                    cchr = c
                elif p < cpos:
                    die(f"{self.path}: not sorted at {cchr}:{p}, after position {cpos}")
                cpos = p
                if cchr == chrom:
                    if p < pos:
                        continue
                    if p == pos:
                        self.add(found, line)
                        continue
                elif cchr < chrom:
                    continue
                self.name, self.held = name, (line, cchr, cpos)
                return found
        except (ValueError, IndexError):
            die(f"{self.path}: malformed row {line[:200].decode(errors='replace')!r}")
        self.exhausted()
        self.held = (None, math.inf, math.inf)
        return found

    def fill(self, table: dict) -> None:
        """Give every key of `table` that the open stream holds its annotation."""
        ci, pi, ri, ai, mi, gi = self.indexes
        # most variants share their pair with many others, so one object per distinct pair
        pairs: dict = {}
        nsplit = max(ci, pi, ri, ai) + 1
        try:
            for line in self.lines:
                f = line.split(TAB, nsplit)
                if TAB.join((f[ci], f[pi], f[ri], f[ai])) in table:
                    f = self.fields(line)
                    key = TAB.join((f[ci], f[pi], f[ri], f[ai]))
                    if table[key] is not None:
                        self.duplicate(line)
                    pair = (f[mi], f[gi])
                    table[key] = pairs.setdefault(pair, pair)
        except IndexError:
            die(f"{self.path}: malformed row {line[:200].decode(errors='replace')!r}")
        self.exhausted()


def plain_integers(path: str, row: int, chrom_field: bytes, pos_field: bytes) -> tuple[int, int]:
    """chr and pos as integers. Another spelling of the same number would compare equal
    here and unequal as text, so the two modes could disagree about it."""
    c, p = int(chrom_field), int(pos_field)
    if c < 1 or not 0 < p <= POS_MASK or chrom_field != b"%d" % c or pos_field != b"%d" % p:
        die(f"{path}: row {row}: chr {chrom_field!r} or pos {pos_field!r} is not a plain positive integer")
    return c, p


def walk(reader: Reader, variants_at, stamp, check) -> int:
    """Drive a position-sorted file: `variants_at` once per position, `stamp` once per row.
    Returns the rows read."""
    ci, pi, ri, ai = reader.layout.indexes[:4]
    nsplit = reader.layout.nsplit
    name = spelled = None
    chrom = pos = 0
    found = NOTHING
    row = 1
    for line in reader.lines:
        row += 1
        f = line.split(TAB, nsplit)
        try:
            if f[pi] != spelled or f[ci] != name:
                c, p = plain_integers(reader.path, row, f[ci], f[pi])
                if (c, p) < (chrom, pos):
                    die(f"{reader.path}: not sorted at {c}:{p} (row {row}), after {chrom}:{pos}")
                name, spelled, chrom, pos = f[ci], f[pi], c, p
                found = variants_at(c, p)
            hit = found.get(f[ri] + TAB + f[ai])
        except (ValueError, IndexError):
            die(f"{reader.path}: row {row} is malformed")
        stamp(f, hit)
        if not row % CHECKPOINT_ROWS:
            check()
    reader.close()
    return row - 1


def merge(reader: Reader, reopen, cons: Consequence, stamper: Stamper, regions: RegionFile, check) -> int:
    if not cons.by_chrom:
        walk(reader, lambda c, p: regions.add(c, p) or NOTHING, lambda f, hit: None, check)
        cons.open_regions(regions)
        reader = reopen()
    return walk(reader, cons.variants_at, stamper.stamp, check)


def lookup(reader: Reader, reopen, cons: Consequence, stamper: Stamper, regions: RegionFile, max_rows: int, check) -> int:
    ci, pi, ri, ai = reader.layout.indexes[:4]
    nsplit = reader.layout.nsplit
    table: dict = {}
    positions = array("q")
    row = 1
    for line in reader.lines:
        row += 1
        f = line.split(TAB, nsplit)
        try:
            key = TAB.join((f[ci], f[pi], f[ri], f[ai]))
            if key not in table:
                c, p = plain_integers(reader.path, row, f[ci], f[pi])
                table[key] = None
                positions.append(c << CHROM_SHIFT | p)
                if len(table) > max_rows:
                    die(f"{reader.path}: more than --max-lookup-rows {max_rows} distinct variants")
        except (ValueError, IndexError):
            die(f"{reader.path}: row {row} is malformed")
        if not row % CHECKPOINT_ROWS:
            check()
    reader.close()
    log(f"{row - 1} rows, {len(table)} distinct variants")

    positions = array("q", sorted(positions))
    if cons.by_chrom:
        for chrom in sorted({k >> CHROM_SHIFT for k in positions}):
            cons.open_chrom(chrom)
            cons.fill(table)
    else:
        for k in positions:
            regions.add(k >> CHROM_SHIFT, k & POS_MASK)
        cons.open_regions(regions)
        cons.fill(table)
    del positions
    check()

    reader = reopen()
    for line in reader.lines:
        f = line.split(TAB, nsplit)
        try:
            hit = table[TAB.join((f[ci], f[pi], f[ri], f[ai]))]
        except (KeyError, IndexError):
            die(f"{reader.path}: changed between the two passes")
        stamper.stamp(f, hit)
    reader.close()
    return row - 1


def annotate(args) -> dict:
    if os.path.exists(args.output) and os.path.samefile(args.input, args.output):
        die("--output is the input")
    check = rss_guard(int(args.max_rss_gb * (1 << 30)))

    def reopen() -> Reader:
        return Reader(args.input, args.threads)

    reader = reopen()
    index = Index(args.input) if reader.bgzf else None
    writer = Writer(args.output, reader.bgzf, args.threads)
    stamper = Stamper(reader, writer.write)
    writer.write(reader.header)
    if args.clear:
        rows = 0
        for line in reader.lines:
            rows += 1
            stamper.stamp(line.split(TAB, reader.layout.nsplit), None)
            if not rows % CHECKPOINT_ROWS:
                check()
        reader.close()
    else:
        access = args.consequence_access
        if access == "auto":
            small = os.path.getsize(args.input) <= args.regions_max_bytes
            access = "regions" if args.mode == "lookup" or small else "scan"
        log(f"{args.mode}, consequence read by {access}")
        cons = Consequence(args.consequence, access == "scan")
        regions = RegionFile(args.output + ".regions.bed")
        if args.mode == "lookup":
            rows = lookup(reader, reopen, cons, stamper, regions, args.max_lookup_rows, check)
        else:
            rows = merge(reader, reopen, cons, stamper, regions, check)
        cons.close()
    if rows != stamper.rows:
        die(f"{args.input}: {rows} rows read, {stamper.rows} written")
    writer.close()
    writer.publish(index, args.threads)
    return {"mode": "clear" if args.clear else args.mode, "rows": stamper.rows, "matched": stamper.matched,
            **transition_counts(stamper.transitions)}


def verify(original: str, stamped: str, threads: int) -> dict:
    a, b = Reader(original, threads), Reader(stamped, threads)
    mi, gi = a.layout.indexes[4:]
    transitions = defaultdict(int)
    rows_a = rows_b = differing = 0
    first = None
    for la, lb in zip_longest(a.lines, b.lines):
        rows_a += la is not None
        rows_b += lb is not None
        if la is None or lb is None:
            continue
        fa, fb = la.rstrip(b"\n").split(TAB), lb.rstrip(b"\n").split(TAB)
        same = len(fa) == len(fb) > max(mi, gi) and la.endswith(b"\n") == lb.endswith(b"\n")
        if same:
            transitions[fa[mi], fa[gi], fb[mi], fb[gi]] += 1
            fa[mi], fa[gi] = fb[mi], fb[gi]
            same = fa == fb
        if not same:
            differing += 1
            first = first or rows_a + 1
    a.close()
    b.close()
    return {
        "mode": "verify", "rows_original": rows_a, "rows_stamped": rows_b,
        "header_identical": a.header == b.header,
        "rows_differing_outside_annotation": differing, "first_differing_row": first,
        **transition_counts(transitions),
        "ok": a.header == b.header and rows_a == rows_b and not differing,
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", help="served file to stamp, bgzip or plain")
    p.add_argument("--output", help="stamped file; a bgzip one gets its index next to it")
    p.add_argument("--consequence", help="bgzip + tabix consequence file")
    p.add_argument("--mode", choices=["merge", "lookup"], default="merge",
                   help="merge for a variant-sorted file, lookup for a gene-indexed one")
    p.add_argument("--clear", action="store_true", help="write NA to both columns, read no consequence file")
    p.add_argument("--verify", nargs=2, metavar=("ORIGINAL", "STAMPED"), help="compare two files instead of stamping")
    p.add_argument("--consequence-access", choices=["auto", "regions", "scan"], default="auto")
    p.add_argument("--regions-max-bytes", type=int, default=32 << 20,
                   help="merge mode reads the consequence file by regions for an input up to this size on disk")
    p.add_argument("--max-lookup-rows", type=int, default=10_000_000,
                   help="lookup mode stops when the input holds more distinct variants than this")
    p.add_argument("--max-rss-gb", type=float, default=3.0, help="stop when resident memory passes this")
    p.add_argument("--threads", type=int, default=min(4, os.cpu_count()), help="bgzip and tabix threads")
    args = p.parse_args()
    if args.verify:
        if args.input or args.output or args.clear:
            p.error("--verify takes only the two files")
    elif not (args.input and args.output and (args.clear or args.consequence)):
        p.error("give --input, --output and one of --consequence, --clear")

    max_rss = int(args.max_rss_gb * (1 << 30))
    # backstop for growth between two checkpoints
    resource.setrlimit(resource.RLIMIT_AS, (4 * max_rss, 4 * max_rss))
    t0 = time.time()
    try:
        summary = verify(*args.verify, args.threads) if args.verify else annotate(args)
    finally:
        for path in SCRATCH:
            if os.path.exists(path):
                os.unlink(path)
    summary["wall_seconds"] = round(time.time() - t0, 2)
    summary["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss >> 10
    print(json.dumps(summary))
    if summary.get("ok") is False:
        sys.exit(1)


if __name__ == "__main__":
    main()
