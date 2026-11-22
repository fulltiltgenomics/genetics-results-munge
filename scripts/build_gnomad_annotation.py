#!/usr/bin/env python3
"""Deduplicate the gnomAD genomes+exomes sites table to one row per variant, and cut the
consequence file that credible-set files are stream-merged against.

Source file
-----------
A position-sorted gnomAD sites TSV (bgzip) holding genome and exome rows, tab separated with
header

    #chr pos ref alt rsids filters AN AF AF_* most_severe gene_most_severe consequences genome_or_exome

`--gnomad` takes that pre-merged file. `--genomes` and `--exomes` instead take one
position-sorted file per source and merge them by position in-process; `genome_or_exome` is
appended from the flag when a file lacks it and is believed as written when a file has it.
Inputs are local paths: nothing here streams from a bucket.

Columns are located by header name. The four key columns must lead the row and
`genome_or_exome` must end it, because the hot loop splits only the leading columns and
reads the source off the tail; any other layout is refused rather than guessed at.

Outputs
-------
Both bgzip + tabix (`-s1 -b2 -e2`), sorted by numeric chromosome, position, ref, alt:

- `--output`: the served sites file, every input column kept, one row per variant;
- `--consequence-output`: `#chr pos ref alt most_severe gene_most_severe`.

In both, `most_severe` and `gene_most_severe` are recomputed from the row's `consequences`
column (next section); the input's own two values survive only on a row with no transcript.

Consequence and gene
--------------------
The input's `most_severe` is VEP's top-level term, the worst over every transcript the
variant touches, so a frameshift in one nonsense-mediated-decay isoform outranks the intron
the gene's canonical transcript sees (15-90885291-CT-C in v4.0). It is replaced by the
canonical-transcript consequence, which is what the gnomAD browser leads with, following
gnomad_methods' `process_consequences` (`worst_csq_for_variant_canonical`) with the
departures listed below. `consequences` is a JSON array with one object per distinct
transcript annotation: `gene_symbol gene_id consequences gene_symbol_source canonical
biotype`.

1. each entry is ranked by the most severe of its terms in `CSQ_ORDER`;
2. the entries are narrowed to the first of these that is not empty: canonical
   (`canonical == 1`, VEP's flag) and protein-coding; canonical; protein-coding; all;
3. the lowest rank in that set gives `most_severe`, and its gene `gene_most_severe`;
4. several genes at that rank: a symbol that is not `Clone_based_ensembl_gene` first, then
   an entry that has a symbol, then the lowest `gene_id`;
5. the gene is written as its symbol, or as its `gene_id` when VEP gives no symbol;
6. a row whose `consequences` is `NA` or `[]` (intergenic, or regulatory only) keeps the
   input's `most_severe` and `gene_most_severe`.

Departures from `process_consequences`: protein-coding genes come before severity (it ranks
a lncRNA's intron above a protein-coding gene's upstream variant); a variant that reaches no
canonical transcript falls back to the others instead of going without a consequence; ties
are broken (it leaves them to sort order); and LOFTEE's HC/LC/OS score adjustment is not
applied, because the export does not carry `lof`. With one canonical transcript per gene
that adjustment can only change which of two genes is named, and only where both carry a
loss-of-function term: 61 of 6.18M rows of the v4.1.1 vep115 file over 15:88-93M,
19:40-45M and 1:150-155M, all with the same term in both genes. OS on a term that is not
loss-of-function is the case that count cannot see.

A term missing from `CSQ_ORDER`, an entry with no terms and JSON that does not parse each
stop the run, naming the variant.

Which row survives
------------------
A variant seen in both sources keeps its exome row only when the exome `AN` is strictly
greater than the genome `AN`; a tie keeps the genome row. `AN` is the only statistic
compared. A non-integer `AN` on a variant that has both rows stops the run: there is no
ruling for it.

Chromosome coding
-----------------
The `chr` prefix is stripped case-insensitively and the letter chromosomes are written as
`sumstat_utils.CHROM_CODES` has them — X 23, Y 24, M/MT 26 — which is the one coding every
writer in this repo uses and the one results-api's tabix reader, which serves this file, is
being moved to. A contig still non-numeric after the mapping is dropped and counted.

Format assumptions (checked against the delivered v4.0 file)
-----------------------------------------------------------
- contigs are spelled 1-22, X, Y, in that order, with no mitochondrion;
- a variant present in both sources is two rows at one position, in either order and not
  necessarily adjacent when the position is multi-allelic;
- `canonical` is 1 or null, and no gene has two canonical entries on one row (v4.1.1 vep115,
  the three regions above). The v4.0 file also carries RefSeq entries with their own
  canonical flag; they share the Ensembl entry's symbol, so they add ties, not genes.

Contig order and the guards
---------------------------
Each contig's rows are compressed into their own part file and the parts are concatenated
in numeric order at the end, so an input that orders contigs differently from their numeric
codes still yields a sorted file without a global sort. Outputs from an earlier run are
removed at the start and the new ones appear only when the run succeeds; a failed run
leaves `*.part.*` files behind and nothing a reader could mistake for a result.

The run exits non-zero, naming the position, on: a position lower than the one before it
on the same contig, a contig that comes back after another one started, more than two rows
for one variant, and two rows from the same source for one variant. Only one position's
rows are held at a time; `--max-rss-gb` is checked as the stream advances and enforced a
second time, on address space, through RLIMIT_AS.

Usage
-----
    python scripts/build_gnomad_annotation.py \\
        --gnomad /mnt/disks/data/gnomad/gnomad.genomes.exomes.v4.0.sites.v2.tsv.bgz \\
        --output /mnt/disks/data/gnomad/out/gnomad.v4.0.sites.dedup.tsv.bgz \\
        --consequence-output /mnt/disks/data/gnomad/out/gnomad.v4.0.consequence.tsv.bgz
"""

import argparse
import fcntl
import heapq
import json
import os
import resource
import shutil
import subprocess
import sys
import time

from sumstat_utils import CHROM_CODES

TAB = b"\t"
KEY_COLUMNS = [b"#chr", b"pos", b"ref", b"alt"]
SOURCE_COLUMN = b"genome_or_exome"
CONSEQUENCE_HEADER = b"#chr\tpos\tref\talt\tmost_severe\tgene_most_severe\n"
CONSEQUENCE_ROW = b"%b\t%b\t%b\t%b\t%b\t%b\n"
# severity order, most severe first: gnomad_methods' gnomad.utils.vep.CSQ_ORDER, copied because
# importing it pulls in hail. A term VEP adds later is absent here and stops the run
CSQ_ORDER = (
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained",
    "frameshift_variant", "stop_lost", "start_lost", "initiator_codon_variant",
    "transcript_amplification", "inframe_insertion", "inframe_deletion", "missense_variant",
    "protein_altering_variant", "splice_donor_5th_base_variant", "splice_region_variant",
    "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant", "start_retained_variant", "stop_retained_variant",
    "synonymous_variant", "coding_sequence_variant", "coding_transcript_variant",
    "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant", "non_coding_exon_variant", "intron_variant",
    "NMD_transcript_variant", "non_coding_transcript_variant", "nc_transcript_variant",
    "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification",
    "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification",
    "feature_elongation", "regulatory_region_variant", "feature_truncation",
    "intergenic_variant", "sequence_variant",
)
CSQ_RANK = {term: rank for rank, term in enumerate(CSQ_ORDER)}
# which transcripts decided a row, in the order they are preferred
TIERS = ("canonical_protein_coding", "canonical_other", "noncanonical_protein_coding", "noncanonical_other")
NO_TRANSCRIPT = len(TIERS)
# `consequences` of a row with no transcript: hail writes a missing array as NA
NO_TRANSCRIPT_VALUES = (b"NA", b"[]")
# neighbouring variants in one gene repeat the same JSON; the cache is emptied at this many
# distinct values so that it stays far below --max-rss-gb
CACHE_ENTRIES = 1 << 16
# every bgzip stream ends with this empty block; all but the last are cut when parts are
# concatenated, because older BGZF readers stop at the first one they meet
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
CHECKPOINT_ROWS = 1 << 24
PAGE_SIZE = os.sysconf("SC_PAGE_SIZE")
# the source column's two values, as the byte a row ends in before its newline
GENOME, EXOME = ord("g"), ord("e")


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def die(msg: str) -> None:
    sys.exit(f"build_gnomad_annotation: {msg}")


def rss_bytes() -> int:
    with open("/proc/self/statm") as f:
        return int(f.read().split()[1]) * PAGE_SIZE


class Layout:
    """Column positions read out of a header line."""

    def __init__(self, header: bytes, path: str):
        cols = header.rstrip(b"\n").split(TAB)
        if cols[:4] != KEY_COLUMNS:
            die(f"{path}: header must start with {b' '.join(KEY_COLUMNS).decode()}")
        self.has_source = SOURCE_COLUMN in cols
        if self.has_source and cols[-1] != SOURCE_COLUMN:
            die(f"{path}: {SOURCE_COLUMN.decode()} must be the last column")
        try:
            self.an, self.most_severe, self.gene, self.consequences = (
                cols.index(c) for c in (b"AN", b"most_severe", b"gene_most_severe", b"consequences")
            )
        except ValueError as e:
            die(f"{path}: missing column, {e}")
        self.columns = cols if self.has_source else cols + [SOURCE_COLUMN]
        # a row is split only this far; the source column always follows, so the last piece
        # of the split is never one of these
        self.nsplit = max(self.an, self.most_severe, self.gene, self.consequences) + 1


class ChromCoder:
    """Contig name -> numeric code as bytes, or None for a contig that has none."""

    def __init__(self):
        self.special = {name.encode(): code.encode() for name, code in CHROM_CODES.items()}
        self.cache: dict[bytes, bytes | None] = {}

    def __call__(self, name: bytes) -> bytes | None:
        if name not in self.cache:
            base = name.upper().removeprefix(b"CHR")
            base = self.special.get(base, base)
            self.cache[name] = base if base.isdigit() else None
        return self.cache[name]


def pick_consequence(consequences: bytes):
    """(most_severe, gene, tier) decided by the transcripts in a row's `consequences`, or
    None for a row that has none. The rule is the module docstring's; raises KeyError on a
    term outside CSQ_ORDER and ValueError on an entry without terms or on bad JSON."""
    if consequences in NO_TRANSCRIPT_VALUES:
        return None
    best = best_key = None
    for entry in json.loads(consequences):
        tier = 2 * (entry["canonical"] != 1) + (entry["biotype"] != "protein_coding")
        key = (
            tier,
            min(CSQ_RANK[term] for term in entry["consequences"]),
            entry["gene_symbol_source"] == "Clone_based_ensembl_gene",
            entry["gene_symbol"] is None,
            entry["gene_id"] or "",
        )
        if best_key is None or key < best_key:
            best, best_key = entry, key
    if best is None:
        raise ValueError("empty consequences array")
    gene = best["gene_symbol"] or best["gene_id"] or "NA"
    return CSQ_ORDER[best_key[1]].encode(), gene.encode(), best_key[0]


def widen(pipe) -> None:
    """Grow a pipe to the largest size an unprivileged process may ask for: at the default
    every buffer-sized read or write is cut into many small ones, each a context switch."""
    try:
        with open("/proc/sys/fs/pipe-max-size") as f:
            fcntl.fcntl(pipe, fcntl.F_SETPIPE_SZ, int(f.read()))
    except OSError:
        pass


def open_input(path: str, threads: int):
    proc = subprocess.Popen(["bgzip", "-dc", f"-@{threads}", path], stdout=subprocess.PIPE, bufsize=1 << 20)
    widen(proc.stdout)
    header = proc.stdout.readline()
    return proc, Layout(header, path)


def keyed(lines, source: bytes | None, coder: ChromCoder, stats: dict):
    """(chromosome code, position, line) for merging two streams, with the source column
    appended when `source` is given. Uncoded contigs are dropped here rather than in the
    main loop, as there is no rank to merge them by."""
    suffix = None if source is None else TAB + source + b"\n"
    for line in lines:
        chrom, pos, _ = line.split(TAB, 2)
        code = coder(chrom)
        if code is None:
            stats["dropped_contig_rows"] += 1
            continue
        yield int(code), int(pos), line if suffix is None else line[:-1] + suffix


def open_sources(args, coder: ChromCoder, stats: dict):
    """-> (subprocesses, layout, iterator over data lines carrying the source column)"""
    if args.gnomad:
        proc, layout = open_input(args.gnomad, args.threads)
        if not layout.has_source:
            die(f"{args.gnomad}: a single input needs the {SOURCE_COLUMN.decode()} column")
        return [proc], layout, proc.stdout
    procs, layouts, streams = [], [], []
    for path, source in ((args.genomes, b"g"), (args.exomes, b"e")):
        proc, layout = open_input(path, max(1, args.threads // 2))
        procs.append(proc)
        layouts.append(layout)
        streams.append(keyed(proc.stdout, None if layout.has_source else source, coder, stats))
    if layouts[0].columns != layouts[1].columns:
        die("--genomes and --exomes have different columns")
    return procs, layouts[0], (row[2] for row in heapq.merge(*streams))


class PartWriter:
    """One bgzip stream per contig, concatenated in numeric contig order by `assemble`."""

    def __init__(self, path: str, threads: int, level: int):
        self.path = path
        self.threads = threads
        self.cmd = ["bgzip", "-c", f"-@{threads}", "-l", str(level)]
        self.parts: dict[int, str] = {}
        self.proc = None

    def start(self, code: bytes):
        """Close the running part and return the write method of a new one."""
        self.close()
        part = f"{self.path}.part.{int(code):02d}"
        self.parts[int(code)] = part
        with open(part, "wb") as fh:
            self.proc = subprocess.Popen(self.cmd, stdin=subprocess.PIPE, stdout=fh, bufsize=1 << 20)
        widen(self.proc.stdin)
        return self.proc.stdin.write

    def close(self) -> None:
        if self.proc is not None:
            self.proc.stdin.close()
            if self.proc.wait() != 0:
                die(f"bgzip failed writing a part of {self.path}")
            self.proc = None

    def assemble(self, header: bytes) -> None:
        self.close()
        head = subprocess.run(self.cmd, input=header, stdout=subprocess.PIPE, check=True).stdout
        with open(self.path + ".assembling", "wb") as out:
            out.write(head[: -len(BGZF_EOF)])
            for code in sorted(self.parts):
                part = self.parts[code]
                body = os.path.getsize(part) - len(BGZF_EOF)
                with open(part, "rb") as f:
                    while body > 0:
                        buf = f.read(min(body, 8 << 20))
                        out.write(buf)
                        body -= len(buf)
                    if f.read() != BGZF_EOF:
                        die(f"{part} does not end in a BGZF EOF block")
                os.remove(part)
            out.write(BGZF_EOF)
        os.replace(self.path + ".assembling", self.path)
        subprocess.run(["tabix", "-f", f"-@{self.threads}", "-s1", "-b2", "-e2", self.path], check=True)


def dedup(lines, layout: Layout, coder: ChromCoder, sites: PartWriter, cons: PartWriter, max_rss: int, stats: dict):
    """Stream the rows through, writing one row per variant; -> the key of the last one."""
    nsplit, i_an, i_ms, i_gene = layout.nsplit, layout.an, layout.most_severe, layout.gene
    i_cons = layout.consequences
    picked: dict[bytes, tuple | None] = {}
    tiers = [0] * (len(TIERS) + 1)
    term_changed = gene_changed = 0
    n_fields = nsplit + 1
    write_sites = write_cons = None
    cur_chr = cur_pos = code = out_code = None
    pos_int = 0
    rewrite = 0
    seen: set[bytes] = set()
    # the held position's rows, as its first row and a list of any others: most positions
    # have one row, and building a list per row is measurable at this row count
    first = first_line = None
    rest: list = []
    n = out_rows = pairs = exome_wins = 0
    checkpoint = CHECKPOINT_ROWS
    last = None

    def where(f) -> str:
        return b":".join(f[:4]).decode(errors="replace")

    def resolve(rows) -> list:
        """The surviving row of each variant at one position, in (ref, alt) order."""
        nonlocal pairs, exome_wins
        by_variant: dict = {}
        for row in rows:
            by_variant.setdefault((row[0][2], row[0][3]), []).append(row)
        winners = []
        for _, same in sorted(by_variant.items()):
            if len(same) > 2:
                die(f"{len(same)} rows for one variant at {where(same[0][0])}")
            if len(same) == 2:
                a, b = same
                if a[1][-2] == b[1][-2]:
                    die(f"two rows from the same source at {where(a[0])}")
                genome, exome = (a, b) if b[1][-2] == EXOME else (b, a)
                try:
                    exome_won = int(exome[0][i_an]) > int(genome[0][i_an])
                except ValueError:
                    die(f"non-integer AN on a variant with both rows at {where(a[0])}")
                pairs += 1
                exome_wins += exome_won
                same = [exome if exome_won else genome]
            winners.append(same[0])
        return winners

    def flush() -> None:
        nonlocal n, out_rows, checkpoint, last, first, rest, term_changed, gene_changed
        if rest:
            n += len(rest) + 1
            rows = resolve([(first, first_line), *rest])
            rest = []
        else:
            n += 1
            rows = ((first, first_line),)
        for f, line in rows:
            if line[-2] != EXOME and line[-2] != GENOME:
                die(f"{SOURCE_COLUMN.decode()} is neither g nor e at {where(f)}")
            consequences = f[i_cons]
            try:
                pick = picked[consequences]
            except KeyError:
                try:
                    pick = pick_consequence(consequences)
                except (KeyError, ValueError, TypeError) as e:
                    die(f"cannot read consequences at {where(f)}: {e!r}")
                if len(picked) >= CACHE_ENTRIES:
                    picked.clear()
                picked[consequences] = pick
            if pick is None:
                tiers[NO_TRANSCRIPT] += 1
            else:
                term, gene, tier = pick
                tiers[tier] += 1
                if term != f[i_ms] or gene != f[i_gene]:
                    if term != f[i_ms]:
                        term_changed += 1
                    else:
                        gene_changed += 1
                    f[i_ms], f[i_gene] = term, gene
                    line = TAB.join(f)
            write_sites(code + line[rewrite:] if rewrite else line)
            write_cons(CONSEQUENCE_ROW % (code, f[1], f[2], f[3], f[i_ms], f[i_gene]))
        out_rows += len(rows)
        last = f
        first = None
        if n >= checkpoint:
            checkpoint += CHECKPOINT_ROWS
            rss = rss_bytes()
            log(f"  {n:,} rows in, {out_rows:,} out, at {where(f)}, rss {rss >> 20} MB")
            if rss > max_rss:
                die(f"resident memory {rss >> 20} MB is over --max-rss-gb at {where(f)}")

    for line in lines:
        f = line.split(TAB, nsplit)
        if len(f) != n_fields:
            die(f"short row at {where(f)}")
        if f[1] == cur_pos and f[0] == cur_chr:
            rest.append((f, line))
            continue
        if first is not None:
            flush()
        if f[0] != cur_chr:
            cur_chr, cur_pos, pos_int = f[0], None, 0
            code = coder(cur_chr)
            if code is not None:
                if code in seen:
                    die(f"contig {code.decode()} resumes at {where(f)} after another contig started")
                seen.add(code)
                out_code = code
                rewrite = 0 if code == cur_chr else len(cur_chr)
                write_sites, write_cons = sites.start(code), cons.start(code)
                log(f"contig {cur_chr.decode()} -> {code.decode()}")
        if code is None:
            stats["dropped_contig_rows"] += 1
            continue
        p = int(f[1])
        if p < pos_int:
            die(f"position goes backwards at {where(f)}, after {pos_int}")
        cur_pos, pos_int = f[1], p
        first, first_line = f, line
    if first is not None:
        flush()
    stats.update(rows_in=n + stats["dropped_contig_rows"], rows_out=out_rows, duplicates_removed=n - out_rows,
                 variants_in_both=pairs, exome_wins=exome_wins, genome_wins=pairs - exome_wins,
                 most_severe_changed=term_changed, gene_only_changed=gene_changed,
                 **{f"decided_by_{name}": count for name, count in zip((*TIERS, "no_transcript"), tiers)})
    return last and (out_code, *last[1:4])


def self_check(path: str, key: tuple) -> None:
    """The last variant written must come back from the index as exactly one row."""
    code, pos, ref, alt = (k.decode() for k in key)
    r = subprocess.run(["tabix", path, f"{code}:{pos}-{pos}"], capture_output=True, text=True, check=True)
    hits = [l for l in r.stdout.splitlines() if l.split("\t")[:4] == [code, pos, ref, alt]]
    if len(hits) != 1:
        die(f"self-check: {path} returned {len(hits)} rows for {code}:{pos}:{ref}:{alt}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--gnomad", help="pre-merged genomes+exomes sites TSV (bgzip), position-sorted")
    p.add_argument("--genomes", help="genomes-only sites TSV (bgzip); give together with --exomes")
    p.add_argument("--exomes", help="exomes-only sites TSV (bgzip); give together with --genomes")
    p.add_argument("--output", required=True, help="deduplicated sites file (.bgz); the .tbi is written next to it")
    p.add_argument("--consequence-output", required=True, help="consequence file (.bgz); the .tbi is written next to it")
    p.add_argument("--threads", type=int, default=os.cpu_count(), help="bgzip threads, for decompression and for the sites output")
    p.add_argument("--compress-level", type=int, default=6, help="bgzip level of both outputs")
    p.add_argument("--max-rss-gb", type=float, default=2.0,
                   help="exit when this process's resident memory passes it; address space is capped at four times this")
    p.add_argument("--min-free-gb", type=float,
                   help="free space required in each output directory (default: 1.5x the input size)")
    args = p.parse_args()
    inputs = [args.gnomad] if args.gnomad else [args.genomes, args.exomes]
    if not all(inputs) or (args.gnomad and (args.genomes or args.exomes)):
        p.error("give either --gnomad, or both --genomes and --exomes")

    max_rss = int(args.max_rss_gb * (1 << 30))
    # address space, not resident memory, so it sits well above the ceiling: this is the
    # backstop for growth between two checkpoints, and bgzip's threads inherit it
    resource.setrlimit(resource.RLIMIT_AS, (4 * max_rss, 4 * max_rss))

    need = args.min_free_gb * (1 << 30) if args.min_free_gb is not None else 1.5 * sum(os.path.getsize(i) for i in inputs)
    for out in (args.output, args.consequence_output):
        out_dir = os.path.dirname(os.path.abspath(out))
        os.makedirs(out_dir, exist_ok=True)
        free = shutil.disk_usage(out_dir).free
        if free < need:
            die(f"{out_dir} has {free / (1 << 30):.1f} GB free, {need / (1 << 30):.1f} GB required")
        # an earlier result must not outlive a run that fails
        for stale in (out, out + ".tbi"):
            if os.path.exists(stale):
                os.remove(stale)

    stats = {"dropped_contig_rows": 0}
    coder = ChromCoder()
    procs, layout, lines = open_sources(args, coder, stats)
    sites = PartWriter(args.output, args.threads, args.compress_level)
    cons = PartWriter(args.consequence_output, 1, args.compress_level)
    last_key = dedup(lines, layout, coder, sites, cons, max_rss, stats)
    for proc in procs:
        if proc.wait() != 0:
            die("decompressing an input failed")
    if last_key is None:
        die("no rows on a coded contig in the input")
    log("concatenating parts and indexing")
    sites.assemble(TAB.join(layout.columns) + b"\n")
    cons.assemble(CONSEQUENCE_HEADER)
    self_check(args.output, last_key)
    self_check(args.consequence_output, last_key)
    for name, value in stats.items():
        log(f"{name}: {value:,}")
    log("done")


if __name__ == "__main__":
    main()
