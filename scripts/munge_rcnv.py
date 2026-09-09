#!/usr/bin/env python3
"""Munge the Collins et al. 2022 rare-CNV dosage sensitivity map into suite tables.

Source (Zenodo record 6347673, v0.2 2022-03-11, CC-BY 4.0):
  Collins et al., "A cross-disorder dosage sensitivity map of the human genome",
  Cell 2022, 185(16):3041-3055, doi:10.1016/j.cell.2022.06.036.

  --product scores -> Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz
      18,641 rows, header `#gene pHaplo pTriplo`, one row per autosomal protein-coding
      gene of Gencode v19. No coordinates in the file at all.

`--product` exists because the same Zenodo record also ships the gene association
sumstats, the sliding-window sumstats and the gene feature matrix. Only `scores` is
implemented; the others are separate work and are not stubbed here.

WHAT WAS READ OFF THE BYTES (not taken from the paper):
  - 18,641 rows, no duplicate gene symbols, no header comment block beyond line 1.
  - pHaplo/pTriplo are plain probabilities in [0,1], full double precision, no NA token.
  - pHaplo >= 0.86 selects 2,987 genes and pTriplo >= 0.94 selects 1,559, matching the
    paper's abstract exactly. Those are therefore the published thresholds and they are
    written into the output as boolean columns so no consumer has to re-derive them.
  - The gene symbols are Gencode v19 (GRCh37-era). All 18,641 are present in the
    `gene_name_19` column of the gencode mapping file, so the symbol -> ENSG step never
    misses; what misses is the ENSG -> *current* symbol step, for genes Gencode later
    dropped or left unnamed.

NO COORDINATES ARE ADDED. The scores are a per-gene property, the source file carries no
positions, and adding GRCh37 ones (or lifting them) would create a build-dependent column
where the data has none. The table is keyed by gene symbol and is build-independent by
construction, so it joins to the suite's GRCh38 gene views by symbol / ENSG.

SYMBOL RESOLUTION (--gencode-mapping, --hgnc):
  1. Gencode: `gene_name_19` -> `ensg` in gencode_gene_name_mapping_49-45-43-39-35-32-19.tsv
     (built by scripts/create_gene_name_mapping_across_gencode_versions.py), then the newest
     non-NA `gene_name_*` for that ENSG is the current symbol. Gencode writes the ENSG id
     itself as `gene_name` for genes it carries but does not name; those are not symbols and
     are skipped when picking the current name.
  2. Where two ENSGs share one v19 symbol, the candidates are ordered: one whose current
     symbol IS the v19 symbol beats one Gencode renamed, which beats a clone-style
     placeholder (RP11-*, AC0*.1, ...), which beats one Gencode leaves unnamed; the ENSG id
     only breaks what is still tied, so a rerun produces the same file. Ordering on the ENSG
     id alone would hand a gene its own name's clone id -- TUBB3 has two, and the clone
     locus has the smaller ENSG.
  3. HGNC fallback for the ENSGs Gencode no longer names: the v19 symbol is looked up
     (case-insensitively -- HGNC spells the C*orf* class in mixed case) among approved
     symbols, then `prev_symbol`, then `alias_symbol` of the HGNC complete set -- but never
     onto a symbol another ENSG already holds, because HGNC records a merged locus as a
     `prev_symbol` and following that link makes the join key ambiguous (see
     resolve_fallback).
  4. Anything neither resolves keeps its v19 symbol, with `ensembl_gene_id` still filled from
     Gencode. Most are clone-based placeholder names (RP11-*, AC0*.1, CTA-*) that never had
     an approved symbol, so there is nothing to update them to.
  The counts of each route are printed on every run, never assumed. What is asserted is
  what an assertion can actually catch: the input file's shape (row count, threshold
  counts) and the output's key integrity (ensembl_gene_id unique and never NA, and a
  ceiling on duplicate symbols -- the one check a mapping regression cannot slip past).
  See docs/rcnv-dosage-sensitivity.md.

`symbol_gencode_v19` is kept alongside `symbol` because the sibling rCNV products (gene
association sumstats, sliding windows) are keyed on the v19 symbol, and because a published
result must stay traceable to the identifier the paper actually used.

Output (one bgzipped TSV, no tabix index -- there is nothing to index):
  symbol  symbol_gencode_v19  ensembl_gene_id  phaplo  ptriplo  haploinsufficient  triplosensitive

Usage:
  python3 scripts/munge_rcnv.py --product scores --download --output out/collins_rcnv_2022_dosage_sensitivity.tsv.gz
  scripts/munge_rcnv.sh            # produce locally
  scripts/munge_rcnv.sh --stage    # produce and publish to both buckets
"""

import argparse
import csv
import gzip
import io
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

ZENODO_RECORD = "6347673"
SCORES_FILE = "Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz"
SCORES_URL = f"https://zenodo.org/records/{ZENODO_RECORD}/files/{SCORES_FILE}?download=1"

# the mapping inputs are already staged in the daly bucket; --hgnc-url is the public
# alternative for a machine without access to it
MAPPING_BUCKET = "gs://daly-genetics-results/mapping_files"
GENCODE_MAPPING_FILE = "gencode_gene_name_mapping_49-45-43-39-35-32-19.tsv"
HGNC_FILE = "hgnc_complete_set.txt"
HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
)

# descending, and must stay a subset of the versions in GENCODE_MAPPING_FILE's name --
# the newest column that names the gene is the current symbol
GENCODE_VERSIONS = ["49", "45", "43", "39", "35", "32"]

# published thresholds; they reproduce the abstract's 2,987 / 1,559 exactly (see EXPECTED)
PHAPLO_THRESHOLD = 0.86
PTRIPLO_THRESHOLD = 0.94

# input-file integrity only: these three are invariant to every symbol-mapping decision
# below, so they catch a truncated or re-released Zenodo file and nothing else
EXPECTED_ROWS = 18641
EXPECTED_HAPLOINSUFFICIENT = 2987
EXPECTED_TRIPLOSENSITIVE = 1559

# the mapping's own ceiling, taken from a run: symbols carrying more than one row. Each of
# the 10 pairs is one ENSG Gencode renamed onto this symbol and one ENSG Gencode no longer
# names at all (gene_name_49 is NA), which falls through to the unmapped route and keeps its
# v19 spelling -- the same string the other row was renamed to. The unmapped route never
# consults `claimed`, because inventing a name for a gene Gencode itself declines to name
# would be worse than a duplicate. Raising this number wants the printed collision list
# checked against the doc before it moves.
MAX_DUPLICATE_SYMBOLS = 10

COLUMNS = [
    "symbol",
    "symbol_gencode_v19",
    "ensembl_gene_id",
    "phaplo",
    "ptriplo",
    "haploinsufficient",
    "triplosensitive",
]

RESOURCE = "rcnv"
DATASET = "collins_rcnv_2022"
GCS = {
    "scores": (
        f"gs://finngen-commons/results_api_data/{RESOURCE}/{DATASET}/"
        f"{DATASET}_dosage_sensitivity.tsv.gz",
        f"gs://daly-genetics-results/{RESOURCE}/{DATASET}/"
        f"{DATASET}_dosage_sensitivity.tsv.gz",
    ),
}

# Gencode writes the bare ENSG id as gene_name for genes it carries but does not name
_ENSG_AS_NAME = re.compile(r"^ENSG\d+$")
# clone/contig-derived placeholder names: real strings, but never what a gene is called
_CLONE_NAME = re.compile(r"^(RP\d+-|AC\d|AL\d|AP\d|CT[ABD]-|Z\d)")
NA = "NA"


def fetch(url: str, dest: Path) -> Path:
    """Download `url` to `dest` unless it is already cached there."""
    if dest.exists():
        print(f"  cached: {dest}", file=sys.stderr)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  downloading {url}", file=sys.stderr)
    tmp = dest.with_suffix(dest.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "genetics-results-munge"})
    try:
        with urllib.request.urlopen(request, timeout=120) as response, tmp.open("wb") as out:
            shutil.copyfileobj(response, out)
    except (urllib.error.URLError, TimeoutError) as err:
        tmp.unlink(missing_ok=True)
        raise SystemExit(
            f"could not fetch {url}: {err}\n"
            f"download it by hand and put it at {dest}, then rerun without --download"
        ) from err
    tmp.rename(dest)
    return dest


def fetch_gcs(gcs_path: str, dest: Path) -> Path:
    if dest.exists():
        print(f"  cached: {dest}", file=sys.stderr)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  copying {gcs_path}", file=sys.stderr)
    subprocess.run(["gcloud", "storage", "cp", gcs_path, str(dest)], check=True)
    return dest


def read_scores(path: Path) -> list[tuple[str, str, str]]:
    """Read the Zenodo scores TSV as (gene, pHaplo, pTriplo) with the values verbatim.

    The probabilities are passed through as the source wrote them: rounding them would
    put the file's own boolean columns at risk of disagreeing with a consumer that
    recomputes the threshold from the rounded number.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if header[:3] != ["#gene", "pHaplo", "pTriplo"]:
            raise SystemExit(f"{path}: unexpected header {header}, expected #gene pHaplo pTriplo")
        rows = [(r[0], r[1], r[2]) for r in reader if r]
    genes = {r[0] for r in rows}
    if len(genes) != len(rows):
        raise SystemExit(f"{path}: {len(rows) - len(genes)} duplicate gene symbols")
    return rows


def read_gencode_mapping(path: Path) -> dict[str, list[tuple[str, str | None]]]:
    """{gene_name_19: [(ensg, current_symbol_or_None), ...]} from the mapping TSV."""
    by_v19: dict[str, list[tuple[str, str | None]]] = defaultdict(list)
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for version in GENCODE_VERSIONS + ["19"]:
            if f"gene_name_{version}" not in (reader.fieldnames or []):
                raise SystemExit(
                    f"{path} has no gene_name_{version} column; got {reader.fieldnames}"
                )
        for row in reader:
            v19 = row["gene_name_19"]
            if v19 == NA:
                continue
            current = next(
                (
                    row[f"gene_name_{v}"]
                    for v in GENCODE_VERSIONS
                    if row[f"gene_name_{v}"] != NA and not _ENSG_AS_NAME.match(row[f"gene_name_{v}"])
                ),
                None,
            )
            by_v19[v19].append((row["ensg"], current))
    if not by_v19:
        raise SystemExit(f"{path} carried no gene_name_19 entries")
    return by_v19


def read_hgnc(path: Path) -> tuple[dict[str, str], dict[str, list[str]], dict[str, list[str]]]:
    """(UPPER -> approved symbol, UPPER prev_symbol -> [approved], UPPER alias -> [approved]).

    Keyed in upper case because the two sources disagree on capitalisation of the same
    name: Gencode v19 writes C2ORF15, HGNC records it as C2orf15, and an exact-case
    lookup silently misses the record instead of resolving it. The value keeps HGNC's own
    spelling, which is what goes into the output.
    """
    approved: dict[str, str] = {}
    prev: dict[str, list[str]] = defaultdict(list)
    alias: dict[str, list[str]] = defaultdict(list)
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("status") != "Approved":
                continue
            symbol = row["symbol"]
            approved.setdefault(symbol.upper(), symbol)
            for field, index in (("prev_symbol", prev), ("alias_symbol", alias)):
                for name in (row.get(field) or "").split("|"):
                    if name and name.upper() != symbol.upper():
                        index[name.upper()].append(symbol)
    if not approved:
        raise SystemExit(f"{path} carried no approved symbols")
    return approved, prev, alias


def pick_gencode(
    v19: str, by_v19: dict[str, list[tuple[str, str | None]]]
) -> tuple[str, str | None]:
    """(ensg, current_symbol_or_None) for one Gencode v19 symbol.

    Where two ENSGs share the v19 symbol the locus that kept the name wins, then one
    Gencode renamed, then a clone-style placeholder, then one Gencode leaves unnamed; the
    ENSG id only breaks what is still tied. Preferring the smaller ENSG before any of that
    hands a gene the clone id sitting on its symbol -- TUBB3 -> AC092143.1 was six such.
    """
    candidates = by_v19.get(v19, [])
    if not candidates:
        return NA, None

    def preference(candidate: tuple[str, str | None]) -> tuple[int, str]:
        ensg, current = candidate
        if not current:
            rank = 3
        elif current == v19:
            rank = 0
        elif _CLONE_NAME.match(current):
            rank = 2
        else:
            rank = 1
        return rank, ensg

    return min(candidates, key=preference)


def hgnc_candidates(
    v19: str,
    approved: dict[str, str],
    prev: dict[str, list[str]],
    alias: dict[str, list[str]],
) -> list[tuple[str, str]]:
    """[(symbol, route)] HGNC offers for one v19 symbol, best route first."""
    key = v19.upper()
    candidates: list[tuple[str, str]] = []
    if key in approved:
        candidates.append((approved[key], "hgnc_approved"))
    for index, route in ((prev, "hgnc_prev"), (alias, "hgnc_alias")):
        candidates.extend((symbol, route) for symbol in sorted(index.get(key, [])))
    return candidates


def resolve_fallback(
    v19: str,
    approved: dict[str, str],
    prev: dict[str, list[str]],
    alias: dict[str, list[str]],
    claimed: dict[str, set[str]],
) -> tuple[str, str]:
    """(current_symbol, route) for a gene Gencode no longer names.

    `claimed` is {symbol: {ensg}} for every symbol already settled. HGNC records a merged or
    withdrawn locus by listing its old symbol as a `prev_symbol` of the surviving gene, so
    following that link renames a gene onto a symbol another ENSG in this same file already
    holds -- C2orf48 -> RRM2, C16orf47 -> ZFHX3, C17orf47 -> SEPTIN4. Both rows keep their
    own pHaplo/pTriplo, so the merge would make the join key ambiguous rather than current.
    Those keep the v19 symbol instead; `ensembl_gene_id` still identifies the gene exactly.

    The guard applies to every route, the approved one included: a symbol another gene keeps
    unchanged is just as taken as one Gencode assigned, and skipping the check there is what
    let MRC1L1 -> MRC1 land on top of MRC1's own row.
    """
    for symbol, route in hgnc_candidates(v19, approved, prev, alias):
        if symbol in claimed:
            continue
        return symbol, route
    return v19, "unmapped"


def munge_scores(
    scores_path: Path, gencode_path: Path, hgnc_path: Path, output: Path
) -> list[list[str]]:
    scores = read_scores(scores_path)
    by_v19 = read_gencode_mapping(gencode_path)
    approved, prev, alias = read_hgnc(hgnc_path)

    # symbols that keep their name are settled before any renaming route runs, so the
    # fallback can see every symbol already taken and by which ENSG. Both non-renaming
    # routes belong in that pass: Gencode's own current name, and the v19 symbol HGNC still
    # approves. Leaving the latter out is what let a prev_symbol rename land on a symbol
    # its own gene was keeping.
    picks = {v19: pick_gencode(v19, by_v19) for v19, _, _ in scores}
    claimed: dict[str, set[str]] = defaultdict(set)
    settled: dict[str, tuple[str, str]] = {}
    for v19, (ensg, current) in picks.items():
        if current:
            settled[v19] = (current, "gencode")
            claimed[current].add(ensg)
    # unlike the renaming pass below, this loop runs in file order, not sorted -- two v19
    # symbols differing only in case, both left unnamed by Gencode, resolving to the same
    # approved symbol would make the winner depend on Zenodo's row order; no such pair
    # exists in this input.
    for v19, (ensg, current) in picks.items():
        if current:
            continue
        approved_symbol = approved.get(v19.upper())
        if approved_symbol and approved_symbol not in claimed:
            settled[v19] = (approved_symbol, "hgnc_approved")
            claimed[approved_symbol].add(ensg)

    # the renaming routes run last and claim as they go, so two genes HGNC lists as previous
    # symbols of the same survivor cannot both take it (GATSL1 and GATSL2 -> CASTOR2). Sorted
    # rather than file order so which one wins does not depend on how Zenodo sorted the file.
    for v19 in sorted(picks):
        if v19 not in settled:
            symbol, route = resolve_fallback(v19, approved, prev, alias, claimed)
            settled[v19] = (symbol, route)
            if route != "unmapped":
                claimed[symbol].add(picks[v19][0])

    rows: list[list[str]] = []
    routes: Counter[str] = Counter()
    unmapped: list[str] = []
    for v19, phaplo, ptriplo in scores:
        ensg, _ = picks[v19]
        symbol, route = settled[v19]
        routes[route] += 1
        if route == "unmapped":
            unmapped.append(v19)
        rows.append(
            [
                symbol,
                v19,
                ensg,
                phaplo,
                ptriplo,
                "true" if float(phaplo) >= PHAPLO_THRESHOLD else "false",
                "true" if float(ptriplo) >= PTRIPLO_THRESHOLD else "false",
            ]
        )
    rows.sort(key=lambda r: (r[0], r[1]))

    haplo = sum(1 for r in rows if r[5] == "true")
    triplo = sum(1 for r in rows if r[6] == "true")
    report(rows, routes, unmapped, haplo, triplo)

    # these three check the input file only -- they are invariant to every mapping decision
    # above, so a symbol landing on the wrong locus passes all of them
    if len(rows) != EXPECTED_ROWS:
        raise SystemExit(f"expected {EXPECTED_ROWS} rows, got {len(rows)}")
    if haplo != EXPECTED_HAPLOINSUFFICIENT:
        raise SystemExit(
            f"expected {EXPECTED_HAPLOINSUFFICIENT} haploinsufficient genes "
            f"(pHaplo >= {PHAPLO_THRESHOLD}), got {haplo}"
        )
    if triplo != EXPECTED_TRIPLOSENSITIVE:
        raise SystemExit(
            f"expected {EXPECTED_TRIPLOSENSITIVE} triplosensitive genes "
            f"(pTriplo >= {PTRIPLO_THRESHOLD}), got {triplo}"
        )

    # what the mapping itself is checked on. ensembl_gene_id is the primary key, so it has
    # to be present and unique; the duplicate-symbol ceiling is the tell for a fallback that
    # renames a gene onto a symbol another one already holds.
    ensgs = [r[2] for r in rows]
    missing = sum(1 for e in ensgs if e == NA)
    if missing:
        raise SystemExit(
            f"{missing} rows have no ensembl_gene_id; every v19 symbol is supposed to be in "
            f"the mapping file's gene_name_19 column, so this means a mapping file that no "
            f"longer covers Gencode v19"
        )
    if len(set(ensgs)) != len(ensgs):
        repeated = sorted(e for e, n in Counter(ensgs).items() if n > 1)
        raise SystemExit(
            f"ensembl_gene_id is not unique ({len(repeated)} repeated, e.g. {repeated[:5]}); "
            f"it is this table's primary key, so two v19 symbols resolving to one ENSG is a "
            f"pick_gencode regression"
        )
    duplicate_symbols = sorted(s for s, n in Counter(r[0] for r in rows).items() if n > 1)
    if len(duplicate_symbols) > MAX_DUPLICATE_SYMBOLS:
        raise SystemExit(
            f"{len(duplicate_symbols)} symbols carry more than one row, above the "
            f"{MAX_DUPLICATE_SYMBOLS} this mapping is known to produce: "
            f"{duplicate_symbols}. That means either a rename landed on a symbol a "
            f"still-named gene already holds, or a fallback route bypassed the `claimed` "
            f"guard -- check the printed list against the doc's known 10 before raising "
            f"MAX_DUPLICATE_SYMBOLS"
        )

    write_bgzip(rows, output)
    return rows


def report(
    rows: list[list[str]],
    routes: Counter[str],
    unmapped: list[str],
    haplo: int,
    triplo: int,
) -> None:
    hgnc = sum(v for k, v in routes.items() if k.startswith("hgnc_"))
    print("", file=sys.stderr)
    print(f"rows:                          {len(rows)}", file=sys.stderr)
    print(f"symbol via gencode:            {routes['gencode']}", file=sys.stderr)
    print(f"symbol via HGNC:               {hgnc}", file=sys.stderr)
    for route in ("hgnc_approved", "hgnc_prev", "hgnc_alias"):
        print(f"  {route[5:]:<28} {routes[route]}", file=sys.stderr)
    print(f"unmapped (kept as v19 symbol): {routes['unmapped']}", file=sys.stderr)
    print(f"symbol changed from v19:       {sum(1 for r in rows if r[0] != r[1])}", file=sys.stderr)
    print(f"ensembl_gene_id missing:       {sum(1 for r in rows if r[2] == NA)}", file=sys.stderr)
    collisions = [s for s, n in Counter(r[0] for r in rows).items() if n > 1]
    print(f"symbols carrying >1 row:       {len(collisions)}", file=sys.stderr)
    if collisions:
        print(f"  {sorted(collisions)}", file=sys.stderr)
    print(f"haploinsufficient (>= {PHAPLO_THRESHOLD}):    {haplo}", file=sys.stderr)
    print(f"triplosensitive   (>= {PTRIPLO_THRESHOLD}):    {triplo}", file=sys.stderr)
    print(f"unmapped symbols: {sorted(unmapped)}", file=sys.stderr)


def write_bgzip(rows: list[list[str]], output: Path) -> None:
    """One bgzipped TSV, no index.

    The output family writers in sumstat_utils/peak_utils all build a tabix index over
    coordinates. This product has none, so it writes its own bgzip pipe: the file is a
    BigQuery load file only.
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as fh:
        proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=fh)
        with io.TextIOWrapper(proc.stdin, "utf-8", newline="") as pipe:
            writer = csv.writer(pipe, delimiter="\t", lineterminator="\n")
            writer.writerow(COLUMNS)
            writer.writerows(rows)
        if proc.wait() != 0:
            raise SystemExit("bgzip failed")
    print(f"\nwrote {output}", file=sys.stderr)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--product",
        choices=["scores"],
        required=True,
        help="which Zenodo product to munge; only the dosage-sensitivity scores today",
    )
    parser.add_argument("--scores", type=Path, help=f"local {SCORES_FILE} (default: <cache-dir>/{SCORES_FILE})")
    parser.add_argument(
        "--gencode-mapping",
        type=Path,
        help=f"gencode gene name mapping TSV (default: <cache-dir>/{GENCODE_MAPPING_FILE})",
    )
    parser.add_argument("--hgnc", type=Path, help=f"HGNC complete set TSV (default: <cache-dir>/{HGNC_FILE})")
    parser.add_argument("--cache-dir", type=Path, default=Path("data/rcnv"), help="where downloaded inputs are cached")
    parser.add_argument(
        "--download",
        action="store_true",
        help="fetch any missing input into --cache-dir (Zenodo for the scores, the daly "
        "mapping_files bucket for the gencode mapping and HGNC set); OFF by default",
    )
    parser.add_argument("--hgnc-url", default=HGNC_URL, help="public HGNC source, used with --download when gcloud is unavailable")
    parser.add_argument("--output", type=Path, help="output .tsv.gz (default derived from --product)")
    parser.add_argument("--stage", action="store_true", help="upload the output to both profile buckets; OFF by default")
    parser.add_argument("--gcs-finngen", help="override the finngen destination")
    parser.add_argument("--gcs-daly", help="override the daly destination")
    return parser.parse_args()


def resolve_inputs(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    scores = args.scores or args.cache_dir / SCORES_FILE
    gencode = args.gencode_mapping or args.cache_dir / GENCODE_MAPPING_FILE
    hgnc = args.hgnc or args.cache_dir / HGNC_FILE
    if args.download:
        fetch(SCORES_URL, scores)
        fetch_gcs(f"{MAPPING_BUCKET}/{GENCODE_MAPPING_FILE}", gencode)
        if not hgnc.exists():
            try:
                fetch_gcs(f"{MAPPING_BUCKET}/{HGNC_FILE}", hgnc)
            except (subprocess.CalledProcessError, FileNotFoundError):
                print("  bucket copy failed, falling back to genenames.org", file=sys.stderr)
                fetch(args.hgnc_url, hgnc)
    for path, flag in ((scores, "--scores"), (gencode, "--gencode-mapping"), (hgnc, "--hgnc")):
        if not path.exists():
            raise SystemExit(f"{path} not found (pass {flag}, or run with --download)")
    return scores, gencode, hgnc


def main() -> None:
    args = parse_args()
    if shutil.which("bgzip") is None:
        raise SystemExit("bgzip not found on PATH (run inside the munge image)")

    scores, gencode, hgnc = resolve_inputs(args)
    output = args.output or Path(f"{DATASET}_dosage_sensitivity.tsv.gz")
    munge_scores(scores, gencode, hgnc, output)

    if args.stage:
        finngen, daly = GCS[args.product]
        print("Staging to GCS (both buckets) ...", file=sys.stderr)
        # the two buckets live in different projects and a host commonly reaches only one of
        # them, so a failure on the first must not decide whether the second is even tried
        failures = []
        for dest in (args.gcs_finngen or finngen, args.gcs_daly or daly):
            result = subprocess.run(["gcloud", "storage", "cp", str(output), dest])
            if result.returncode == 0:
                print(f"  uploaded {dest}", file=sys.stderr)
            else:
                failures.append(dest)
                print(f"  FAILED (exit {result.returncode}) {dest}", file=sys.stderr)
        if failures:
            raise SystemExit(f"{len(failures)} of 2 destinations failed: {', '.join(failures)}")
    else:
        print("  --stage not set: skipping GCS upload", file=sys.stderr)


if __name__ == "__main__":
    main()
