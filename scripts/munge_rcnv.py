#!/usr/bin/env python3
"""Munge the Collins et al. 2022 rare-CNV dosage sensitivity map into suite tables.

Source (Zenodo record 6347673, v0.2 2022-03-11, CC-BY 4.0):
  Collins et al., "A cross-disorder dosage sensitivity map of the human genome",
  Cell 2022, 185(16):3041-3055, doi:10.1016/j.cell.2022.06.036.

  --product scores -> Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz
      18,641 rows, header `#gene pHaplo pTriplo`, one row per autosomal protein-coding
      gene of Gencode v19. No coordinates in the file at all.

  --product genes -> Collins_rCNV_2022.gene_association_sumstats.tar.gz
      108 tabixed BEDs (54 phenotypes x DEL/DUP), 17,263 rows each -- the same gene *set* in
      every file (row order is not identical: two genes tied on GRCh37 (chr, start) come out
      in a different relative order between a phenotype's DEL and DUP file -- confirmed on
      HP0000118), one row per autosomal protein-coding gene of Gencode v19. 21 source columns
      per the bundled README; phenotype code and CNV type come from the file name, not a
      column. GRCh37 chr/start/end are dropped: a row's coordinates come from
      `gene_annotations_v` at query time, the same way the `scores` product carries none.
      Overall 65.2% of rows are all-NA past `case_freq`/`control_freq` (mean 11,248 rows per
      file, ranging 591 to 16,570 depending on phenotype -- there is no single per-file
      figure). `meta_lnOR` onward is NA wherever the meta-analysis produced no estimate,
      which correlates with no CNV ever observed (`case_freq = control_freq = 0`) but is not
      implied by `n_nominal_cohorts`: 449,103 rows have `n_nominal_cohorts = 0` with a
      non-NA beta, and 17,412 NA rows have `n_nominal_cohorts >= 1`. Filtering on
      `n_nominal_cohorts` does not select the analysable rows -- filter on `beta IS NOT
      NULL` (or `mlog10p`) instead. NA rows are kept rather than dropped, so "tested, no
      meta-analysis" stays distinguishable from "gene absent from this file".

`--product` exists because the same Zenodo record also ships the gene association
sumstats, the sliding-window sumstats and the gene feature matrix. `scores` and `genes`
are implemented; the sliding-window sumstats and the feature matrix are separate work and
are not stubbed here.

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

`--product genes` resolves symbols through the exact same functions
(`read_gencode_mapping`, `read_hgnc`, `pick_gencode`, `resolve_fallback`, wrapped by the
shared `resolve_symbols`) -- one v19->current-symbol resolution for both products, run
separately on each product's own set of v19 symbols (17,263 for genes, a subset of the
scores' 18,641), so the two products can pick different winners for a symbol that is
ambiguous only because of who else is in the set, without ever forking the resolution logic
itself.

Output, scores (one bgzipped TSV, no tabix index -- there is nothing to index):
  symbol  symbol_gencode_v19  ensembl_gene_id  phaplo  ptriplo  haploinsufficient  triplosensitive

Output, genes (one bgzipped TSV, long format, no tabix index):
  dataset  phenotype  cnv_type  symbol  symbol_gencode_v19  ensembl_gene_id
  n_nominal_cohorts  top_cohort  cohorts_excluded  case_freq  control_freq
  beta  beta_lower  beta_upper  z  mlog10p  mlog10_fdr_q
  beta_secondary  beta_lower_secondary  beta_upper_secondary  z_secondary
  mlog10p_secondary  mlog10_fdr_q_secondary
See docs/rcnv-dosage-sensitivity.md for the full source-column mapping table and the
NA-row contract. Per the repo's Statistics invariants, `beta`/`beta_lower`/`beta_upper`
(and their `_secondary` twins) are formatted `:.3e` and `mlog10p`/`mlog10_fdr_q` are
rounded to 4 decimals; `z` has no such house rule and is passed through verbatim. NA stays
the literal string `NA` throughout -- unlike the scores' probabilities, these values carry
no boolean threshold downstream, so there is no verbatim-vs-rounded disagreement to protect
against.

Usage:
  python3 scripts/munge_rcnv.py --product scores --download --output out/collins_rcnv_2022_dosage_sensitivity.tsv.gz
  python3 scripts/munge_rcnv.py --product genes --download --output out/collins_rcnv_2022_gene_associations.tsv.gz
  scripts/munge_rcnv.sh                       # produce locally (PRODUCT=scores by default)
  scripts/munge_rcnv.sh --stage               # produce and publish to both buckets
  PRODUCT=genes scripts/munge_rcnv.sh --stage # gene associations, produce and publish
"""

import argparse
import csv
import gzip
import io
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

ZENODO_RECORD = "6347673"
SCORES_FILE = "Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz"
SCORES_URL = f"https://zenodo.org/records/{ZENODO_RECORD}/files/{SCORES_FILE}?download=1"

GENE_ASSOC_DIR_NAME = "Collins_rCNV_2022.gene_association_sumstats"
GENE_ASSOC_TAR = f"{GENE_ASSOC_DIR_NAME}.tar.gz"
GENE_ASSOC_URL = f"https://zenodo.org/records/{ZENODO_RECORD}/files/{GENE_ASSOC_TAR}?download=1"

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

# house names for the 21 source columns (`gene` and the GRCh37 chr/start/end are handled
# separately -- gene becomes symbol_gencode_v19, chr/start/end are dropped); order matches
# GENE_SOURCE_HEADER from n_nominal_cohorts onward
GENE_SOURCE_HEADER = [
    "#chr", "start", "end", "gene",
    "n_nominal_cohorts", "top_cohort", "cohorts_excluded_from_meta",
    "case_freq", "control_freq",
    "meta_lnOR", "meta_lnOR_lower", "meta_lnOR_upper", "meta_z",
    "meta_neg_log10_p", "meta_neg_log10_fdr_q",
    "meta_lnOR_secondary", "meta_lnOR_lower_secondary", "meta_lnOR_upper_secondary",
    "meta_z_secondary", "meta_neg_log10_p_secondary", "meta_neg_log10_fdr_q_secondary",
]

GENE_COLUMNS = [
    "dataset", "phenotype", "cnv_type", "symbol", "symbol_gencode_v19", "ensembl_gene_id",
    "n_nominal_cohorts", "top_cohort", "cohorts_excluded", "case_freq", "control_freq",
    "beta", "beta_lower", "beta_upper", "z", "mlog10p", "mlog10_fdr_q",
    "beta_secondary", "beta_lower_secondary", "beta_upper_secondary", "z_secondary",
    "mlog10p_secondary", "mlog10_fdr_q_secondary",
]

# the file name, not a column, carries phenotype and CNV type
GENE_FILE_RE = re.compile(
    r"^(?P<phenotype>HP\d+|UNKNOWN)\.rCNV\.(?P<cnv_type>DEL|DUP)\."
    r"gene_association\.meta_analysis\.stats\.bed\.gz$"
)

# input-file integrity for --product genes: 54 phenotypes x DEL/DUP, 17,263 genes each
EXPECTED_GENE_FILES = 108
EXPECTED_GENE_PHENOTYPES = 54
EXPECTED_GENES_PER_FILE = 17263
EXPECTED_GENE_ROWS = EXPECTED_GENE_FILES * EXPECTED_GENES_PER_FILE

# the gene set is identical across all 108 files (row order is not -- see munge_genes)

# the literal Zenodo file-name spelling, for the `dataset` column inside the long-format
# gene-association output; distinct from the lowercase DATASET id below, which names GCS
# paths and the BigQuery dataset/table
DATASET_LABEL = "Collins_rCNV_2022"

RESOURCE = "rcnv"
DATASET = "collins_rcnv_2022"
GCS = {
    "scores": (
        f"gs://finngen-commons/results_api_data/{RESOURCE}/{DATASET}/"
        f"{DATASET}_dosage_sensitivity.tsv.gz",
        f"gs://daly-genetics-results/{RESOURCE}/{DATASET}/"
        f"{DATASET}_dosage_sensitivity.tsv.gz",
    ),
    "genes": (
        f"gs://finngen-commons/results_api_data/{RESOURCE}/{DATASET}/"
        f"{DATASET}_gene_associations.tsv.gz",
        f"gs://daly-genetics-results/{RESOURCE}/{DATASET}/"
        f"{DATASET}_gene_associations.tsv.gz",
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


def fetch_gene_assoc(cache_dir: Path) -> Path:
    """Download and unpack the gene-association tar.gz into cache_dir; return the unpacked dir.

    Unverified against a live Zenodo download -- this host cannot reach Zenodo, so this was
    exercised only against a pre-fetched unpacked copy passed via --gene-assoc-dir. If the
    tar's own layout does not unpack to cache_dir/GENE_ASSOC_DIR_NAME, this raises with the
    tar's location so a working host can extract it by hand and pass --gene-assoc-dir.
    """
    target = cache_dir / GENE_ASSOC_DIR_NAME
    if target.exists():
        return target
    tar_path = cache_dir / GENE_ASSOC_TAR
    fetch(GENE_ASSOC_URL, tar_path)
    print(f"  extracting {tar_path}", file=sys.stderr)
    with tarfile.open(tar_path) as tf:
        # the Dockerfile pins python:3.13-slim, which has the `filter` argument (PEP 706)
        tf.extractall(cache_dir, filter="data")
    if not target.exists():
        raise SystemExit(
            f"expected {tar_path} to unpack into {target}; check its actual layout and "
            f"rerun with --gene-assoc-dir pointing at the unpacked BEDs"
        )
    return target


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


def resolve_symbols(
    v19_symbols,
    by_v19: dict[str, list[tuple[str, str | None]]],
    approved: dict[str, str],
    prev: dict[str, list[str]],
    alias: dict[str, list[str]],
) -> dict[str, tuple[str, str, str]]:
    """{v19_symbol: (current_symbol, ensembl_gene_id, route)} for one set of v19 symbols.

    The one symbol->ENSG->current-symbol resolution both --product scores and --product
    genes run -- see pick_gencode/resolve_fallback for the ordering rules. Run per-product on
    that product's own set of v19 symbols (not shared across products) because `claimed`
    collisions depend on who else is in the set; the routes and functions are identical, the
    winners for an ambiguous symbol need not be.
    """
    # symbols that keep their name are settled before any renaming route runs, so the
    # fallback can see every symbol already taken and by which ENSG. Both non-renaming
    # routes belong in that pass: Gencode's own current name, and the v19 symbol HGNC still
    # approves. Leaving the latter out is what let a prev_symbol rename land on a symbol
    # its own gene was keeping.
    picks = {v19: pick_gencode(v19, by_v19) for v19 in v19_symbols}
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

    return {v19: (settled[v19][0], picks[v19][0], settled[v19][1]) for v19 in picks}


def munge_scores(
    scores_path: Path, gencode_path: Path, hgnc_path: Path, output: Path
) -> list[list[str]]:
    scores = read_scores(scores_path)
    by_v19 = read_gencode_mapping(gencode_path)
    approved, prev, alias = read_hgnc(hgnc_path)

    resolved = resolve_symbols({v19 for v19, _, _ in scores}, by_v19, approved, prev, alias)

    rows: list[list[str]] = []
    routes: Counter[str] = Counter()
    unmapped: list[str] = []
    for v19, phaplo, ptriplo in scores:
        symbol, ensg, route = resolved[v19]
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

    write_bgzip(rows, COLUMNS, output)
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


def discover_gene_assoc_files(gene_assoc_dir: Path) -> list[Path]:
    """Sorted list of the 108 phenotype x DEL/DUP BEDs, filtered to the expected name shape."""
    matched = sorted(p for p in gene_assoc_dir.glob("*.bed.gz") if GENE_FILE_RE.match(p.name))
    if len(matched) != EXPECTED_GENE_FILES:
        found = sorted(p.name for p in gene_assoc_dir.glob("*.bed.gz"))
        raise SystemExit(
            f"{gene_assoc_dir}: expected {EXPECTED_GENE_FILES} gene-association bed.gz "
            f"files matching <phenotype>.rCNV.<DEL|DUP>.gene_association...bed.gz, found "
            f"{len(matched)} of {len(found)} .bed.gz files present: {found[:10]}"
        )
    return matched


def read_gene_assoc_file(path: Path) -> list[list[str]]:
    """Read one tabixed gene-association BED as its 21 raw columns, header and width verified."""
    with gzip.open(path, "rt") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if header != GENE_SOURCE_HEADER:
            raise SystemExit(f"{path}: unexpected header {header}")
        rows = []
        for row in reader:
            if not row:
                continue
            if len(row) != len(GENE_SOURCE_HEADER):
                raise SystemExit(
                    f"{path}:{reader.line_num}: expected {len(GENE_SOURCE_HEADER)} columns, "
                    f"got {len(row)}"
                )
            rows.append(row)
        return rows


def fmt_exp(value: str) -> str:
    """`:.3e`, matching the repo's beta-formatting invariant; NA passes through unchanged."""
    return value if value == NA else f"{float(value):.3e}"


def fmt_round4(value: str) -> str:
    """Round to 4 decimals, matching the repo's mlog10p invariant; NA passes through unchanged."""
    return value if value == NA else f"{round(float(value), 4)}"


def munge_genes(
    gene_assoc_dir: Path, gencode_path: Path, hgnc_path: Path, output: Path
) -> list[list[str]]:
    files = discover_gene_assoc_files(gene_assoc_dir)
    by_v19 = read_gencode_mapping(gencode_path)
    approved, prev, alias = read_hgnc(hgnc_path)

    phenotypes: set[str] = set()
    cnv_types: set[str] = set()
    reference_genes: list[str] | None = None
    all_v19: set[str] = set()
    parsed: list[tuple[str, str, list[list[str]]]] = []
    total_rows = 0
    for i, path in enumerate(files):
        match = GENE_FILE_RE.match(path.name)
        phenotype, cnv_type = match["phenotype"], match["cnv_type"]
        phenotypes.add(phenotype)
        cnv_types.add(cnv_type)
        data = read_gene_assoc_file(path)
        if len(data) != EXPECTED_GENES_PER_FILE:
            raise SystemExit(
                f"{path}: expected {EXPECTED_GENES_PER_FILE} rows, got {len(data)}"
            )
        genes_here = [row[3] for row in data]
        if reference_genes is None:
            reference_genes = genes_here
        # sorted, not positional: adjacent rows tied on GRCh37 (chr, start) are ordered
        # inconsistently between DEL and DUP files for the same phenotype (a source sort
        # instability, confirmed on HP0000118 -- e.g. APITD1 and PMF1-BGLAP swap by one
        # position), so row order is not actually identical across files. The gene *set*
        # is what matters here: this output drops chr/start/end entirely, so a local
        # reorder around a coordinate tie is invisible downstream.
        elif sorted(genes_here) != sorted(reference_genes):
            raise SystemExit(
                f"{path}: gene set differs from {files[0].name} -- the README claims the "
                f"same 17,263 genes in every phenotype file"
            )
        all_v19.update(genes_here)
        total_rows += len(data)
        parsed.append((phenotype, cnv_type, data))

    if len(phenotypes) != EXPECTED_GENE_PHENOTYPES:
        raise SystemExit(
            f"expected {EXPECTED_GENE_PHENOTYPES} distinct phenotypes, got "
            f"{len(phenotypes)}: {sorted(phenotypes)}"
        )
    if cnv_types != {"DEL", "DUP"}:
        raise SystemExit(f"expected cnv_type in {{DEL, DUP}} only, got {sorted(cnv_types)}")
    if total_rows != EXPECTED_GENE_ROWS:
        raise SystemExit(
            f"expected {EXPECTED_GENE_ROWS} rows ({EXPECTED_GENE_FILES} files x "
            f"{EXPECTED_GENES_PER_FILE} genes), got {total_rows}"
        )

    resolved = resolve_symbols(all_v19, by_v19, approved, prev, alias)

    rows: list[list[str]] = []
    for phenotype, cnv_type, data in parsed:
        for row in data:
            v19 = row[3]
            symbol, ensg, _route = resolved[v19]
            rows.append(
                [
                    DATASET_LABEL, phenotype, cnv_type, symbol, v19, ensg,
                    row[4], row[5], row[6], row[7], row[8],
                    fmt_exp(row[9]), fmt_exp(row[10]), fmt_exp(row[11]),
                    row[12], fmt_round4(row[13]), fmt_round4(row[14]),
                    fmt_exp(row[15]), fmt_exp(row[16]), fmt_exp(row[17]),
                    row[18], fmt_round4(row[19]), fmt_round4(row[20]),
                ]
            )

    na_rows = sum(1 for row in rows if row[11] == NA)
    ensgs_missing = sum(1 for row in rows if row[5] == NA)
    report_genes(rows, resolved, na_rows)

    if ensgs_missing:
        raise SystemExit(
            f"{ensgs_missing} rows have no ensembl_gene_id; every v19 symbol is supposed to "
            f"be in the mapping file's gene_name_19 column"
        )

    write_bgzip(rows, GENE_COLUMNS, output)
    return rows


def report_genes(
    rows: list[list[str]], resolved: dict[str, tuple[str, str, str]], na_rows: int
) -> None:
    routes: Counter[str] = Counter(route for _, _, route in resolved.values())
    hgnc = sum(v for k, v in routes.items() if k.startswith("hgnc_"))
    print("", file=sys.stderr)
    print(f"rows:                          {len(rows)}", file=sys.stderr)
    print(f"distinct v19 symbols:          {len(resolved)}", file=sys.stderr)
    print(f"symbol via gencode:            {routes['gencode']}", file=sys.stderr)
    print(f"symbol via HGNC:               {hgnc}", file=sys.stderr)
    for route in ("hgnc_approved", "hgnc_prev", "hgnc_alias"):
        print(f"  {route[5:]:<28} {routes[route]}", file=sys.stderr)
    print(f"unmapped (kept as v19 symbol): {routes['unmapped']}", file=sys.stderr)
    print(f"NA rows (no meta-analysis):    {na_rows} ({na_rows / len(rows):.4%} of rows)", file=sys.stderr)


def write_bgzip(rows: list[list[str]], columns: list[str], output: Path) -> None:
    """One bgzipped TSV, no index.

    The output family writers in sumstat_utils/peak_utils all build a tabix index over
    coordinates. Both rCNV products have none, so this writes its own bgzip pipe: the file
    is a BigQuery load file only.
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as fh:
        proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=fh)
        with io.TextIOWrapper(proc.stdin, "utf-8", newline="") as pipe:
            writer = csv.writer(pipe, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            writer.writerows(rows)
        if proc.wait() != 0:
            raise SystemExit("bgzip failed")
    print(f"\nwrote {output}", file=sys.stderr)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--product",
        choices=["scores", "genes"],
        required=True,
        help="which Zenodo product to munge: dosage-sensitivity scores, or gene-association sumstats",
    )
    parser.add_argument("--scores", type=Path, help=f"local {SCORES_FILE} (default: <cache-dir>/{SCORES_FILE})")
    parser.add_argument(
        "--gene-assoc-dir",
        type=Path,
        help=f"local unpacked {GENE_ASSOC_DIR_NAME}/ (108 bed.gz) (default: <cache-dir>/{GENE_ASSOC_DIR_NAME})",
    )
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
        help="fetch any missing input into --cache-dir (Zenodo for the scores/gene-assoc tar, "
        "the daly mapping_files bucket for the gencode mapping and HGNC set); OFF by default",
    )
    parser.add_argument("--hgnc-url", default=HGNC_URL, help="public HGNC source, used with --download when gcloud is unavailable")
    parser.add_argument("--output", type=Path, help="output .tsv.gz (default derived from --product)")
    parser.add_argument("--stage", action="store_true", help="upload the output to both profile buckets; OFF by default")
    parser.add_argument("--gcs-finngen", help="override the finngen destination")
    parser.add_argument("--gcs-daly", help="override the daly destination")
    return parser.parse_args()


def resolve_inputs(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    """(product-specific input, gencode mapping, hgnc set) -- the mapping inputs are shared."""
    gencode = args.gencode_mapping or args.cache_dir / GENCODE_MAPPING_FILE
    hgnc = args.hgnc or args.cache_dir / HGNC_FILE
    if args.download:
        fetch_gcs(f"{MAPPING_BUCKET}/{GENCODE_MAPPING_FILE}", gencode)
        if not hgnc.exists():
            try:
                fetch_gcs(f"{MAPPING_BUCKET}/{HGNC_FILE}", hgnc)
            except (subprocess.CalledProcessError, FileNotFoundError):
                print("  bucket copy failed, falling back to genenames.org", file=sys.stderr)
                fetch(args.hgnc_url, hgnc)
    for path, flag in ((gencode, "--gencode-mapping"), (hgnc, "--hgnc")):
        if not path.exists():
            raise SystemExit(f"{path} not found (pass {flag}, or run with --download)")

    if args.product == "scores":
        scores = args.scores or args.cache_dir / SCORES_FILE
        if args.download:
            fetch(SCORES_URL, scores)
        if not scores.exists():
            raise SystemExit(f"{scores} not found (pass --scores, or run with --download)")
        return scores, gencode, hgnc

    if args.gene_assoc_dir is not None:
        if not args.gene_assoc_dir.exists():
            raise SystemExit(f"--gene-assoc-dir {args.gene_assoc_dir} not found")
        return args.gene_assoc_dir, gencode, hgnc

    gene_assoc_dir = args.cache_dir / GENE_ASSOC_DIR_NAME
    if not gene_assoc_dir.exists() and args.download:
        gene_assoc_dir = fetch_gene_assoc(args.cache_dir)
    if not gene_assoc_dir.exists():
        raise SystemExit(
            f"{gene_assoc_dir} not found (pass --gene-assoc-dir, or run with --download)"
        )
    return gene_assoc_dir, gencode, hgnc


def main() -> None:
    args = parse_args()
    if shutil.which("bgzip") is None:
        raise SystemExit("bgzip not found on PATH (run inside the munge image)")

    primary, gencode, hgnc = resolve_inputs(args)
    output_suffix = "dosage_sensitivity" if args.product == "scores" else "gene_associations"
    output = args.output or Path(f"{DATASET}_{output_suffix}.tsv.gz")
    if args.product == "scores":
        munge_scores(primary, gencode, hgnc, output)
    else:
        munge_genes(primary, gencode, hgnc, output)

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
