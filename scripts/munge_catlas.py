#!/usr/bin/env python3
"""Munge the Zhang et al. 2021 Cell "Catlas" body-wide human snATAC atlas into the canonical `open_chromatin` TSV.

Source (public):
  - Paper: Zhang, Hou, Chiou, ..., Ren, "A single-cell atlas of chromatin accessibility in the
    human genome", Cell 2021, 184(24):5985-6001. doi:10.1016/j.cell.2021.10.024
  - ~1.3M nuclei from 30 adult tissues -> ~1.15M candidate cis-regulatory elements (cCREs)
    across 222 cell types. Genome build hg38.
  - Portal: http://catlas.org/  (adult body-wide "human tissues" resource;
    mirror: https://decoder-genetics.wustl.edu/catlasv1/catlas_hub/). Provides a union cCRE BED,
    per-cell-type cCRE BEDs, a cell-type-by-cCRE accessibility matrix, signal-track bigwigs, and
    supplementary tables (nuclei, cell types, cCRE-gene links).

>>> TO-VERIFY (format assumptions) <<<
The Catlas portal is a JavaScript app and the tissue-level downloads were "coming soon" at
scrape time; the exact download URLs and per-file column layouts could not be fetched
programmatically. The inputs this script consumes and the layouts it assumes are documented
below. The user MUST verify these against the real files and adjust the --*-col flags (no code
change needed for column-name differences).

INPUT 1 — per-cell-type cCRE accessibility matrix  (--ccre-matrix, REQUIRED)
  Assumed: a TSV where each row is one cCRE and there is one column per Catlas cell type holding
  an accessibility value. The Catlas portal states signal tracks are "normalized signal per
  million reads" (CPM), so the default --score-type is cpm. cCRE coordinates identified either by
  three explicit columns (default names chrom/start/end, hg38, BED 0-based half-open) OR by a
  single id column like "chr1:1000-2000" / "chr1-1000-2000" / "chr1_1000_2000" (--id-col).
  Any column that is not a coordinate/id column is treated as a cell-type value column.
  ALTERNATIVE the user may hit instead: Catlas also publishes a BINARY cell-type-by-cCRE
  membership matrix (0/1 whether a cCRE is called in that cell type). For that pass
  --score-type presence (score is emitted empty and rows with value>0 are kept).

INPUT 2 — cell-type -> tissue mapping  (--celltype-tissue-map, OPTIONAL but RECOMMENDED)
  Catlas is body-wide, so tissue is NOT a constant. Assumed: a TSV with a cell-type column and a
  tissue column (and an optional uberon_id column). This is the RELIABLE way to assign tissue —
  build it from the Catlas cell-type supplementary table (which lists the tissue of origin for
  each of the 222 cell types). Columns are configurable via --map-*-col.
  FALLBACK when no map is given (or a cell type is absent from it): _fallback_tissue() scans the
  cell-type label for an organ/tissue keyword (e.g. "Cardiac ..." -> heart, "Alveolar ..." ->
  lung). Many Catlas labels do NOT encode tissue, so unmatched labels get tissue="unknown".
  >>> TO-VERIFY: the keyword heuristic is best-effort; prefer the explicit map. <<<

INPUT 3 — cCRE -> gene links  (--ccre-gene-links, OPTIONAL)
  Assumed: a TSV linking each cCRE to one or more target genes (Catlas derives cCRE-gene links by
  correlating accessibility with expression). Columns (configurable): a cCRE id/coords, a gene
  symbol, and an Ensembl gene id. Multiple genes per cCRE are aggregated into comma-separated
  lists, with symbol and id paired POSITIONALLY so the i-th symbol matches the i-th id.

INPUT 4 — cell-type metadata  (--cell-metadata, OPTIONAL)
  Assumed: a TSV mapping cell_type -> supporting nucleus count, used to fill `n_cells`.
  Absent -> n_cells left empty.

Genome build: hg38 (per the paper and the epic design) — NO liftOver.

OUTPUT — canonical `open_chromatin` long TSV, 18 columns IN THIS ORDER (tab-separated):
  chrom, start, end, peak_id, dataset, cell_type, tissue, life_stage, condition, assay,
  score, score_type, n_cells, cell_ontology_id, uberon_id, target_gene, target_gene_id, version
The header's first token is written as `#chrom` (comment-prefixed) following this repo's
convention (cf. munge_hpa.py `#dataset`, sumstats `#chr`): it makes the header a tabix meta
line so `tabix -p bed` skips it while the logical column name stays `chrom`.

Field rules for Catlas (see task genetics-results-suite-bzl.15):
  chrom            NUMERIC, no "chr" prefix (1..22, X->23, Y->24, M/MT->25); REQUIRED by the
                   tabix indexing contract.
  start,end        cCRE interval, hg38, BED 0-based half-open (kept verbatim from source).
  peak_id          f"{chrom}-{start}-{end}" with the numeric chrom (e.g. "23-100-200").
  dataset          "catlas_open_chromatin"  (drives resource mapping to catlas).
  cell_type        verbatim Catlas cell-type label (the matrix column name).
  tissue           per cell type, from --celltype-tissue-map, else keyword fallback, else "unknown".
  life_stage       DERIVED per cell type: "fetal" when the cell-type label (or the map's optional
                   life_stage column) marks it fetal (matches /fetal/i, e.g. "Fetal_*"), else
                   "adult". Catlas mixes adult and fetal cell types. Force all rows with --life-stage.
  condition        "unknown" (default for atlases).
  assay            "snATAC".
  score/score_type CPM accessibility value + "cpm"  (or empty + "presence" for a binary matrix).
  n_cells          from --cell-metadata if given, else empty.
  cell_ontology_id empty (not trivially mappable from the free-text Catlas label).
  uberon_id        from the map's uberon column if present, else derived from tissue via a small
                   built-in UBERON table; empty for unknown/unmapped tissues.
  target_gene(_id) from --ccre-gene-links when present, else empty.
  version          "2021".

MEMORY-BOUNDED FULL RUN (built in — no external driver needed):
  The full atlas is a ~1.2M-cCRE x 222-cell-type matrix; a single dense in-RAM unpivot is ~266M
  rows and OOMs a 31GB box. The normal matrix path instead processes the cell-type COLUMNS in
  --col-batch-size batches (default 25): each batch reads only its columns, melts+build_output's,
  and streams its 18-col body rows to an on-disk temp (--tmpdir); after all batches the concatenated
  body is header-prepended, externally sorted (LC_ALL=C sort -k1,1 -k2,2n, scratch in --tmpdir),
  bgzipped and tabix-indexed once. Output is content-identical to a single-shot run (rows differ
  only in pre-sort order, which the sort fixes). So a full 222-cell-type run is ONE command:
    python3 scripts/munge_catlas.py --ccre-matrix ccre_matrix.presence.tsv --score-type presence \
      --celltype-tissue-map celltype_tissue_map.tsv --tmpdir /mnt/disks/data/oc_munge/catlas/_tmp \
      --output catlas_open_chromatin.tsv.gz

INDEXING CONTRACT (from the epic design, surfaced by the results-api open_chromatin review):
  open_chromatin files MUST be INTERVAL-indexed: `tabix -p bed` (i.e. -s1 -b2 -e3, distinct
  start/end columns). Do NOT point-index (-b2 -e2) or the API variant-overlap fast path would
  silently miss peaks whose interval contains pos. Pipeline: sort -k1,1 -k2,2n -> bgzip -> tabix.

STAGING (off by default): with --stage the .tsv.gz and .tsv.gz.tbi are uploaded to BOTH
  gs://finngen-commons/results_api_data/open_chromatin/catlas/catlas_open_chromatin.tsv.gz
  gs://daly-genetics-results/open_chromatin/catlas/catlas_open_chromatin.tsv.gz
  These paths must match the results-api profile config files. Without --stage nothing is uploaded.

Local validation without the full dataset or any upload:
  python3 scripts/munge_catlas.py --sample
"""

import argparse
import re
import subprocess
import tempfile
from pathlib import Path

import polars as pl

DATASET = "catlas_open_chromatin"
RESOURCE = "catlas"
VERSION = "2021"

CONDITION = "unknown"
ASSAY = "snATAC"

# canonical open_chromatin column order; first header token is comment-prefixed for tabix
OUTPUT_COLUMNS = [
    "chrom", "start", "end", "peak_id", "dataset", "cell_type", "tissue", "life_stage",
    "condition", "assay", "score", "score_type", "n_cells", "cell_ontology_id",
    "uberon_id", "target_gene", "target_gene_id", "version",
]

GCS_FINNGEN = f"gs://finngen-commons/results_api_data/open_chromatin/{RESOURCE}/{DATASET}.tsv.gz"
GCS_DALY = f"gs://daly-genetics-results/open_chromatin/{RESOURCE}/{DATASET}.tsv.gz"

# "chr1:1000-2000", "chr1-1000-2000", "chr1_1000_2000"; the "chr" prefix is optional and
# gets stripped in _coords_from_id so ids like "1:1000-2000" are not dropped.
_ID_RE = re.compile(r"^(?:chr)?([^:_\-]+)[:_\-](\d+)[_\-](\d+)$")


def _numeric_chrom(expr: pl.Expr) -> pl.Expr:
    """Strip any 'chr' prefix and map to a numeric-string chrom (X=23, Y=24, M/MT=25)."""
    e = expr.cast(pl.Utf8).str.replace(r"^chr", "")
    up = e.str.to_uppercase()
    return (
        pl.when(up == "X").then(pl.lit("23"))
        .when(up == "Y").then(pl.lit("24"))
        .when(up.is_in(["M", "MT"])).then(pl.lit("25"))
        .otherwise(e)
    )


# canonical primary-assembly chromosomes as numeric strings (1..22, X=23, Y=24, M/MT=25). the
# platform is primary-assembly only and loads chrom as INT64, so non-canonical hg38 contigs
# (alt/random/scaffold/unplaced/Un_*) must be DROPPED here or they break the BigQuery chr load.
CANONICAL_CHROMS = frozenset(str(c) for c in range(1, 26))


def _drop_noncanonical(df: pl.DataFrame) -> pl.DataFrame:
    """Keep only rows whose (already numeric) chrom is a canonical primary chromosome."""
    return df.filter(pl.col("chrom").is_in(CANONICAL_CHROMS))


# best-effort keyword -> harmonized tissue for the fallback path when no explicit
# cell-type -> tissue map is provided. Matched on lowercased cell-type label, longest key first.
# >>> TO-VERIFY: heuristic only; the Catlas cell-type supplementary table is authoritative. <<<
_TISSUE_KEYWORDS = {
    "cardiac": "heart",
    "cardiomyocyte": "heart",
    "myocard": "heart",
    "alveolar": "lung",
    "pulmonary": "lung",
    "bronchial": "lung",
    "hepatocyte": "liver",
    "hepatic": "liver",
    "pancreatic": "pancreas",
    "islet": "pancreas",
    "acinar": "pancreas",
    "colon": "colon",
    "colonic": "colon",
    "enterocyte": "intestine",
    "intestinal": "intestine",
    "gastric": "stomach",
    "esophageal": "esophagus",
    "mammary": "breast",
    "cortical": "brain",
    "cerebellar": "brain",
    "astrocyte": "brain",
    "oligodendrocyte": "brain",
    "neuron": "brain",
    "renal": "kidney",
    "nephron": "kidney",
    "adrenal": "adrenal gland",
    "thyroid": "thyroid gland",
    "uterine": "uterus",
    "ovarian": "ovary",
    "prostate": "prostate",
    "skeletal muscle": "muscle",
    "muscle satellite": "muscle",
    "keratinocyte": "skin",
    "melanocyte": "skin",
    "corneal": "eye",
    "retinal": "eye",
    "schwann": "nerve",
}

# minimal harmonized-tissue -> UBERON mapping so uberon_id is trivially fillable from tissue.
# extend as needed; unmapped tissues leave uberon_id empty.
_TISSUE_UBERON = {
    "heart": "UBERON:0000948",
    "lung": "UBERON:0002048",
    "liver": "UBERON:0002107",
    "pancreas": "UBERON:0001264",
    "colon": "UBERON:0001155",
    "intestine": "UBERON:0000160",
    "stomach": "UBERON:0000945",
    "esophagus": "UBERON:0001043",
    "breast": "UBERON:0000310",
    "brain": "UBERON:0000955",
    "kidney": "UBERON:0002113",
    "adrenal gland": "UBERON:0002369",
    "thyroid gland": "UBERON:0002046",
    "uterus": "UBERON:0000995",
    "ovary": "UBERON:0000992",
    "prostate": "UBERON:0002367",
    "muscle": "UBERON:0001134",
    "skin": "UBERON:0002097",
    "eye": "UBERON:0000970",
    "nerve": "UBERON:0001021",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ccre-matrix", help="cCRE x cell-type accessibility matrix TSV (see header)")
    p.add_argument("--id-col", help="single cCRE id column to parse coords from (e.g. chr1:1000-2000)")
    p.add_argument("--chrom-col", default="chrom", help="coordinate column: chrom (default: chrom)")
    p.add_argument("--start-col", default="start", help="coordinate column: start (default: start)")
    p.add_argument("--end-col", default="end", help="coordinate column: end (default: end)")
    p.add_argument("--score-type", default="cpm", help="score_type token (default: cpm; use 'presence' for a binary matrix)")
    p.add_argument("--min-score", type=float, default=0.0, help="keep cell-type entries with value > min-score (default: 0.0)")

    p.add_argument("--col-batch-size", type=int, default=25,
                   help="process the matrix's cell-type COLUMNS in batches of this many, streaming each "
                        "batch's long rows to an on-disk temp then external-sorting, so a full 222-cell-type "
                        "matrix never triggers a dense in-RAM unpivot that OOMs (default: 25)")
    p.add_argument("--tmpdir", help="directory for the on-disk body + external sort scratch "
                                    "(default: system temp; point at a big disk for the full matrix)")

    p.add_argument("--celltype-tissue-map", help="optional cell_type -> tissue (+ optional uberon_id) TSV")
    p.add_argument("--map-celltype-col", default="cell_type", help="map: cell_type column (default: cell_type)")
    p.add_argument("--map-tissue-col", default="tissue", help="map: tissue column (default: tissue)")
    p.add_argument("--map-uberon-col", default="uberon_id", help="map: optional uberon_id column (default: uberon_id)")
    p.add_argument("--map-lifestage-col", default="life_stage", help="map: optional life_stage column (default: life_stage)")
    p.add_argument("--life-stage", help="force life_stage for ALL rows (default: derived per cell type; fetal if label/map marks it)")

    p.add_argument("--ccre-gene-links", help="optional cCRE->gene links TSV")
    p.add_argument("--links-id-col", help="links: cCRE id column to parse coords from")
    p.add_argument("--links-chrom-col", default="chrom", help="links: chrom column (default: chrom)")
    p.add_argument("--links-start-col", default="start", help="links: start column (default: start)")
    p.add_argument("--links-end-col", default="end", help="links: end column (default: end)")
    p.add_argument("--links-gene-col", default="gene_name", help="links: gene symbol column (default: gene_name)")
    p.add_argument("--links-geneid-col", default="gene_id", help="links: Ensembl gene id column (default: gene_id)")

    p.add_argument("--cell-metadata", help="optional cell_type -> n_cells TSV")
    p.add_argument("--meta-celltype-col", default="cell_type", help="metadata: cell_type column (default: cell_type)")
    p.add_argument("--meta-ncells-col", default="n_cells", help="metadata: n_cells column (default: n_cells)")

    p.add_argument("--output", help="output .tsv.gz path (default: ./catlas_open_chromatin.tsv.gz)")
    p.add_argument("--stage", action="store_true", help="upload .tsv.gz + .tbi to BOTH GCS buckets (OFF by default)")
    p.add_argument("--gcs-finngen", default=GCS_FINNGEN, help="finngen GCS destination")
    p.add_argument("--gcs-daly", default=GCS_DALY, help="daly GCS destination")

    p.add_argument("--sample", "--dry-run", dest="sample", action="store_true",
                   help="run the transform on tiny synthetic input, print first rows, validate tabix; no upload")
    return p.parse_args()


def _coords_from_id(df: pl.DataFrame, id_col: str) -> pl.DataFrame:
    """Parse chrom/start/end from a single id column like chr1:1000-2000 (chrom -> numeric)."""
    parsed = df.select(
        pl.col(id_col).str.extract_groups(_ID_RE.pattern).alias("_g")
    ).unnest("_g")
    return _drop_noncanonical(df.with_columns(
        _numeric_chrom(parsed["1"]).alias("chrom"),
        parsed["2"].cast(pl.Int64).alias("start"),
        parsed["3"].cast(pl.Int64).alias("end"),
    ))


def _normalize_coords(df: pl.DataFrame, chrom_col: str, start_col: str, end_col: str) -> pl.DataFrame:
    """Rename explicit coordinate columns to chrom/start/end and enforce numeric chrom + int coords."""
    df = df.rename({chrom_col: "chrom", start_col: "start", end_col: "end"})
    return _drop_noncanonical(df.with_columns(
        _numeric_chrom(pl.col("chrom")).alias("chrom"),
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
    ))


def _normalize_matrix_coords(df: pl.DataFrame, args: argparse.Namespace) -> tuple[pl.DataFrame, list[str]]:
    """Normalize coords to numeric chrom/start/end (dropping non-canonical contigs) and return the
    frame plus the original coordinate/id source columns present (to exclude from the cell columns)."""
    if args.id_col:
        df = _coords_from_id(df, args.id_col)
        coord_source_cols = [args.id_col]
    else:
        df = _normalize_coords(df, args.chrom_col, args.start_col, args.end_col)
        coord_source_cols = [args.chrom_col, args.start_col, args.end_col]
    coord_source_cols = [c for c in coord_source_cols if c in df.columns]
    return df, coord_source_cols


def _unpivot_to_long(df: pl.DataFrame, cell_cols: list[str], args: argparse.Namespace) -> pl.DataFrame:
    """Unpivot a wide (already coord-normalized) frame to one row per accessible (cCRE, cell_type),
    keeping entries with value > min-score and emitting the `score` column per score_type."""
    long = df.unpivot(
        index=["chrom", "start", "end"], on=cell_cols,
        variable_name="cell_type", value_name="_value",
    )
    long = long.with_columns(pl.col("_value").cast(pl.Float64, strict=False))
    long = long.filter(pl.col("_value") > args.min_score)

    if args.score_type == "presence":
        long = long.with_columns(pl.lit(None, dtype=pl.Utf8).alias("score"))
    else:
        long = long.with_columns(pl.col("_value").alias("score"))
    return long.drop("_value")


def load_matrix(args: argparse.Namespace) -> pl.DataFrame:
    """Read the cCRE x cell-type matrix and unpivot to long: one row per (cCRE, cell_type).

    In-RAM path used by --sample and small matrices; the full 222-cell-type matrix is processed by
    build_matrix_batched() instead, which avoids the dense in-RAM unpivot.
    """
    df = pl.read_csv(args.ccre_matrix, separator="\t", infer_schema_length=10_000)
    df, coord_source_cols = _normalize_matrix_coords(df, args)
    cell_cols = [c for c in df.columns if c not in {"chrom", "start", "end", *coord_source_cols}]
    if not cell_cols:
        raise ValueError("no cell-type value columns found in --ccre-matrix (check coordinate/id flags)")
    return _unpivot_to_long(df, cell_cols, args)


def build_matrix_batched(args: argparse.Namespace, output_path: str) -> int:
    """Memory-bounded matrix transform: process the cell-type COLUMNS in --col-batch-size batches,
    melting+building each batch and streaming its 18-col body rows to an on-disk temp, then externally
    sort -> bgzip -> tabix once. Avoids the dense (n_cCRE x n_cell_type) in-RAM unpivot that OOMs on a
    full 222-cell-type matrix. build_output (tissue/uberon/life_stage derivation, gene links, metadata,
    canonical 18-col order) is reused UNCHANGED per batch, so output content matches a single-shot run
    (rows differ only in pre-sort order, which the external sort fixes)."""
    header_cols = pl.read_csv(args.ccre_matrix, separator="\t", n_rows=0).columns
    coord_cols = [args.id_col] if args.id_col else [args.chrom_col, args.start_col, args.end_col]
    coord_cols = [c for c in coord_cols if c in header_cols]
    cell_cols = [c for c in header_cols if c not in {"chrom", "start", "end", *coord_cols}]
    if not cell_cols:
        raise ValueError("no cell-type value columns found in --ccre-matrix (check coordinate/id flags)")

    tmp_root = Path(args.tmpdir) if args.tmpdir else None
    if tmp_root is not None:
        tmp_root.mkdir(parents=True, exist_ok=True)
    tmpdir = Path(tempfile.mkdtemp(dir=tmp_root))
    body = tmpdir / "body.tsv"
    batch_size = max(1, args.col_batch_size)
    total_rows = 0
    print(f"  {len(cell_cols)} cell-type columns; batching {batch_size} at a time (tmpdir={tmpdir})")
    with open(body, "wb") as bfh:
        for start in range(0, len(cell_cols), batch_size):
            batch = cell_cols[start:start + batch_size]
            sub = pl.read_csv(args.ccre_matrix, separator="\t", columns=[*coord_cols, *batch],
                              infer_schema_length=10_000)
            sub, coord_source_cols = _normalize_matrix_coords(sub, args)
            long = _unpivot_to_long(sub, batch, args)
            out = build_output(long, args)
            assert out.columns == OUTPUT_COLUMNS, f"column order mismatch: {out.columns}"
            bfh.write(out.write_csv(separator="\t", include_header=False, null_value="NA").encode())
            total_rows += out.height
            print(f"    cols {start}-{start + len(batch) - 1}: +{out.height} rows (cum {total_rows})")

    _sort_index_bgzip(body, output_path, tmp_root)
    print(f"  wrote {total_rows} rows -> {output_path}")
    print(f"  indexed {output_path}.tbi (tabix -p bed / -s1 -b2 -e3)")
    return total_rows


def load_gene_links(args: argparse.Namespace) -> pl.DataFrame:
    """Read cCRE->gene links and aggregate to one comma-joined gene list per cCRE key."""
    df = pl.read_csv(args.ccre_gene_links, separator="\t", infer_schema_length=10_000)
    if args.links_id_col:
        df = _coords_from_id(df, args.links_id_col)
    else:
        df = _normalize_coords(df, args.links_chrom_col, args.links_start_col, args.links_end_col)

    df = df.with_columns(_peak_key().alias("_key"))
    gene = pl.col(args.links_gene_col).cast(pl.Utf8) if args.links_gene_col in df.columns else pl.lit(None, dtype=pl.Utf8)
    geneid = pl.col(args.links_geneid_col).cast(pl.Utf8) if args.links_geneid_col in df.columns else pl.lit(None, dtype=pl.Utf8)
    # normalize nulls to empty strings so a missing symbol or id still keeps its pair and
    # emits an empty token on the missing side, keeping the two comma-joined lists aligned
    df = df.with_columns(
        gene.fill_null("").alias("_gene"),
        geneid.fill_null("").alias("_geneid"),
    )

    # aggregate the (symbol, id) pairs JOINTLY: build a struct per link, dedup the pairs, and
    # sort the structs (lexicographic by _gene first, then _geneid) so the i-th token of
    # target_gene always corresponds to the i-th token of target_gene_id
    agg = (
        df.group_by("_key")
        .agg(pl.struct("_gene", "_geneid").unique().sort().alias("_pairs"))
        .with_columns(
            pl.col("_pairs").list.eval(pl.element().struct.field("_gene")).list.join(",").alias("target_gene"),
            pl.col("_pairs").list.eval(pl.element().struct.field("_geneid")).list.join(",").alias("target_gene_id"),
        )
        .drop("_pairs")
    )
    # empty-cell -> NA convention: a joined list with no real token (only commas/whitespace, e.g.
    # when the source provides gene symbols but no Ensembl ids) becomes NA, not "" / "," / ",,".
    # applied symmetrically; rows that DO have real tokens keep their positional pairing intact.
    return agg.with_columns(
        *[
            pl.when(pl.col(c).str.replace_all(r"[,\s]", "").str.len_chars() == 0)
            .then(pl.lit(None, dtype=pl.Utf8))
            .otherwise(pl.col(c))
            .alias(c)
            for c in ("target_gene", "target_gene_id")
        ]
    )


def _peak_key() -> pl.Expr:
    return pl.format("{}-{}-{}", pl.col("chrom"), pl.col("start"), pl.col("end"))


def _derive_life_stage(cell_type: str, provided: str | None) -> str:
    """fetal when the map's life_stage value (authoritative) or the cell-type label marks it
    fetal (matches /fetal/i, e.g. a "Fetal_*" label); else adult."""
    val = provided if provided else cell_type
    return "fetal" if re.search(r"fetal", val, re.IGNORECASE) else "adult"


def _fallback_tissue(cell_type: str) -> str:
    """Best-effort tissue from a Catlas cell-type label; 'unknown' when no keyword matches."""
    label = cell_type.lower()
    # longest keyword first so e.g. "skeletal muscle" wins over a bare "muscle" substring
    for kw in sorted(_TISSUE_KEYWORDS, key=len, reverse=True):
        if kw in label:
            return _TISSUE_KEYWORDS[kw]
    return "unknown"


def derive_tissue_map(cell_types: list[str], args: argparse.Namespace) -> pl.DataFrame:
    """Build a cell_type -> (tissue, uberon_id, life_stage) table.

    Priority per cell type: explicit --celltype-tissue-map, else keyword fallback, else 'unknown'.
    uberon_id: the map's uberon column if provided, else derived from the harmonized tissue.
    life_stage: --life-stage override, else map life_stage column, else the cell-type label
    (fetal if it marks fetal, e.g. "Fetal_*"), else adult.
    """
    provided_tissue: dict[str, str] = {}
    provided_uberon: dict[str, str] = {}
    provided_life: dict[str, str] = {}
    if args.celltype_tissue_map:
        m = pl.read_csv(args.celltype_tissue_map, separator="\t", infer_schema_length=10_000)
        has_uberon = args.map_uberon_col in m.columns
        has_life = args.map_lifestage_col in m.columns
        for row in m.iter_rows(named=True):
            ct = str(row[args.map_celltype_col])
            provided_tissue[ct] = str(row[args.map_tissue_col])
            if has_uberon and row.get(args.map_uberon_col) not in (None, ""):
                provided_uberon[ct] = str(row[args.map_uberon_col])
            if has_life and row.get(args.map_lifestage_col) not in (None, ""):
                provided_life[ct] = str(row[args.map_lifestage_col])

    rows = []
    for ct in cell_types:
        tissue = provided_tissue.get(ct) or _fallback_tissue(ct)
        uberon = provided_uberon.get(ct) or _TISSUE_UBERON.get(tissue)
        life = args.life_stage or _derive_life_stage(ct, provided_life.get(ct))
        rows.append((ct, tissue, uberon, life))
    return pl.DataFrame(
        rows,
        schema={"cell_type": pl.Utf8, "tissue": pl.Utf8, "uberon_id": pl.Utf8, "life_stage": pl.Utf8},
        orient="row",
    )


def build_output(long: pl.DataFrame, args: argparse.Namespace) -> pl.DataFrame:
    """Assemble the 18 canonical columns in order from the long (cCRE, cell_type) table."""
    df = long.with_columns(_peak_key().alias("peak_id"))

    # per-cell-type tissue + uberon (body-wide atlas: tissue is not a constant)
    tissue_map = derive_tissue_map(df["cell_type"].unique().to_list(), args)
    df = df.join(tissue_map, on="cell_type", how="left")

    if args.ccre_gene_links:
        links = load_gene_links(args)
        df = df.with_columns(_peak_key().alias("_key")).join(links, on="_key", how="left").drop("_key")
    else:
        df = df.with_columns(
            pl.lit(None, dtype=pl.Utf8).alias("target_gene"),
            pl.lit(None, dtype=pl.Utf8).alias("target_gene_id"),
        )

    if args.cell_metadata:
        meta = pl.read_csv(args.cell_metadata, separator="\t", infer_schema_length=10_000)
        meta = meta.select(
            pl.col(args.meta_celltype_col).cast(pl.Utf8).alias("cell_type"),
            pl.col(args.meta_ncells_col).cast(pl.Utf8).alias("n_cells"),
        )
        df = df.join(meta, on="cell_type", how="left")
    else:
        df = df.with_columns(pl.lit(None, dtype=pl.Utf8).alias("n_cells"))

    df = df.with_columns(
        pl.lit(DATASET).alias("dataset"),
        pl.lit(CONDITION).alias("condition"),
        pl.lit(ASSAY).alias("assay"),
        pl.lit(args.score_type).alias("score_type"),
        pl.lit(None, dtype=pl.Utf8).alias("cell_ontology_id"),
        pl.lit(VERSION).alias("version"),
    )
    # note: load_matrix always emits a "score" column (null for presence matrices), so no fallback needed

    return df.select(OUTPUT_COLUMNS)


def _sort_index_bgzip(body: Path, output_path: str, sort_tmp: Path | None = None) -> None:
    """External LC_ALL=C sort -k1,1 -k2,2n of an on-disk body TSV, prepend the header, bgzip, and
    tabix -p bed (interval index). sort_tmp keeps the sort's scratch off a tiny root disk when set."""
    header = "#" + "\t".join(OUTPUT_COLUMNS)
    tflag = f"-T {shell_quote(str(sort_tmp))} " if sort_tmp is not None else ""
    # prepend header AFTER sorting the body (header must stay on top, not be sorted)
    pipeline = (
        f'( printf "%s\\n" {shell_quote(header)}; '
        f'LC_ALL=C sort {tflag}-k1,1 -k2,2n {shell_quote(str(body))} ) | bgzip -c > {shell_quote(output_path)}'
    )
    subprocess.run(pipeline, shell=True, check=True, executable="/bin/bash")
    subprocess.run(["tabix", "-f", "-p", "bed", output_path], check=True)


def write_open_chromatin(df: pl.DataFrame, output_path: str) -> None:
    """sort -k1,1 -k2,2n -> bgzip -> tabix -p bed (interval index), missing values as "NA"."""
    tmpdir = tempfile.mkdtemp()
    body = Path(tmpdir) / "body.tsv"
    # write body without header; every empty/missing cell serialized as the literal "NA"
    df.write_csv(body, separator="\t", include_header=False, null_value="NA")
    _sort_index_bgzip(body, output_path)
    print(f"  wrote {df.height} rows -> {output_path}")
    print(f"  indexed {output_path}.tbi (tabix -p bed / -s1 -b2 -e3)")


def shell_quote(s: str) -> str:
    return "'" + s.replace("'", "'\\''") + "'"


def upload_to_gcs(local_path: str, gcs_path: str) -> None:
    subprocess.run(["gcloud", "storage", "cp", local_path, gcs_path], check=True)
    print(f"  uploaded {gcs_path}")
    tbi = local_path + ".tbi"
    if Path(tbi).exists():
        subprocess.run(["gcloud", "storage", "cp", tbi, gcs_path + ".tbi"], check=True)
        print(f"  uploaded {gcs_path}.tbi")


def _synthetic_inputs(tmpdir: Path) -> tuple[str, str, str, str]:
    """Write a tiny synthetic matrix + tissue map + gene links + metadata mimicking the layouts.

    Exercises: a multi-tissue cell-type set (heart/liver/lung/unknown, hitting all three tissue
    paths: explicit map, keyword fallback, and unmapped->unknown) and a multi-gene cCRE.
    """
    matrix = tmpdir / "ccre_matrix.tsv"
    # "Fetal_Hepatocyte" exercises the per-cell-type life_stage derivation (-> fetal); the other
    # labels are adult. Catlas mixes adult and fetal cell types, so life_stage is NOT a constant.
    matrix.write_text(
        "chrom\tstart\tend\tCardiomyocyte\tHepatocyte\tAlveolar Type 2\tMystery Cell\tFetal_Hepatocyte\n"
        "chr1\t100100\t100600\t0.0\t3.5\t0.0\t1.0\t2.5\n"
        "chr1\t200200\t200700\t7.2\t0.0\t1.1\t0.0\t0.0\n"
        "chr2\t50050\t50550\t0.0\t0.0\t4.4\t0.0\t0.0\n"
        "chrX\t900900\t901400\t2.0\t0.0\t0.0\t0.0\t0.0\n"
        # non-canonical hg38 contig (alt/random/scaffold) — MUST be dropped by the canonical filter
        "chr1_KI270706v1_random\t300300\t300800\t5.0\t0.0\t0.0\t0.0\t0.0\n"
    )
    # explicit map covers heart/liver (with uberon for heart); Alveolar Type 2 -> lung via
    # keyword fallback; "Mystery Cell" is unmapped and has no keyword -> tissue "unknown"
    ct_map = tmpdir / "celltype_tissue_map.tsv"
    ct_map.write_text(
        "cell_type\ttissue\tuberon_id\n"
        "Cardiomyocyte\theart\tUBERON:0000948\n"
        "Hepatocyte\tliver\t\n"
    )
    links = tmpdir / "ccre_gene_links.tsv"
    links.write_text(
        "chrom\tstart\tend\tgene_name\tgene_id\n"
        "chr1\t200200\t200700\tSAMD11\tENSG00000187634\n"
        "chr1\t200200\t200700\tNOC2L\tENSG00000188976\n"
        "chr2\t50050\t50550\tSOX11\tENSG00000176887\n"
        # symbols but NO Ensembl ids -> target_gene keeps the symbols, target_gene_id -> NA
        "chr1\t100100\t100600\tFOXA1\t\n"
        "chr1\t100100\t100600\tGATA3\t\n"
    )
    meta = tmpdir / "cell_metadata.tsv"
    meta.write_text(
        "cell_type\tn_cells\n"
        "Cardiomyocyte\t45000\n"
        "Hepatocyte\t30000\n"
        "Alveolar Type 2\t22000\n"
        "Mystery Cell\t500\n"
        "Fetal_Hepatocyte\t15000\n"
    )
    return str(matrix), str(ct_map), str(links), str(meta)


def run_sample() -> None:
    print("=== SAMPLE / DRY-RUN: synthetic input, no GCS upload ===")
    tmpdir = Path(tempfile.mkdtemp())
    matrix, ct_map, links, meta = _synthetic_inputs(tmpdir)

    args = parse_args()
    args.ccre_matrix = matrix
    args.id_col = None
    args.chrom_col, args.start_col, args.end_col = "chrom", "start", "end"
    args.celltype_tissue_map = ct_map
    args.map_celltype_col, args.map_tissue_col, args.map_uberon_col = "cell_type", "tissue", "uberon_id"
    args.map_lifestage_col = "life_stage"
    args.life_stage = None
    args.ccre_gene_links = links
    args.links_id_col = None
    args.links_chrom_col, args.links_start_col, args.links_end_col = "chrom", "start", "end"
    args.links_gene_col, args.links_geneid_col = "gene_name", "gene_id"
    args.cell_metadata = meta
    args.meta_celltype_col, args.meta_ncells_col = "cell_type", "n_cells"
    args.score_type = "cpm"
    args.min_score = 0.0

    long = load_matrix(args)
    out = build_output(long, args)

    assert out.columns == OUTPUT_COLUMNS, f"column order mismatch: {out.columns}"
    print(f"  output has {len(out.columns)} columns in canonical order: OK")

    # tissue derivation covers all three paths
    tissue_by_ct = {r["cell_type"]: r["tissue"] for r in out.select("cell_type", "tissue").unique().iter_rows(named=True)}
    assert tissue_by_ct["Cardiomyocyte"] == "heart", tissue_by_ct        # from explicit map
    assert tissue_by_ct["Hepatocyte"] == "liver", tissue_by_ct           # from explicit map
    assert tissue_by_ct["Alveolar Type 2"] == "lung", tissue_by_ct       # from keyword fallback
    assert tissue_by_ct["Mystery Cell"] == "unknown", tissue_by_ct       # unmapped -> unknown
    print(f"  tissue derivation (map + fallback + unknown): OK  {tissue_by_ct}")

    # life_stage derived per cell type: a Fetal_* label -> fetal, other labels -> adult (fix)
    life_by_ct = {r["cell_type"]: r["life_stage"] for r in out.select("cell_type", "life_stage").unique().iter_rows(named=True)}
    assert life_by_ct["Fetal_Hepatocyte"] == "fetal", life_by_ct
    assert life_by_ct["Cardiomyocyte"] == "adult", life_by_ct
    assert life_by_ct["Hepatocyte"] == "adult", life_by_ct
    print(f"  life_stage derivation (fetal vs adult): OK  {life_by_ct}")

    # chrom / peak_id are numeric with chrX -> 23
    chroms = set(out["chrom"].to_list())
    assert chroms <= {"1", "2", "23"}, chroms
    assert "23" in chroms, f"chrX should map to 23; got {chroms}"
    xrow = out.filter(pl.col("chrom") == "23").row(0, named=True)
    assert xrow["peak_id"] == "23-900900-901400", xrow["peak_id"]
    print(f"  numeric chrom (chrX -> 23): OK  chroms={sorted(chroms)}  peak_id={xrow['peak_id']}")

    # non-canonical hg38 contigs (alt/random/scaffold/Un) are DROPPED (Fix A) or they break the
    # INT64 chr load; the "chr1_KI270706v1_random" synthetic cCRE must not survive while chrX does
    assert not any("KI270706" in p for p in out["peak_id"].to_list()), "non-canonical contig leaked"
    assert chroms == {"1", "2", "23"}, f"expected canonical chroms only; got {chroms}"
    print("  non-canonical contig (chr1_KI270706v1_random) dropped: OK")

    # multi-gene cCRE: symbol and id lists must be paired positionally (SAMD11<->ENSG..634)
    multi = out.filter(pl.col("peak_id") == "1-200200-200700").row(0, named=True)
    assert multi["target_gene"] == "NOC2L,SAMD11", multi["target_gene"]
    assert multi["target_gene_id"] == "ENSG00000188976,ENSG00000187634", multi["target_gene_id"]
    genes = multi["target_gene"].split(",")
    ids = multi["target_gene_id"].split(",")
    pairs = dict(zip(genes, ids))
    assert pairs["SAMD11"] == "ENSG00000187634" and pairs["NOC2L"] == "ENSG00000188976", pairs
    print(f"  gene/id positional pairing: OK  {pairs}")

    # symbols-but-no-ids link: target_gene keeps the symbols, target_gene_id renders as NA (fix)
    sym_only = out.filter(pl.col("peak_id") == "1-100100-100600").row(0, named=True)
    assert sym_only["target_gene"] == "FOXA1,GATA3", sym_only["target_gene"]
    assert sym_only["target_gene_id"] is None, sym_only["target_gene_id"]
    print(f"  symbols-without-ids -> target_gene={sym_only['target_gene']!r}, target_gene_id=NA: OK")

    print("\nfirst output rows:")
    with pl.Config(tbl_cols=-1, tbl_width_chars=260, fmt_str_lengths=60):
        print(out.head())

    out_path = str(tmpdir / f"{DATASET}.sample.tsv.gz")
    write_open_chromatin(out, out_path)

    # empty cells serialize as the literal "NA"; every row still has 18 columns
    import gzip
    with gzip.open(out_path, "rt") as fh:
        body_lines = [ln for ln in fh if not ln.startswith("#")]
    assert all(len(ln.rstrip("\n").split("\t")) == 18 for ln in body_lines), "not 18 columns"
    assert any("\tNA\t" in ln or ln.rstrip("\n").endswith("\tNA") for ln in body_lines), "no NA emitted"
    print(f"  empty cells rendered as NA, 18 columns per row: OK ({len(body_lines)} rows)")

    # the symbols-without-ids cCRE serializes target_gene_id as the literal NA (not "" or ",")
    na_fields = next(ln.rstrip("\n").split("\t") for ln in body_lines if "1-100100-100600" in ln)
    assert na_fields[15] == "FOXA1,GATA3" and na_fields[16] == "NA", na_fields[15:17]
    print(f"  serialized symbols/NA pair: target_gene={na_fields[15]!r}, target_gene_id={na_fields[16]!r}: OK")

    # validate the interval index answers overlap queries (numeric seqnames, incl. chrX->23)
    q = subprocess.run(["tabix", out_path, "1:200300-200400"], capture_output=True, text=True, check=True)
    print("\ntabix overlap query 1:200300-200400 ->")
    print(q.stdout.rstrip() or "  (no rows)")
    assert q.stdout.strip(), "interval query returned nothing — indexing is wrong"
    qx = subprocess.run(["tabix", out_path, "23:901000-901100"], capture_output=True, text=True, check=True)
    print("tabix overlap query 23:901000-901100 (chrX) ->")
    print(qx.stdout.rstrip() or "  (no rows)")
    assert qx.stdout.strip(), "chrX (23) interval query returned nothing — indexing is wrong"
    print("\n=== SAMPLE OK (no upload performed) ===")


def main() -> None:
    args = parse_args()
    if args.sample:
        run_sample()
        return

    if not args.ccre_matrix:
        raise SystemExit("--ccre-matrix is required (or use --sample). See --help.")

    output_path = args.output or f"./{DATASET}.tsv.gz"
    print(f"Reading matrix {args.ccre_matrix} (column-batched, --col-batch-size={args.col_batch_size}) ...")
    # memory-bounded: process cell-type columns in batches to an on-disk temp then external-sort, so
    # a full 222-cell-type matrix never triggers the dense in-RAM unpivot that OOMs. This is the normal
    # path for both small and full matrices (the batch loop degenerates to a single batch when the
    # matrix has <= --col-batch-size cell-type columns).
    build_matrix_batched(args, output_path)

    if args.stage:
        print("Staging to GCS (both buckets) ...")
        upload_to_gcs(output_path, args.gcs_finngen)
        upload_to_gcs(output_path, args.gcs_daly)
    else:
        print("  --stage not set: skipping GCS upload")
    print("Done.")


if __name__ == "__main__":
    main()
