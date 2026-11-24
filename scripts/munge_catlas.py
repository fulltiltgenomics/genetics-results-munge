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
  chrom            NUMERIC, no "chr" prefix (1..22, X->23, Y->24, M/MT->26); REQUIRED by the
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
import tempfile
from pathlib import Path

import polars as pl

from peak_utils import (
    OPEN_CHROMATIN_COLUMNS,
    build_matrix_output,
    check_sample_chroms,
    check_sample_symbols_only,
    coords_from_id,
    external_sort_bgzip_index,
    normalize_coords,
    open_chromatin_gcs_paths,
    upload_to_gcs,
    write_and_check_sample,
)

DATASET = "catlas_open_chromatin"
RESOURCE = "catlas"
VERSION = "2021"

CONDITION = "unknown"
ASSAY = "snATAC"

# canonical open_chromatin column order; first header token is comment-prefixed for tabix
OUTPUT_COLUMNS = OPEN_CHROMATIN_COLUMNS

GCS_FINNGEN, GCS_DALY = open_chromatin_gcs_paths(RESOURCE, DATASET)


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


def _normalize_matrix_coords(df: pl.DataFrame, args: argparse.Namespace) -> tuple[pl.DataFrame, list[str]]:
    """Normalize coords to numeric chrom/start/end (dropping non-canonical contigs) and return the
    frame plus the original coordinate/id source columns present (to exclude from the cell columns)."""
    if args.id_col:
        df = coords_from_id(df, args.id_col)
        coord_source_cols = [args.id_col]
    else:
        df = normalize_coords(df, args.chrom_col, args.start_col, args.end_col)
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

    external_sort_bgzip_index(body, output_path, OUTPUT_COLUMNS, sort_tmp=tmp_root)
    print(f"  wrote {total_rows} rows -> {output_path}")
    return total_rows


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
    # per-cell-type tissue + uberon (body-wide atlas: tissue is not a constant)
    tissue_map = derive_tissue_map(long["cell_type"].unique().to_list(), args)
    df = long.join(tissue_map, on="cell_type", how="left")
    # load_matrix always emits a "score" column (null for presence matrices), so no fallback needed
    return build_matrix_output(
        df, args,
        dataset=DATASET, condition=CONDITION, assay=ASSAY, score_type=args.score_type,
        cell_ontology_id=None, version=VERSION,
    )


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

    check_sample_chroms(out)

    # multi-gene cCRE: symbol and id lists must be paired positionally (SAMD11<->ENSG..634)
    multi = out.filter(pl.col("peak_id") == "1-200200-200700").row(0, named=True)
    assert multi["target_gene"] == "NOC2L,SAMD11", multi["target_gene"]
    assert multi["target_gene_id"] == "ENSG00000188976,ENSG00000187634", multi["target_gene_id"]
    genes = multi["target_gene"].split(",")
    ids = multi["target_gene_id"].split(",")
    pairs = dict(zip(genes, ids))
    assert pairs["SAMD11"] == "ENSG00000187634" and pairs["NOC2L"] == "ENSG00000188976", pairs
    print(f"  gene/id positional pairing: OK  {pairs}")

    check_sample_symbols_only(out)
    write_and_check_sample(out, str(tmpdir / f"{DATASET}.sample.tsv.gz"), tbl_width_chars=260)


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
