#!/usr/bin/env python3
"""Munge the EpiMap (Boix et al. 2021 Nature) ChromHMM 18-state calls into the canonical `open_chromatin` TSV.

Source (public):
  - Paper: Boix, James, Park, Meuleman, Kellis, "Regulatory genomic circuitry of human disease
    loci by integrative epigenomics", Nature 2021, 590(7845):300-307. doi:10.1038/s41586-020-03145-z
  - EpiMap = 833 biosamples, integrative (observed + imputed) epigenomes. We take ONLY its
    value-add ChromHMM 18-state regulatory segmentation (the "observed_aux_18" model), in hg38.
  - Broad downloads: https://personal.broadinstitute.org/cboix/epimap/
      per-biosample hg38 ChromHMM CALLS:
        ChromHMM/observed_aux_18_hg38/CALLS/BSS<id>_18_CALLS_segments.bed.gz
      sample metadata:
        metadata/main_metadata_table.tsv  (columns: id, GROUP, ct, SECONDARY, ..., lifestage, ...)

>>> TO-VERIFY (source assumptions) <<<
  - CALLS: use the hg38 model directory `observed_aux_18_hg38/CALLS/` (NOT the hg19 CALLS).
    Build hg38 natively -> NO liftOver. Confirm the directory/build at download.
  - Each CALLS file is a BED9 ChromHMM segmentation for ONE biosample:
      chrom  start  end  state  0  .  thickStart  thickEnd  itemRgb
    The state label is column 4 (1-based; --state-col, default 4) and here appears in NAME form
    (e.g. "EnhA2", "TssA", "Quies"). `normalize_state` also accepts numeric forms (1..18, "E10",
    "10_EnhA2") in case a different export is used. chrom is already chr-prefixed; hg38, BED
    0-based half-open coords are kept verbatim.
  - Do NOT import the EpiMap DNase layer (that is the dropped Meuleman DHS index) — only the
    ChromHMM state CALLS.

STATE MODELING (from task genetics-results-suite-bzl.4 design decision):
  Model = Roadmap/EpiMap 18-state observed_aux_18. Documented canonical state order (1..18):
    1 TssA, 2 TssFlnk, 3 TssFlnkU, 4 TssFlnkD, 5 Tx, 6 TxWk, 7 EnhG1, 8 EnhG2, 9 EnhA1,
    10 EnhA2, 11 EnhWk, 12 ZNF_Rpts, 13 Het, 14 TssBiv, 15 EnhBiv, 16 ReprPC, 17 ReprPCWk, 18 Quies
  INCLUDE (accessible/active regulatory subset) -> emitted as open_chromatin rows:
    active promoter : TssA, TssFlnk, TssFlnkU, TssFlnkD
    active enhancer : EnhA1, EnhA2, EnhG1, EnhG2, EnhWk
    poised/bivalent : TssBiv, EnhBiv
  EXCLUDE (not open peaks; segmentation rows are dropped):
    Tx, TxWk, ZNF_Rpts, Het, ReprPC, ReprPCWk, Quies
  Encoding (NO schema change): the exact state goes into score_type as `chromhmm:<State>`
  (e.g. score_type="chromhmm:EnhA1"); assay="chromHMM"; score is EMPTY/NULL (presence-based —
  ChromHMM is a categorical segmentation, so row existence = that state was called in that
  biosample). Agent-facing regulatory_class grouping (Tss*->promoter, Enh*->enhancer,
  *Biv->bivalent) is documented in the MCP tool description, NOT stored here.

INPUTS
  CALLS files (REQUIRED): positional paths and/or --calls-dir (globbed by --calls-glob).
    One file per biosample. The biosample id (-> cell_type) is parsed from the filename
    (the "BSS<digits>" token, else the filename stem); override a single file with --cell-type.
  --biosample-metadata (OPTIONAL): biosample -> tissue [-> life_stage] TSV. Default column
    names biosample/tissue/life_stage; point them at the raw EpiMap table with
    --meta-biosample-col id --meta-tissue-col SECONDARY --meta-lifestage-col lifestage.
    Absent or missing biosample -> tissue="unknown". life_stage="adult" unless the metadata
    value marks the sample fetal/embryonic (then "fetal").
  --links (OPTIONAL): EpiMap activity-based enhancer->gene links. Joined to segments by exact
    interval key (chrom-start-end) with JOINT symbol/id pairing (mirrors munge_li_brain.py: the
    pairs are deduped and sorted TOGETHER, never independently, so the i-th target_gene token
    matches the i-th target_gene_id token). >>> TO-VERIFY: EpiMap links are keyed by enhancer
    intervals that may not coincide exactly with ChromHMM segment boundaries; supply a links
    file already keyed to the segment intervals, or extend this to an interval-overlap join. <<<

OUTPUT — canonical `open_chromatin` long TSV, 18 columns IN THIS ORDER (tab-separated):
  chrom, start, end, peak_id, dataset, cell_type, tissue, life_stage, condition, assay,
  score, score_type, n_cells, cell_ontology_id, uberon_id, target_gene, target_gene_id, version
The header's first token is written as `#chrom` (comment-prefixed) following this repo's
convention (cf. munge_hpa.py `#dataset`, sumstats `#chr`): it makes the header a tabix meta line
so `tabix -p bed` skips it while the logical column name stays `chrom`.

Field rules for EpiMap:
  chrom            NUMERIC, no "chr" prefix (1..22, X->23, Y->24, M/MT->25); REQUIRED by the
                   tabix indexing contract.
  start,end        ChromHMM segment interval, hg38, BED 0-based half-open (verbatim from source).
  peak_id          f"{chrom}-{start}-{end}" with the numeric chrom (e.g. "23-100-200").
  dataset          "epimap_open_chromatin"  (drives resource mapping to epimap).
  cell_type        verbatim EpiMap biosample id (e.g. "BSS00001").
  tissue           from --biosample-metadata, else "unknown".
  life_stage       "adult" unless the biosample metadata marks it fetal/embryonic (then "fetal").
  condition        "unknown".
  assay            "chromHMM".
  score            EMPTY/NULL (presence-based; ChromHMM is categorical).
  score_type       f"chromhmm:{State}"  (e.g. "chromhmm:EnhA1").
  n_cells          empty (bulk integrative epigenome, not single-cell).
  cell_ontology_id empty.
  uberon_id        derived from tissue via a small built-in UBERON table when trivially mappable,
                   else empty.
  target_gene(_id) from --links when present, else empty.
  version          "2021".

INDEXING CONTRACT (from the epic design, surfaced by the results-api open_chromatin review):
  open_chromatin files MUST be INTERVAL-indexed: `tabix -p bed` (i.e. -s1 -b2 -e3, distinct
  start/end columns). Do NOT point-index (-b2 -e2) or the API variant-overlap fast path would
  silently miss peaks whose interval contains pos. Pipeline: sort -k1,1 -k2,2n -> bgzip -> tabix.

STAGING (off by default): with --stage the .tsv.gz and .tsv.gz.tbi are uploaded to BOTH
  gs://finngen-commons/results_api_data/open_chromatin/epimap/epimap_open_chromatin.tsv.gz
  gs://daly-genetics-results/open_chromatin/epimap/epimap_open_chromatin.tsv.gz
  These paths must match the results-api profile config files. Without --stage nothing is uploaded.

Local validation without the full dataset or any upload:
  python3 scripts/munge_epimap.py --sample
"""

import argparse
import re
import subprocess
import tempfile
from pathlib import Path

import polars as pl

DATASET = "epimap_open_chromatin"
RESOURCE = "epimap"
VERSION = "2021"

CONDITION = "unknown"
ASSAY = "chromHMM"

# canonical open_chromatin column order; first header token is comment-prefixed for tabix
OUTPUT_COLUMNS = [
    "chrom", "start", "end", "peak_id", "dataset", "cell_type", "tissue", "life_stage",
    "condition", "assay", "score", "score_type", "n_cells", "cell_ontology_id",
    "uberon_id", "target_gene", "target_gene_id", "version",
]

GCS_FINNGEN = f"gs://finngen-commons/results_api_data/open_chromatin/{RESOURCE}/{DATASET}.tsv.gz"
GCS_DALY = f"gs://daly-genetics-results/open_chromatin/{RESOURCE}/{DATASET}.tsv.gz"

# EpiMap observed_aux_18 canonical state order (index 1..18); used to resolve numeric labels.
STATE_ORDER = [
    "TssA", "TssFlnk", "TssFlnkU", "TssFlnkD", "Tx", "TxWk",
    "EnhG1", "EnhG2", "EnhA1", "EnhA2", "EnhWk", "ZNF_Rpts",
    "Het", "TssBiv", "EnhBiv", "ReprPC", "ReprPCWk", "Quies",
]

# accessible/active regulatory states kept as open_chromatin rows (all others dropped).
INCLUDED_STATES = {
    "TssA", "TssFlnk", "TssFlnkU", "TssFlnkD",   # active promoter
    "EnhA1", "EnhA2", "EnhG1", "EnhG2", "EnhWk",  # active enhancer
    "TssBiv", "EnhBiv",                            # poised/bivalent
}

_CANON_BY_LOWER = {s.lower(): s for s in STATE_ORDER}
# label spelling variants seen in the wild (e.g. "ZNF/Rpts") -> canonical name
_STATE_ALIASES = {
    "znf/rpts": "ZNF_Rpts",
    "znfrpts": "ZNF_Rpts",
    "znf_rpt": "ZNF_Rpts",
}

# minimal harmonized-tissue -> UBERON mapping so uberon_id is trivially fillable from tissue.
# lookup is case-insensitive; unmapped tissues leave uberon_id empty. Extend as needed.
_TISSUE_UBERON = {
    "brain": "UBERON:0000955",
    "heart": "UBERON:0000948",
    "lung": "UBERON:0002048",
    "liver": "UBERON:0002107",
    "pancreas": "UBERON:0001264",
    "kidney": "UBERON:0002113",
    "blood": "UBERON:0000178",
    "immune": "UBERON:0002405",
    "spleen": "UBERON:0002106",
    "thymus": "UBERON:0002370",
    "muscle": "UBERON:0001134",
    "skin": "UBERON:0002097",
    "intestine": "UBERON:0000160",
    "stomach": "UBERON:0000945",
    "breast": "UBERON:0000310",
    "ovary": "UBERON:0000992",
    "testis": "UBERON:0000473",
    "placenta": "UBERON:0001987",
    "adipose": "UBERON:0001013",
    "bone": "UBERON:0002481",
    "eye": "UBERON:0000970",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("calls", nargs="*", help="per-biosample ChromHMM CALLS BED(.gz) files")
    p.add_argument("--calls-dir", help="directory of CALLS files (globbed by --calls-glob)")
    p.add_argument("--calls-glob", default="BSS*_18_CALLS_segments.bed.gz",
                   help="glob for --calls-dir (default: BSS*_18_CALLS_segments.bed.gz)")
    p.add_argument("--state-col", type=int, default=4,
                   help="1-based column holding the ChromHMM state in the BED (default: 4)")
    p.add_argument("--cell-type", help="override biosample id (only meaningful for a single CALLS file)")

    p.add_argument("--biosample-metadata", help="biosample -> tissue [-> life_stage] TSV")
    p.add_argument("--meta-biosample-col", default="biosample", help="metadata: biosample id column (default: biosample)")
    p.add_argument("--meta-tissue-col", default="tissue", help="metadata: tissue column (default: tissue)")
    p.add_argument("--meta-lifestage-col", default="life_stage", help="metadata: optional life_stage column (default: life_stage)")

    p.add_argument("--links", help="optional enhancer->gene links TSV")
    p.add_argument("--links-id-col", help="links: interval id column to parse coords from (e.g. chr1:1000-2000)")
    p.add_argument("--links-chrom-col", default="chrom", help="links: chrom column (default: chrom)")
    p.add_argument("--links-start-col", default="start", help="links: start column (default: start)")
    p.add_argument("--links-end-col", default="end", help="links: end column (default: end)")
    p.add_argument("--links-gene-col", default="gene_name", help="links: gene symbol column (default: gene_name)")
    p.add_argument("--links-geneid-col", default="gene_id", help="links: Ensembl gene id column (default: gene_id)")

    p.add_argument("--output", help="output .tsv.gz path (default: ./epimap_open_chromatin.tsv.gz)")
    p.add_argument("--stage", action="store_true", help="upload .tsv.gz + .tbi to BOTH GCS buckets (OFF by default)")
    p.add_argument("--gcs-finngen", default=GCS_FINNGEN, help="finngen GCS destination")
    p.add_argument("--gcs-daly", default=GCS_DALY, help="daly GCS destination")

    p.add_argument("--sample", "--dry-run", dest="sample", action="store_true",
                   help="run the transform on tiny synthetic input, print first rows, validate tabix; no upload")
    return p.parse_args()


# "chr1:1000-2000", "chr1-1000-2000", "chr1_1000_2000"; "chr" prefix optional and normalized on.
_ID_RE = re.compile(r"^(?:chr)?([^:_\-]+)[:_\-](\d+)[_\-](\d+)$")
_BIOSAMPLE_RE = re.compile(r"(BSS\d+)")
_NUM_PREFIX_RE = re.compile(r"^[Ee]?\d+[_\-]")
_PURE_NUM_RE = re.compile(r"^[Ee]?(\d+)$")


def normalize_state(raw: str | None) -> str | None:
    """Map a raw ChromHMM state label to its canonical name, or None if unmappable.

    Prefers the NAME token when present (unambiguous) and only falls back to the numeric index,
    so e.g. "EnhA1", "9_EnhA1", "E9_EnhA1" all resolve by name while bare "E10"/"10" resolve to
    STATE_ORDER[10-1]. Handles the "ZNF/Rpts" spelling via _STATE_ALIASES.
    """
    if raw is None:
        return None
    tok = str(raw).strip()
    if not tok:
        return None

    # try the label as a name first (unambiguous), then with any leading "<num>_"/"E<num>_" stripped
    for cand in (tok.replace("/", "_"), _NUM_PREFIX_RE.sub("", tok.replace("/", "_"))):
        low = cand.lower()
        if low in _CANON_BY_LOWER:
            return _CANON_BY_LOWER[low]
        if low in _STATE_ALIASES:
            return _STATE_ALIASES[low]

    # fall back to a bare numeric / "E<num>" index into the canonical 18-state order
    m = _PURE_NUM_RE.match(tok)
    if m:
        idx = int(m.group(1))
        if 1 <= idx <= len(STATE_ORDER):
            return STATE_ORDER[idx - 1]
    return None


def normalize_life_stage(value: str | None) -> str:
    """adult by default; fetal only when the metadata value marks it fetal/embryonic."""
    if value is None:
        return "adult"
    low = str(value).strip().lower()
    if "fetal" in low or "embryo" in low:
        return "fetal"
    return "adult"


def biosample_from_path(path: str, override: str | None) -> str:
    if override:
        return override
    name = Path(path).name
    m = _BIOSAMPLE_RE.search(name)
    return m.group(1) if m else name.split(".")[0]


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


def load_calls(path: str, biosample: str, state_col: int) -> pl.DataFrame:
    """Read one biosample's ChromHMM CALLS BED, filter to included states, return long rows.

    Output columns: chrom, start, end, cell_type(=biosample), state (canonical name).
    Excluded/unmappable states are dropped here.
    """
    df = pl.read_csv(path, separator="\t", has_header=False, comment_prefix="#", infer_schema_length=0)
    cols = df.columns
    if len(cols) < max(4, state_col):
        raise ValueError(f"{path}: expected >= {max(4, state_col)} BED columns, got {len(cols)}")

    df = df.select(
        _numeric_chrom(pl.col(cols[0])).alias("chrom"),
        pl.col(cols[1]).cast(pl.Int64).alias("start"),
        pl.col(cols[2]).cast(pl.Int64).alias("end"),
        pl.col(cols[state_col - 1]).cast(pl.Utf8).alias("_state_raw"),
    )

    # 18 distinct states -> resolve once per unique label, then map the column
    mapping = {r: normalize_state(r) for r in df["_state_raw"].unique().to_list()}
    df = df.with_columns(pl.col("_state_raw").replace_strict(mapping, default=None).alias("state"))
    df = df.filter(pl.col("state").is_in(list(INCLUDED_STATES)))

    return df.select(
        "chrom", "start", "end",
        pl.lit(biosample).alias("cell_type"),
        pl.col("state"),
    )


def _coords_from_id(df: pl.DataFrame, id_col: str) -> pl.DataFrame:
    parsed = df.select(pl.col(id_col).str.extract_groups(_ID_RE.pattern).alias("_g")).unnest("_g")
    return df.with_columns(
        _numeric_chrom(parsed["1"]).alias("chrom"),
        parsed["2"].cast(pl.Int64).alias("start"),
        parsed["3"].cast(pl.Int64).alias("end"),
    )


def _normalize_link_coords(df: pl.DataFrame, chrom_col: str, start_col: str, end_col: str) -> pl.DataFrame:
    df = df.rename({chrom_col: "chrom", start_col: "start", end_col: "end"})
    return df.with_columns(
        _numeric_chrom(pl.col("chrom")).alias("chrom"),
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
    )


def _peak_key() -> pl.Expr:
    return pl.format("{}-{}-{}", pl.col("chrom"), pl.col("start"), pl.col("end"))


def load_gene_links(args: argparse.Namespace) -> pl.DataFrame:
    """Read enhancer->gene links and aggregate to one comma-joined gene list per interval key.

    Symbol and id are paired JOINTLY (dedup + sort the (symbol, id) structs together, never
    independently) so the i-th target_gene token matches the i-th target_gene_id token.
    """
    df = pl.read_csv(args.links, separator="\t", infer_schema_length=10_000)
    if args.links_id_col:
        df = _coords_from_id(df, args.links_id_col)
    else:
        df = _normalize_link_coords(df, args.links_chrom_col, args.links_start_col, args.links_end_col)

    df = df.with_columns(_peak_key().alias("_key"))
    gene = pl.col(args.links_gene_col).cast(pl.Utf8) if args.links_gene_col in df.columns else pl.lit(None, dtype=pl.Utf8)
    geneid = pl.col(args.links_geneid_col).cast(pl.Utf8) if args.links_geneid_col in df.columns else pl.lit(None, dtype=pl.Utf8)
    df = df.with_columns(gene.fill_null("").alias("_gene"), geneid.fill_null("").alias("_geneid"))

    return (
        df.group_by("_key")
        .agg(pl.struct("_gene", "_geneid").unique().sort().alias("_pairs"))
        .with_columns(
            pl.col("_pairs").list.eval(pl.element().struct.field("_gene")).list.join(",").alias("target_gene"),
            pl.col("_pairs").list.eval(pl.element().struct.field("_geneid")).list.join(",").alias("target_gene_id"),
        )
        .drop("_pairs")
    )


def load_biosample_metadata(args: argparse.Namespace) -> dict[str, tuple[str, str]]:
    """biosample -> (tissue, life_stage). tissue verbatim; life_stage normalized to adult/fetal."""
    if not args.biosample_metadata:
        return {}
    m = pl.read_csv(args.biosample_metadata, separator="\t", infer_schema_length=10_000)
    bcol, tcol, lcol = args.meta_biosample_col, args.meta_tissue_col, args.meta_lifestage_col
    if bcol not in m.columns or tcol not in m.columns:
        raise ValueError(f"--biosample-metadata missing '{bcol}' or '{tcol}' (have {m.columns})")
    has_life = lcol in m.columns
    out: dict[str, tuple[str, str]] = {}
    for row in m.iter_rows(named=True):
        bs = str(row[bcol])
        tissue = row[tcol]
        tissue = str(tissue) if tissue not in (None, "") else "unknown"
        life = normalize_life_stage(row[lcol]) if has_life else "adult"
        out[bs] = (tissue, life)
    return out


def build_output(long: pl.DataFrame, args: argparse.Namespace) -> pl.DataFrame:
    """Assemble the 18 canonical columns in order from the long (segment, biosample, state) table."""
    df = long.with_columns(_peak_key().alias("peak_id"))

    # per-biosample tissue + life_stage (+ derived uberon); default unknown/adult when absent
    biomap = load_biosample_metadata(args)
    biosamples = df["cell_type"].unique().to_list()
    rows = []
    for bs in biosamples:
        tissue, life = biomap.get(bs, ("unknown", "adult"))
        uberon = _TISSUE_UBERON.get(tissue.lower())
        rows.append((bs, tissue, life, uberon))
    meta = pl.DataFrame(
        rows,
        schema={"cell_type": pl.Utf8, "tissue": pl.Utf8, "life_stage": pl.Utf8, "uberon_id": pl.Utf8},
        orient="row",
    )
    df = df.join(meta, on="cell_type", how="left")

    if args.links:
        links = load_gene_links(args)
        df = df.with_columns(_peak_key().alias("_key")).join(links, on="_key", how="left").drop("_key")
    else:
        df = df.with_columns(
            pl.lit(None, dtype=pl.Utf8).alias("target_gene"),
            pl.lit(None, dtype=pl.Utf8).alias("target_gene_id"),
        )

    df = df.with_columns(
        pl.lit(DATASET).alias("dataset"),
        pl.lit(CONDITION).alias("condition"),
        pl.lit(ASSAY).alias("assay"),
        pl.lit(None, dtype=pl.Utf8).alias("score"),           # presence-based: score is empty
        ("chromhmm:" + pl.col("state")).alias("score_type"),  # exact state carried here
        pl.lit(None, dtype=pl.Utf8).alias("n_cells"),
        pl.lit(None, dtype=pl.Utf8).alias("cell_ontology_id"),
        pl.lit(VERSION).alias("version"),
    )

    return df.select(OUTPUT_COLUMNS)


def resolve_calls_files(args: argparse.Namespace) -> list[str]:
    files = list(args.calls)
    if args.calls_dir:
        files.extend(sorted(str(p) for p in Path(args.calls_dir).glob(args.calls_glob)))
    return files


def load_all_calls(files: list[str], args: argparse.Namespace) -> pl.DataFrame:
    frames = []
    for path in files:
        bs = biosample_from_path(path, args.cell_type if len(files) == 1 else None)
        frames.append(load_calls(path, bs, args.state_col))
    return pl.concat(frames) if frames else pl.DataFrame()


def write_open_chromatin(df: pl.DataFrame, output_path: str) -> None:
    """sort -k1,1 -k2,2n -> bgzip -> tabix -p bed (interval index), missing values as "NA"."""
    tmpdir = tempfile.mkdtemp()
    body = Path(tmpdir) / "body.tsv"
    # every empty/missing cell serialized as the literal "NA"
    df.write_csv(body, separator="\t", include_header=False, null_value="NA")

    header = "#" + "\t".join(OUTPUT_COLUMNS)
    pipeline = (
        f'( printf "%s\\n" {shell_quote(header)}; '
        f'LC_ALL=C sort -k1,1 -k2,2n {shell_quote(str(body))} ) | bgzip -c > {shell_quote(output_path)}'
    )
    subprocess.run(pipeline, shell=True, check=True, executable="/bin/bash")
    subprocess.run(["tabix", "-f", "-p", "bed", output_path], check=True)
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


def _synthetic_inputs(tmpdir: Path) -> tuple[list[str], str, str]:
    """Write tiny synthetic CALLS files (mixed states + label forms), metadata, and gene links.

    Exercises: multiple biosamples; INCLUDED vs EXCLUDED states; name / numbered-name / numeric /
    "ZNF/Rpts" label forms; metadata tissue + fetal life_stage; a biosample absent from metadata
    (-> unknown/adult); a multi-gene enhancer link (joint pairing).
    """
    # BSS00001: TssA(incl,name), Quies(excl,name), 9_EnhA1(incl,numbered->EnhA1), ZNF/Rpts(excl,alias),
    #           plus a chrX TssA(incl) to exercise the numeric chrX -> 23 mapping
    b1 = tmpdir / "BSS00001_18_CALLS_segments.bed"
    b1.write_text(
        "chr1\t778420\t779420\tTssA\t0\t.\t778420\t779420\t255,0,0\n"
        "chr1\t780420\t826620\tQuies\t0\t.\t780420\t826620\t255,255,255\n"
        "chr1\t903820\t904020\t9_EnhA1\t0\t.\t903820\t904020\t255,195,77\n"
        "chr2\t500000\t500200\tZNF/Rpts\t0\t.\t500000\t500200\t102,205,170\n"
        "chrX\t155000\t155300\tTssA\t0\t.\t155000\t155300\t255,0,0\n"
    )
    # BSS00099: EnhA2(incl,name), 5(excl,numeric->Tx), EnhBiv(incl,name), E18(excl,numeric->Quies)
    b2 = tmpdir / "BSS00099_18_CALLS_segments.bed"
    b2.write_text(
        "chr1\t778420\t779420\tEnhA2\t0\t.\t778420\t779420\t255,195,77\n"
        "chr1\t900000\t900500\t5\t0\t.\t900000\t900500\t0,128,0\n"
        "chr3\t123400\t123900\tEnhBiv\t0\t.\t123400\t123900\t189,183,107\n"
        "chr3\t200000\t201000\tE18\t0\t.\t200000\t201000\t255,255,255\n"
    )
    # BSS00500: TssBiv(incl) — biosample intentionally absent from metadata -> unknown/adult
    b3 = tmpdir / "BSS00500_18_CALLS_segments.bed"
    b3.write_text(
        "chr1\t904020\t904220\tTssBiv\t0\t.\t904020\t904220\t205,92,92\n"
    )

    meta = tmpdir / "biosample_metadata.tsv"
    meta.write_text(
        "biosample\ttissue\tlife_stage\n"
        "BSS00001\tbrain\tadult\n"
        "BSS00099\tblood\tembryonic\n"   # -> life_stage normalized to fetal
    )
    links = tmpdir / "enhancer_gene_links.tsv"
    links.write_text(
        "chrom\tstart\tend\tgene_name\tgene_id\n"
        "chr1\t903820\t904020\tSAMD11\tENSG00000187634\n"   # matches BSS00001 9_EnhA1 segment
        "chr1\t903820\t904020\tNOC2L\tENSG00000188976\n"
    )
    return [str(b1), str(b2), str(b3)], str(meta), str(links)


def run_sample() -> None:
    print("=== SAMPLE / DRY-RUN: synthetic input, no GCS upload ===")
    tmpdir = Path(tempfile.mkdtemp())
    calls, meta, links = _synthetic_inputs(tmpdir)

    args = parse_args()
    args.calls = calls
    args.calls_dir = None
    args.state_col = 4
    args.cell_type = None
    args.biosample_metadata = meta
    args.meta_biosample_col, args.meta_tissue_col, args.meta_lifestage_col = "biosample", "tissue", "life_stage"
    args.links = links
    args.links_id_col = None
    args.links_chrom_col, args.links_start_col, args.links_end_col = "chrom", "start", "end"
    args.links_gene_col, args.links_geneid_col = "gene_name", "gene_id"

    # state normalization: numeric + name + alias forms
    assert normalize_state("TssA") == "TssA"
    assert normalize_state("9_EnhA1") == "EnhA1"      # numbered-name resolves by name
    assert normalize_state("E9") == "EnhA1"           # bare "E<num>" resolves by index
    assert normalize_state("10") == "EnhA2"           # bare numeric resolves by index
    assert normalize_state("ZNF/Rpts") == "ZNF_Rpts"  # spelling alias
    assert normalize_state("18") == "Quies"
    print("  state normalization (name/numbered/numeric/alias): OK")

    long = load_all_calls(calls, args)
    out = build_output(long, args)

    assert out.columns == OUTPUT_COLUMNS, f"column order mismatch: {out.columns}"
    assert len(out.columns) == 18
    print(f"  output has {len(out.columns)} columns in canonical order: OK")

    # excluded states dropped; only included states survive, each as chromhmm:<State>
    kept_states = sorted({st.removeprefix("chromhmm:") for st in out["score_type"].to_list()})
    assert all(s in INCLUDED_STATES for s in kept_states), kept_states
    assert set(kept_states) == {"TssA", "EnhA1", "EnhA2", "EnhBiv", "TssBiv"}, kept_states
    assert out["score_type"].str.starts_with("chromhmm:").all()
    assert out["score"].is_null().all(), "score must be empty (presence-based)"
    assert out["assay"].eq("chromHMM").all()
    print(f"  excluded states dropped; kept (chromhmm:<State>): {kept_states}")
    print("  score empty (presence-based), assay=chromHMM: OK")

    # life_stage: BSS00099 marked embryonic -> fetal; BSS00500 absent -> unknown/adult
    ls = {r["cell_type"]: (r["tissue"], r["life_stage"]) for r in out.select("cell_type", "tissue", "life_stage").unique().iter_rows(named=True)}
    assert ls["BSS00001"] == ("brain", "adult"), ls
    assert ls["BSS00099"] == ("blood", "fetal"), ls
    assert ls["BSS00500"] == ("unknown", "adult"), ls
    print(f"  tissue/life_stage derivation: OK  {ls}")

    # chrom / peak_id are numeric with chrX -> 23
    chroms = set(out["chrom"].to_list())
    assert chroms <= {"1", "2", "3", "23"}, chroms
    assert "23" in chroms, f"chrX should map to 23; got {chroms}"
    xrow = out.filter(pl.col("chrom") == "23").row(0, named=True)
    assert xrow["peak_id"] == "23-155000-155300", xrow["peak_id"]
    print(f"  numeric chrom (chrX -> 23): OK  chroms={sorted(chroms)}  peak_id={xrow['peak_id']}")

    # multi-gene enhancer link joined with joint symbol/id pairing
    linked = out.filter(pl.col("peak_id") == "1-903820-904020").row(0, named=True)
    assert linked["target_gene"] == "NOC2L,SAMD11", linked["target_gene"]
    assert linked["target_gene_id"] == "ENSG00000188976,ENSG00000187634", linked["target_gene_id"]
    pairs = dict(zip(linked["target_gene"].split(","), linked["target_gene_id"].split(",")))
    assert pairs["SAMD11"] == "ENSG00000187634" and pairs["NOC2L"] == "ENSG00000188976", pairs
    print(f"  enhancer->gene joint pairing: OK  {pairs}")

    print("\nfirst output rows:")
    with pl.Config(tbl_cols=-1, tbl_width_chars=260, fmt_str_lengths=40):
        print(out.head(10))

    out_path = str(tmpdir / f"{DATASET}.sample.tsv.gz")
    write_open_chromatin(out, out_path)

    # empty cells serialize as the literal "NA" (score is presence-based -> NA); 18 columns per row
    import gzip
    with gzip.open(out_path, "rt") as fh:
        body_lines = [ln for ln in fh if not ln.startswith("#")]
    assert all(len(ln.rstrip("\n").split("\t")) == 18 for ln in body_lines), "not 18 columns"
    assert any("\tNA\t" in ln for ln in body_lines), "no NA emitted"
    print(f"  empty cells rendered as NA, 18 columns per row: OK ({len(body_lines)} rows)")

    q = subprocess.run(["tabix", out_path, "1:903900-903950"], capture_output=True, text=True, check=True)
    print("\ntabix overlap query 1:903900-903950 ->")
    print(q.stdout.rstrip() or "  (no rows)")
    assert q.stdout.strip(), "interval query returned nothing — indexing is wrong"
    qx = subprocess.run(["tabix", out_path, "23:155100-155200"], capture_output=True, text=True, check=True)
    print("tabix overlap query 23:155100-155200 (chrX) ->")
    print(qx.stdout.rstrip() or "  (no rows)")
    assert qx.stdout.strip(), "chrX (23) interval query returned nothing — indexing is wrong"
    print("\n=== SAMPLE OK (no upload performed) ===")


def main() -> None:
    args = parse_args()
    if args.sample:
        run_sample()
        return

    files = resolve_calls_files(args)
    if not files:
        raise SystemExit("no CALLS files (pass paths and/or --calls-dir, or use --sample). See --help.")

    output_path = args.output or f"./{DATASET}.tsv.gz"
    print(f"Reading {len(files)} biosample CALLS file(s) ...")
    long = load_all_calls(files, args)
    print(f"  {long.height} accessible/active segments across {long['cell_type'].n_unique()} biosamples")

    out = build_output(long, args)
    assert out.columns == OUTPUT_COLUMNS, f"column order mismatch: {out.columns}"

    write_open_chromatin(out, output_path)

    if args.stage:
        print("Staging to GCS (both buckets) ...")
        upload_to_gcs(output_path, args.gcs_finngen)
        upload_to_gcs(output_path, args.gcs_daly)
    else:
        print("  --stage not set: skipping GCS upload")
    print("Done.")


if __name__ == "__main__":
    main()
