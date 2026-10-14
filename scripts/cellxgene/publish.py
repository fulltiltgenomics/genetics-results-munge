#!/usr/bin/env python3
"""CELLxGENE publish step: the cross-dataset reduce from per-dataset worker outputs to a release.

    publish.py --agg <dir | gs://.../agg/v<FORMAT>> --out <dir | gs://.../releases>
               [--release-date YYYY-MM-DD] [--listing datasets.json] [--genes CD3E,ALB,...]
               [--marker-panel] [--gencode TSV] [--workdir DIR]

Inputs
- `--agg`: one subdirectory per dataset_version_id holding `aggregate.py`'s outputs. Only
  subdirectories with `_SUCCESS.json` are read; any other is listed in the manifest as skipped.
  A `_SUCCESS.json` whose `aggregation_format_version` is not `aggregate.FORMAT_VERSION` stops
  the run: a release never mixes formats. A gs:// prefix is copied to `--workdir` first
- roll-up tables `rollup/cell_type_rollup.tsv` and `rollup/tissue_rollup.tsv` (pinned, built by
  `rollup.py`). A native term missing from a table maps to itself and is counted in the manifest
- GENCODE v49 genes (`--gencode`): Ensembl ids are matched with their `.N` version stripped on
  both sides; the source's own id is kept as `source_gene_id`. Features with no v49 match are
  dropped and counted in the manifest (`gencode`), as are non-Ensembl features (spike-ins)
- `--listing`: the curation API datasets JSON, for the collection, title, DOI and citation in
  `datasets.tsv`. Without it those columns are empty

Rules
- PRIMARY metrics, per cell group: `mean_log1p_cp10k` (cell-level mean, zeros included) and
  `frac_expressing` (n_expressing / n_cells); these are what rank. SECONDARY:
  `donor_pseudobulk_mean_cp10k`, the donor-pseudobulk mean: per donor 1e4 * summed counts /
  summed library size, averaged over the donors with a nonzero library
  (`donor_pseudobulk_n_donors`), and NA unless that count is >= `--min-pseudobulk-donors`
  (default 3). It is not on the log scale and does not rank anything
- NA versus 0 (see `aggregate.expression_matrix`): a gene measured by a dataset but with no
  `group_gene` row is 0; a gene the dataset did not measure contributes nothing. A pooled
  value's denominators (n_cells, n_donors, n_datasets, pseudobulk donors) are summed over the
  datasets that MEASURED the gene, so an unmeasured dataset never dilutes a mean towards 0
- suppression: a row (per-dataset group or pooled summary row) is emitted only with
  n_cells >= `--min-cells` (50) and n_donors >= `--min-donors` (3), counted per gene among
  measuring datasets. Rows below are dropped, not shown as NA
- summary: groups with disease `PATO:0000461` (normal), tissue_type `tissue`, cell type not
  `unknown`, and not a spatial assay (a Visium spot is not a cell; `--summary-include-spatial`
  admits them), pooled across datasets at (rolled-up cell type, tissue_general). n_donors is
  the sum over datasets of distinct donor_ids, so a donor published in two datasets counts twice
- tau (Yanai 2005) per gene on the cell-type profile: mean log1p CP10k pooled over all summary
  tissues per rolled-up cell type, over the cell types that pass suppression at that level;
  tau = sum(1 - x_i / max x) / (n - 1). NA when n < 2 or max x = 0. The same value is repeated
  on every summary row of the gene
- per-dataset grain: every group of every dataset (any disease, tissue_type, assay), with the
  rolled-up cell type and tissue_general alongside the native terms. SPARSE: only genes with
  n_expressing >= 1 in the group have a row, as in the worker's output. A gene with no row is 0
  if `measured_genes.tsv.gz` lists it for that dataset, else not measured (NA)

Outputs in `<out>/<release-date>/`, built in a fresh directory and moved into place (replacing
any earlier run of that date); gs:// is built in `--workdir` and uploaded, manifest last. A
`--genes` run writes to `<out>/<release-date>/genes-mode/` instead and refuses a gs:// `--out`
unless `--force-remote-genes`
- `summary.tsv.gz`, `per_dataset.tsv.gz`: bgzipped, sorted by (chrom, gene_start), point-indexed
  `tabix -s2 -b3 -e3` on gene_start like GTEx. chrom is the GENCODE file's own numeric coding
- `bq/summary.tsv.gz`, `bq/per_dataset.tsv.gz`, `bq/groups.tsv.gz`, `bq/datasets.tsv`: BigQuery
  load files (gzip, header row, empty = NULL). The per-dataset grain carries integer
  `dataset_key` and `group_key` into `groups` and `datasets` instead of repeating the labels
- `measured_genes.tsv.gz`: (dataset_version_id, v49 gene_id) for every gene a dataset measured
  (`measured` in its genes.parquet) that maps to v49; read from the worker output, no H5AD
- `manifest.json`: release, inputs, roll-up checksums, gencode drop counts, thresholds, row
  counts, and every other file's bytes and sha256. Each output's column list is the `*_COLUMNS` constant below

Memory: one chromosome at a time. A full run first splits each dataset's `group_gene` and
`group_donor_gene` by chromosome into `--workdir` (one streaming pass per file); each chromosome
then reads only its own pieces, reduced one dataset at a time. `--genes` skips the split and
filters the inputs directly, which is what makes the marker panel take seconds.
"""
import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import aggregate  # noqa: E402
from peak_utils import TabixIndex  # noqa: E402

RELEASE_FORMAT_VERSION = 1
NORMAL = "PATO:0000461"
MIN_CELLS, MIN_DONORS, MIN_PB_DONORS = 50, 3, 3
GENCODE_URI = "gs://daly-genetics-results/mapping_files/gencode.v49.annotation.genes.tsv"
GENE_POINT_INDEX = TabixIndex(["-s", "2", "-b", "3", "-e", "3"], "POINT: tabix -s2 -b3 -e3 (gene_start)")
LABELS = ["tissue", "cell_type", "disease", "assay"]

GENE_COLUMNS = ["gene_id", "chrom", "gene_start", "gene_end", "gene_name"]
SUMMARY_COLUMNS = GENE_COLUMNS + [
    "source_gene_ids", "cell_type_rollup_ontology_term_id", "cell_type_rollup",
    "tissue_general_ontology_term_id", "tissue_general", "n_cells", "n_donors", "n_datasets",
    "mean_log1p_cp10k", "frac_expressing", "donor_pseudobulk_mean_cp10k", "donor_pseudobulk_n_donors", "tau"]
GROUP_COLUMNS = [
    "dataset_id", "tissue_ontology_term_id", "tissue", "tissue_type",
    "tissue_general_ontology_term_id", "tissue_general", "cell_type_ontology_term_id", "cell_type",
    "cell_type_rollup_ontology_term_id", "cell_type_rollup", "disease_ontology_term_id", "disease",
    "assay_ontology_term_id", "assay", "suspension_type", "spatial", "n_cells", "n_donors"]
METRIC_COLUMNS = ["mean_log1p_cp10k", "frac_expressing", "donor_pseudobulk_mean_cp10k"]
PER_DATASET_COLUMNS = GENE_COLUMNS + ["source_gene_id"] + GROUP_COLUMNS + METRIC_COLUMNS
BQ_PER_DATASET_COLUMNS = ["chrom", "gene_start", "gene_id", "gene_name", "dataset_key", "group_key",
                          *METRIC_COLUMNS]
BQ_GROUP_COLUMNS = ["group_key", "dataset_key", "dataset_version_id", *GROUP_COLUMNS]
DATASET_COLUMNS = ["dataset_key", "dataset_id", "dataset_version_id", "collection_id",
                   "collection_name", "title", "collection_doi", "citation", "cells_kept",
                   "n_groups", "n_genes", "n_genes_measured", "n_genes_mapped"]


# ---------------------------------------------------------------- inputs

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# GCS goes through google-cloud-storage, never the gcloud CLI: the worker image the driver runs
# this in has the library and no gcloud binary
def localize_file(uri, dest):
    """A gs:// object downloaded into dest; a local path returned as is."""
    if not uri.startswith("gs://"):
        return uri
    os.makedirs(dest, exist_ok=True)
    bucket, key = aggregate._gcs_bucket(uri)
    path = os.path.join(dest, key.rsplit("/", 1)[-1])
    bucket.blob(key).download_to_filename(path)
    return path


def localize_agg(uri, names, dest):
    """Only the named dataset_version_ids under a gs:// agg prefix, so superseded or withdrawn
    leftovers in the prefix are never read; a local dir returned as is."""
    if not uri.startswith("gs://"):
        return uri
    bucket, prefix = aggregate._gcs_bucket(uri)
    for n in sorted(names):
        for blob in bucket.list_blobs(prefix=f"{prefix}/{n}/"):
            rel = blob.name[len(prefix) + 1:]
            os.makedirs(os.path.dirname(os.path.join(dest, rel)), exist_ok=True)
            blob.download_to_filename(os.path.join(dest, rel))
    os.makedirs(dest, exist_ok=True)
    return dest


def upload(local, final, files):
    """files first, manifest.json last: a release without its manifest is incomplete."""
    bucket, prefix = aggregate._gcs_bucket(final)
    for f in [f for f in files if f != "manifest.json"] + ["manifest.json"]:
        bucket.blob(f"{prefix}/{f}").upload_from_filename(os.path.join(local, f))


def strip_version(col):
    # only the version: a `_PAR_Y` suffix after it names a different gene
    return pl.col(col).str.replace(r"\.\d+", "")


def load_gencode(path):
    g = pl.read_csv(path, separator="\t", infer_schema_length=0).select(
        pl.col("gene_id"), pl.col("chrom").cast(pl.Int64), pl.col("gene_start").cast(pl.Int64),
        pl.col("gene_end").cast(pl.Int64), pl.col("gene_name"))
    g = g.with_columns(ensembl=strip_version("gene_id")).unique("ensembl", keep="first", maintain_order=True)
    return g.sort("chrom", "gene_start", "gene_id").with_row_index("gidx").with_columns(pl.col("gidx").cast(pl.Int32))


def load_rollups(rollup_dir):
    ct = pl.read_csv(f"{rollup_dir}/cell_type_rollup.tsv", separator="\t", infer_schema_length=0).select(
        "cell_type_ontology_term_id", pl.col("rollup_term_id").alias("ct_r_id"),
        pl.col("rollup_label").alias("ct_r_label"))
    ti = pl.read_csv(f"{rollup_dir}/tissue_rollup.tsv", separator="\t", infer_schema_length=0).select(
        "tissue_ontology_term_id", pl.col("tissue_general_ontology_term_id").alias("tg_id"),
        pl.col("tissue_general").alias("tg_label"))
    return ct, ti


def discover(agg_dir, only=None):
    names = sorted(n for n in os.listdir(agg_dir) if os.path.isdir(os.path.join(agg_dir, n)))
    if only:
        names = [n for n in names if n in only]
    done = [n for n in names if os.path.exists(os.path.join(agg_dir, n, "_SUCCESS.json"))]
    for n in done:
        with open(os.path.join(agg_dir, n, "_SUCCESS.json")) as fh:
            v = json.load(fh).get("aggregation_format_version")
        if v != aggregate.FORMAT_VERSION:
            raise SystemExit(f"{n}: aggregation_format_version {v} != {aggregate.FORMAT_VERSION}; "
                             "a release does not mix formats")
    return done, sorted(set(names) - set(done))


class Inputs:
    """Everything small, for all datasets: groups, donors, gene mapping, summary-row keys."""

    def __init__(self, agg_dir, names, gencode, ct, ti, a):
        self.agg_dir, self.names, self.gencode = agg_dir, names, gencode
        self.success, groups, gds, maps, self.gene_stats = [], [], [], [], []
        v49 = gencode.select("ensembl", "gidx", "chrom")
        offset = 0
        for d, n in enumerate(names):
            p = os.path.join(agg_dir, n)
            with open(f"{p}/_SUCCESS.json") as fh:
                self.success.append(json.load(fh))
            g = pl.read_parquet(f"{p}/groups.parquet")
            for c in LABELS:
                if c not in g.columns:
                    g = g.with_columns(pl.lit(None, pl.Utf8).alias(c))
            if "spatial" not in g.columns:
                g = g.with_columns(spatial=pl.col("assay_ontology_term_id").is_in(list(aggregate.SPATIAL_ASSAYS)))
            g = g.with_columns(pl.col("group").cast(pl.Int32), pl.col("n_cells").cast(pl.Int64),
                               pl.col("n_donors").cast(pl.Int64), dataset_key=pl.lit(d, pl.Int32),
                               group_key=(pl.col("group") + offset).cast(pl.Int64),
                               dataset_id=pl.lit(self.success[-1].get("dataset_id")),
                               dataset_version_id=pl.lit(n))
            # group ids can be gapped and need not start at 0 (one worker output has 16..31), so
            # the next dataset starts past this one's largest id, not its row count
            offset += (g["group"].max() + 1) if g.height else 0
            groups.append(g)
            gd = pl.read_parquet(f"{p}/group_donor.parquet").select(
                pl.col("gd").cast(pl.Int32), pl.col("group").cast(pl.Int32), pl.col("donor_id").cast(pl.Utf8),
                pl.col("library_size").cast(pl.Float64))
            gds.append(gd.with_columns(dataset_key=pl.lit(d, pl.Int32)))
            genes = pl.read_parquet(f"{p}/genes.parquet").select(
                pl.col("gene").cast(pl.Int32), pl.col("gene_id").alias("source_gene_id"), "measured")
            genes = genes.with_columns(ensembl=strip_version("source_gene_id"))
            meas = genes.filter(pl.col("measured"))
            ens = meas.filter(pl.col("ensembl").str.starts_with("ENSG"))
            m = ens.join(v49, on="ensembl", how="inner").unique("gidx", keep="first", maintain_order=True)
            self.gene_stats.append({"n_genes": genes.height, "n_genes_measured": meas.height,
                                    "non_ensembl": meas.height - ens.height,
                                    "unmatched": sorted(set(ens["ensembl"]) - set(v49["ensembl"])),
                                    "n_genes_mapped": m.height})
            maps.append(m.select("gene", "gidx", "chrom", "source_gene_id").with_columns(dataset_key=pl.lit(d, pl.Int32)))
        g = pl.concat(groups, how="diagonal_relaxed")
        if g["group_key"].n_unique() != g.height:
            raise SystemExit("duplicate group_key across datasets")
        g = g.join(ct, on="cell_type_ontology_term_id", how="left").join(ti, on="tissue_ontology_term_id", how="left")
        self.unmapped = {"cell_type": sorted(set(g.filter(pl.col("ct_r_id").is_null())["cell_type_ontology_term_id"])),
                         "tissue": sorted(set(g.filter(pl.col("tg_id").is_null())["tissue_ontology_term_id"]))}
        g = g.with_columns(
            cell_type_rollup_ontology_term_id=pl.coalesce("ct_r_id", "cell_type_ontology_term_id"),
            cell_type_rollup=pl.coalesce("ct_r_label", "cell_type"),
            tissue_general_ontology_term_id=pl.coalesce("tg_id", "tissue_ontology_term_id"),
            tissue_general=pl.coalesce("tg_label", "tissue"))
        gd = pl.concat(gds)
        n_pbd = gd.filter(pl.col("library_size") > 0).group_by("dataset_key", "group").agg(n_pbd=pl.len())
        g = g.join(n_pbd, on=["dataset_key", "group"], how="left").with_columns(pl.col("n_pbd").fill_null(0))
        g = g.with_columns(
            emit=(pl.col("n_cells") >= a.min_cells) & (pl.col("n_donors") >= a.min_donors),
            summary=(pl.col("disease_ontology_term_id") == NORMAL) & (pl.col("tissue_type") == "tissue")
            & (pl.col("cell_type_ontology_term_id") != "unknown")
            & (pl.lit(a.summary_include_spatial) | ~pl.col("spatial")))
        # summary row = (rolled-up cell type, tissue_general); ct = the cell-type level used for tau
        el = g.filter("summary")
        rows = el.group_by("cell_type_rollup_ontology_term_id", "tissue_general_ontology_term_id").agg(
            pl.col("cell_type_rollup").sort().first(), pl.col("tissue_general").sort().first()).sort(
            "cell_type_rollup_ontology_term_id", "tissue_general_ontology_term_id").with_row_index("row")
        cts = rows.select("cell_type_rollup_ontology_term_id").unique().sort(
            "cell_type_rollup_ontology_term_id").with_row_index("ct")
        self.rows = rows.with_columns(pl.col("row").cast(pl.Int32)).join(
            cts.with_columns(pl.col("ct").cast(pl.Int32)), on="cell_type_rollup_ontology_term_id")
        g = g.join(self.rows.select("cell_type_rollup_ontology_term_id", "tissue_general_ontology_term_id", "row", "ct"),
                   on=["cell_type_rollup_ontology_term_id", "tissue_general_ontology_term_id"], how="left")
        g = g.with_columns(pl.when("summary").then(pl.col("row")).alias("row"),
                           pl.when("summary").then(pl.col("ct")).alias("ct"))
        self.groups, self.gd, self.map = g, gd, pl.concat(maps)
        # per (dataset, level key): cells, distinct donors, donor libraries for the pseudobulk
        gdl = gd.join(g.filter("summary").select("dataset_key", "group", "row", "ct", "n_cells"),
                      on=["dataset_key", "group"])
        self.level = {}
        for key in ("row", "ct"):
            cells = g.filter("summary").group_by("dataset_key", key).agg(n_cells=pl.col("n_cells").sum())
            lib = gdl.group_by("dataset_key", key, "donor_id").agg(lib=pl.col("library_size").sum())
            don = lib.group_by("dataset_key", key).agg(n_donors=pl.len(), n_pbd=(pl.col("lib") > 0).sum())
            self.level[key] = (cells.join(don, on=["dataset_key", key], how="left").fill_null(0), lib)

    def path(self, d, name):
        return os.path.join(self.agg_dir, self.names[d], name)


# ---------------------------------------------------------------- per-chromosome split

KINDS = {"group_gene": ["group", "gene", "n_expressing", "sum_log1p_cp10k"],
         "group_donor_gene": ["gd", "gene", "counts"]}


def split_by_chrom(inp, stage):
    """One streaming pass per input file: rows of mapped genes, one Parquet file per chromosome."""
    for d in range(len(inp.names)):
        m = inp.map.filter(pl.col("dataset_key") == d)
        n_var = inp.gene_stats[d]["n_genes"]
        chrom_of = np.full(n_var, -1, np.int64)
        chrom_of[m["gene"].to_numpy()] = m["chrom"].to_numpy()
        for kind, cols in KINDS.items():
            writers = {}
            try:
                for b in pq.ParquetFile(inp.path(d, f"{kind}.parquet")).iter_batches(1 << 22, columns=cols):
                    ch = chrom_of[b.column("gene").to_numpy()]
                    for c in np.unique(ch[ch >= 0]):
                        sub = b.filter(pa.array(ch == c))
                        if c not in writers:
                            os.makedirs(f"{stage}/{kind}/{c}", exist_ok=True)
                            writers[c] = pq.ParquetWriter(f"{stage}/{kind}/{c}/{d}.parquet", sub.schema)
                        writers[c].write_batch(sub)
            finally:
                for w in writers.values():
                    w.close()


def read_part(inp, d, kind, chrom, genes, stage):
    if stage is None:
        return pl.scan_parquet(inp.path(d, f"{kind}.parquet")).select(KINDS[kind]).filter(
            pl.col("gene").is_in(genes.to_list())).collect()
    p = f"{stage}/{kind}/{chrom}/{d}.parquet"
    if not os.path.exists(p):
        return pl.DataFrame(schema={c: (pl.Float64 if c in ("sum_log1p_cp10k", "counts") else pl.Int32)
                                    for c in KINDS[kind]})
    return pl.read_parquet(p)


# ---------------------------------------------------------------- reduce

def pooled(num, pbn, measured, inp, key, a):
    """Pool one level (row or ct) across the datasets that measured each gene."""
    cells, _ = inp.level[key]
    den = measured.join(cells, on="dataset_key").group_by(key, "gidx").agg(
        pl.col("n_cells").sum(), pl.col("n_donors").sum(), pl.col("n_pbd").sum(),
        n_datasets=pl.col("dataset_key").n_unique())
    num = num.group_by(key, "gidx").agg(pl.col("sum_log").sum(), pl.col("n_expr").sum())
    pbn = pbn.group_by(key, "gidx").agg(pl.col("pb_sum").sum())
    s = den.join(num, on=[key, "gidx"], how="left").join(pbn, on=[key, "gidx"], how="left").with_columns(
        pl.col("sum_log", "n_expr", "pb_sum").fill_null(0))
    return s.with_columns(
        mean_log1p_cp10k=pl.col("sum_log") / pl.col("n_cells"),
        frac_expressing=pl.col("n_expr") / pl.col("n_cells"),
        donor_pseudobulk_mean_cp10k=pl.when(pl.col("n_pbd") >= a.min_pseudobulk_donors).then(pl.col("pb_sum") / pl.col("n_pbd")),
        pass_=(pl.col("n_cells") >= a.min_cells) & (pl.col("n_donors") >= a.min_donors))


def tau_of(ct):
    x = ct.filter("pass_").group_by("gidx").agg(
        n=pl.len(), mx=pl.col("mean_log1p_cp10k").max(),
        s=(1 - pl.col("mean_log1p_cp10k") / pl.col("mean_log1p_cp10k").max()).sum())
    return x.select("gidx", tau=pl.when((pl.col("n") >= 2) & (pl.col("mx") > 0)).then(pl.col("s") / (pl.col("n") - 1)))


def reduce_chrom(inp, chrom, a, stage, gidx_filter=None):
    """(summary rows, per-dataset rows, cell-type profile) for one chromosome."""
    groups = inp.groups
    num = {"row": [], "ct": []}
    pbn = {"row": [], "ct": []}
    per_ds, measured = [], []
    for d in range(len(inp.names)):
        m = inp.map.filter((pl.col("dataset_key") == d) & (pl.col("chrom") == chrom))
        if gidx_filter is not None:
            m = m.filter(pl.col("gidx").is_in(gidx_filter.to_list()))
        if m.is_empty():
            continue
        measured.append(m.select("dataset_key", "gidx", "source_gene_id"))
        genes = m["gene"]
        gm = m.select("gene", "gidx")
        gg = read_part(inp, d, "group_gene", chrom, genes, stage).join(gm, on="gene")
        gdg = read_part(inp, d, "group_donor_gene", chrom, genes, stage).join(gm, on="gene").join(
            inp.gd.filter(pl.col("dataset_key") == d).drop("dataset_key"), on="gd")
        gd_ = groups.filter(pl.col("dataset_key") == d)

        cp = gdg.filter(pl.col("library_size") > 0).with_columns(cp=1e4 * pl.col("counts") / pl.col("library_size"))
        pbg = cp.group_by("group", "gidx").agg(pb_sum=pl.col("cp").sum())
        em = gd_.filter("emit").select("group", "group_key", "n_cells", "n_pbd")
        per_ds.append(gg.join(em, on="group").join(pbg, on=["group", "gidx"], how="left").join(
            m.select("gidx", "source_gene_id"), on="gidx").select(
            "gidx", "group_key", "source_gene_id",
            mean_log1p_cp10k=pl.col("sum_log1p_cp10k") / pl.col("n_cells"),
            frac_expressing=pl.col("n_expressing") / pl.col("n_cells"),
            donor_pseudobulk_mean_cp10k=pl.when(pl.col("n_pbd") >= a.min_pseudobulk_donors).then(
                pl.col("pb_sum").fill_null(0) / pl.col("n_pbd"))))

        sm = gd_.filter("summary").select("group", "row", "ct")
        if sm.is_empty():
            continue
        ggs = gg.join(sm, on="group")
        gds = gdg.join(sm, on="group")
        for key in ("row", "ct"):
            num[key].append(ggs.group_by(key, "gidx").agg(
                sum_log=pl.col("sum_log1p_cp10k").sum(), n_expr=pl.col("n_expressing").sum()).with_columns(
                dataset_key=pl.lit(d, pl.Int32)))
            lib = inp.level[key][1].filter(pl.col("dataset_key") == d).drop("dataset_key")
            pbn[key].append(gds.group_by(key, "donor_id", "gidx").agg(c=pl.col("counts").sum()).join(
                lib, on=[key, "donor_id"]).filter(pl.col("lib") > 0).group_by(key, "gidx").agg(
                pb_sum=(1e4 * pl.col("c") / pl.col("lib")).sum()))
    if not measured:
        return None
    measured = pl.concat(measured)
    out = {}
    for key in ("row", "ct"):
        nm = pl.concat(num[key], how="diagonal_relaxed") if num[key] else pl.DataFrame(
            schema={key: pl.Int32, "gidx": pl.Int32, "sum_log": pl.Float64, "n_expr": pl.Int64})
        pb = pl.concat(pbn[key], how="diagonal_relaxed") if pbn[key] else pl.DataFrame(
            schema={key: pl.Int32, "gidx": pl.Int32, "pb_sum": pl.Float64})
        out[key] = pooled(nm, pb, measured.select("dataset_key", "gidx"), inp, key, a)
    gene = inp.gencode.select("gidx", *GENE_COLUMNS)
    src = measured.group_by("gidx").agg(source_gene_ids=pl.col("source_gene_id").unique().sort().str.join(","))
    ct = out["ct"].join(inp.rows.select("ct", "cell_type_rollup_ontology_term_id", "cell_type_rollup").unique("ct"), on="ct")
    summary = out["row"].filter("pass_").join(inp.rows, on="row").join(tau_of(out["ct"]), on="gidx", how="left").join(
        gene, on="gidx").join(src, on="gidx").rename({"n_pbd": "donor_pseudobulk_n_donors"}).sort(
        "gene_start", "gene_id", "cell_type_rollup_ontology_term_id", "tissue_general_ontology_term_id")
    pdr = pl.concat(per_ds, how="vertical_relaxed").join(gene, on="gidx").join(groups.select(
        "group_key", "dataset_key", *GROUP_COLUMNS), on="group_key").sort("gene_start", "gene_id", "group_key")
    return summary, pdr, ct.join(gene, on="gidx")


def fmt(df):
    return df.with_columns(pl.col("mean_log1p_cp10k", "frac_expressing").round(4),
                           pl.col("donor_pseudobulk_mean_cp10k").round_sig_figs(4),
                           *([pl.col("tau").round(4)] if "tau" in df.columns else []))


# ---------------------------------------------------------------- writing

class Pipe:
    """A compressor subprocess writing `path`, fed rows one chromosome at a time.

    Its own pipe rather than `peak_utils`' writers: a full release's per-dataset rows never fit
    one DataFrame, and an uncompressed body on disk would be several times the output.
    """

    def __init__(self, cmd, path, header):
        self.fh = open(path, "wb")
        self.proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=self.fh)
        self.proc.stdin.write(header.encode())

    def close(self):
        self.proc.stdin.close()
        rc = self.proc.wait()
        self.fh.close()
        if rc != 0:
            raise subprocess.CalledProcessError(rc, self.proc.args)


# gzip -n: no name or mtime in the header, so re-running the same inputs is byte-identical.
# BigQuery's gzip support is single-member gzip; bgzip's multi-member BGZF is kept for tabix
GZIP = ["gzip", "-n", "-6", "-c"]
BGZIP = ["bgzip", "-@", "4", "-c"]


class Release:
    def __init__(self, local, a):
        self.local = local
        os.makedirs(f"{local}/bq", exist_ok=True)
        self.cols = {"summary": (SUMMARY_COLUMNS, SUMMARY_COLUMNS), "per_dataset": (PER_DATASET_COLUMNS, BQ_PER_DATASET_COLUMNS)}
        self.tab = {k: Pipe(BGZIP, f"{local}/{k}.tsv.gz", "#" + "\t".join(c[0]) + "\n") for k, c in self.cols.items()}
        self.bq = {k: Pipe(GZIP, f"{local}/bq/{k}.tsv.gz", "\t".join(c[1]) + "\n") for k, c in self.cols.items()}
        self.rows = {k: 0 for k in self.cols}

    def append(self, kind, df):
        tab, bqc = self.cols[kind]
        df = fmt(df)
        df.select(tab).write_csv(self.tab[kind].proc.stdin, separator="\t", include_header=False, null_value="NA")
        df.select(bqc).write_csv(self.bq[kind].proc.stdin, separator="\t", include_header=False, null_value="")
        self.rows[kind] += df.height

    def close(self):
        for k in self.cols:
            self.tab[k].close()
            self.bq[k].close()
            subprocess.run(["tabix", "-f"] + GENE_POINT_INDEX.args + [f"{self.local}/{k}.tsv.gz"], check=True)


def datasets_table(inp, listing):
    meta = {x["dataset_version_id"]: x for x in listing or []}
    rows = []
    for d, n in enumerate(inp.names):
        s, m, gs = inp.success[d], meta.get(n, {}), inp.gene_stats[d]
        rows.append({"dataset_key": d, "dataset_id": s.get("dataset_id"), "dataset_version_id": n,
                     "collection_id": m.get("collection_id"), "collection_name": m.get("collection_name"),
                     "title": m.get("title"), "collection_doi": m.get("collection_doi"),
                     "citation": m.get("citation"), "cells_kept": s.get("cells_kept"),
                     "n_groups": s.get("n_groups"), "n_genes": gs["n_genes"],
                     "n_genes_measured": gs["n_genes_measured"], "n_genes_mapped": gs["n_genes_mapped"]})
    return pl.DataFrame(rows, schema_overrides={c: pl.Utf8 for c in DATASET_COLUMNS[1:8]}).select(DATASET_COLUMNS)


def resolve_genes(gencode, genes):
    want = [g.strip() for g in genes if g.strip()]
    ids = set(want)
    hit = gencode.filter(pl.col("gene_name").is_in(ids) | pl.col("ensembl").is_in(ids) | pl.col("gene_id").is_in(ids))
    missing = sorted(ids - set(hit["gene_name"]) - set(hit["ensembl"]) - set(hit["gene_id"]))
    return hit["gidx"], missing


def run(a):
    workdir = a.workdir or tempfile.mkdtemp(prefix="cxg-publish-")
    os.makedirs(workdir, exist_ok=True)
    only = set(x for x in (a.datasets or "").split(",") if x)
    if a.datasets_file:
        with open(a.datasets_file) as fh:
            only |= {x.strip() for x in fh if x.strip()}
    if a.agg.startswith("gs://") and not only:
        raise SystemExit("a gs:// --agg needs --datasets or --datasets-file: the prefix also holds "
                         "superseded and withdrawn versions")
    agg_dir = localize_agg(a.agg, only, f"{workdir}/agg")
    gencode_path = localize_file(a.gencode, f"{workdir}/ref")
    names, skipped = discover(agg_dir, only or None)
    # a requested version with no output at all is as incomplete as one with no marker
    skipped = sorted(set(skipped) | (only - set(names) - set(skipped)))
    if not names:
        raise SystemExit(f"no completed datasets under {a.agg}")
    gencode = load_gencode(gencode_path)
    ct, ti = load_rollups(a.rollup_dir)
    inp = Inputs(agg_dir, names, gencode, ct, ti, a)

    gidx_filter, missing = None, []
    if a.genes:
        gidx_filter, missing = resolve_genes(gencode, a.genes.split(","))
    remote = a.out.startswith("gs://")
    final = f"{a.out.rstrip('/')}/{a.release_date}"
    if a.genes:
        # a partial release must never land where the real one is read from
        # a sibling of the release, never inside it: the release manifest would not list it
        final += "-genes-mode"
        if remote and not a.force_remote_genes:
            raise SystemExit("--genes with a gs:// --out writes a partial release beside the real one; "
                             "use a local --out, or --force-remote-genes")
    # built in a fresh directory and moved into place, so nothing from an earlier run survives
    parent = f"{workdir}/release" if remote else os.path.dirname(final)
    os.makedirs(parent, exist_ok=True)
    local = tempfile.mkdtemp(prefix=f".{a.release_date}-", dir=parent)
    stage = None
    if gidx_filter is None:
        stage = f"{workdir}/stage"
        shutil.rmtree(stage, ignore_errors=True)
        split_by_chrom(inp, stage)

    rel = Release(local, a)
    profiles, summaries = [], []
    chroms = sorted(set(inp.map["chrom"].to_list()))
    for c in chroms:
        r = reduce_chrom(inp, c, a, stage, gidx_filter)
        if r is None:
            continue
        summary, pdr, prof = r
        rel.append("summary", summary)
        rel.append("per_dataset", pdr)
        if a.marker_panel:
            profiles.append(prof)
            summaries.append(summary)
        print(f"chrom {c}: summary {summary.height} per_dataset {pdr.height}", file=sys.stderr, flush=True)
    rel.close()
    if stage:
        shutil.rmtree(stage, ignore_errors=True)

    groups = Pipe(GZIP, f"{local}/bq/groups.tsv.gz", "\t".join(BQ_GROUP_COLUMNS) + "\n")
    inp.groups.select(BQ_GROUP_COLUMNS).sort("group_key").write_csv(
        groups.proc.stdin, separator="\t", include_header=False, null_value="")
    groups.close()
    listing = None
    if a.listing:
        with open(a.listing) as fh:
            listing = json.load(fh)
    measured = Pipe(GZIP, f"{local}/measured_genes.tsv.gz", "dataset_version_id\tgene_id\n")
    inp.map.join(inp.gencode.select("gidx", "gene_id"), on="gidx").select(
        pl.col("dataset_key").map_elements(lambda d: names[d], return_dtype=pl.Utf8).alias("dataset_version_id"),
        "gene_id").sort("dataset_version_id", "gene_id").write_csv(
        measured.proc.stdin, separator="\t", include_header=False)
    measured.close()
    datasets_table(inp, listing).write_csv(f"{local}/bq/datasets.tsv", separator="\t", null_value="")

    unmatched = sorted(set().union(*(set(s["unmatched"]) for s in inp.gene_stats)))
    manifest = {
        "release_date": a.release_date, "release_format_version": RELEASE_FORMAT_VERSION,
        "aggregation_format_version": aggregate.FORMAT_VERSION,
        "input_dataset_version_ids": names, "skipped_incomplete": skipped,
        "rollup_sha256": {n: sha256(f"{a.rollup_dir}/{n}") for n in
                          ("cell_type_rollup.tsv", "tissue_rollup.tsv", "metadata.json")
                          if os.path.exists(f"{a.rollup_dir}/{n}")},
        "rollup_unmapped": inp.unmapped,
        "gencode": {"source": a.gencode, "sha256": sha256(gencode_path),
                    "dropped_no_v49_match": len(unmatched), "dropped_ids": unmatched[:1000],
                    "dropped_non_ensembl_features": sum(s["non_ensembl"] for s in inp.gene_stats),
                    "per_dataset_dropped": {n: len(s["unmatched"]) for n, s in zip(names, inp.gene_stats)}},
        "thresholds": {"min_cells": a.min_cells, "min_donors": a.min_donors,
                       "min_pseudobulk_donors": a.min_pseudobulk_donors,
                       "summary_include_spatial": a.summary_include_spatial},
        "genes_mode": {"requested": a.genes.split(","), "not_in_gencode": missing} if a.genes else None,
        "rows": {**rel.rows, "groups": inp.groups.height, "datasets": len(names),
                 "summary_rows_defined": inp.rows.height},
        "files": {}, "index": GENE_POINT_INDEX.label,
    }
    panel = marker_panel(pl.concat(profiles), pl.concat(summaries), a, local) if a.marker_panel else None
    for root, _, files in os.walk(local):
        for f in files:
            p = os.path.join(root, f)
            manifest["files"][os.path.relpath(p, local)] = {"bytes": os.path.getsize(p), "sha256": sha256(p)}
    manifest["files"] = dict(sorted(manifest["files"].items()))
    with open(f"{local}/manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=1, sort_keys=True)
        fh.write("\n")
    if remote:
        upload(local, final, sorted(manifest["files"]) + ["manifest.json"])
        shutil.rmtree(local)
    else:
        os.chmod(local, 0o775)
        if os.path.isdir(final):
            shutil.rmtree(final)
        os.replace(local, final)
    return manifest, panel


# ---------------------------------------------------------------- marker panel

def marker_panel(profile, summary, a, local):
    """The proving ground's ~40-gene panel against the shipped roll-up: % of covered markers
    whose top cell type is the expected type or a descendant of it.

    Two levels: `celltype` = the tau profile (pooled over tissues, as the proving ground
    measured), `summary` = the top shipped summary row (cell type x tissue_general).
    """
    import proving_ground as pg
    parents, names = pg.load_cl(a.cl_obo)
    prof = profile.filter("pass_")
    res = []
    for marker, exp in pg.PANEL.items():
        rec = {"marker": marker, "expected": f"{exp} {names.get(exp, '')}"}
        for lvl, tab in (("celltype", prof), ("summary", summary)):
            t = tab.filter(pl.col("gene_name") == marker)
            rec[f"{lvl}_covered"] = any(exp in pg.ancestors(x, parents) for x in t["cell_type_rollup_ontology_term_id"])
            for metric in ("mean_log1p_cp10k", "frac_expressing", "donor_pseudobulk_mean_cp10k"):
                s = t.filter(pl.col(metric).is_not_null()).sort(metric, "cell_type_rollup_ontology_term_id",
                                                                 descending=[True, False])
                k = f"{lvl}_{metric}"
                if s.is_empty():
                    rec[k + "_top"], rec[k + "_ok"] = "absent", False
                    continue
                top = s.row(0, named=True)
                where = f" / {top['tissue_general']}" if lvl == "summary" else ""
                rec[k + "_top"] = f"{top['cell_type_rollup_ontology_term_id']} {top['cell_type_rollup']}{where} ({top[metric]:.3g})"
                rec[k + "_ok"] = exp in pg.ancestors(top["cell_type_rollup_ontology_term_id"], parents)
        res.append(rec)
    res = pl.DataFrame(res)
    out = {"n_markers": res.height}
    for lvl in ("celltype", "summary"):
        cov = res.filter(pl.col(f"{lvl}_covered"))
        out[f"{lvl}_n_covered"] = cov.height
        for metric in ("mean_log1p_cp10k", "frac_expressing", "donor_pseudobulk_mean_cp10k"):
            out[f"{lvl}_{metric}_pct"] = round(100 * cov[f"{lvl}_{metric}_ok"].mean(), 1) if cov.height else None
    res.write_csv(f"{local}/marker_panel.tsv", separator="\t")
    with open(f"{local}/marker_panel.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out))
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--agg", required=True, help="dir or gs:// prefix holding <dataset_version_id>/ worker outputs")
    p.add_argument("--out", required=True, help="releases root, local or gs://; output goes to <out>/<release-date>/")
    p.add_argument("--release-date", default=datetime.date.today().isoformat(),
                   help="release name: the directory under --out and manifest's release_date, used as given")
    p.add_argument("--datasets", help="comma-separated dataset_version_ids to publish (default for a local "
                   "--agg: every completed one; a gs:// --agg needs this or --datasets-file)")
    p.add_argument("--datasets-file", help="file with one dataset_version_id per line; adds to --datasets")
    p.add_argument("--listing", help="curation API datasets JSON, for datasets.tsv metadata")
    p.add_argument("--gencode", default=GENCODE_URI)
    p.add_argument("--rollup-dir", default=str(HERE / "rollup"))
    p.add_argument("--workdir")
    p.add_argument("--genes", help="comma-separated gene names or Ensembl ids: compute only these")
    p.add_argument("--force-remote-genes", action="store_true",
                   help="allow --genes with a gs:// --out (output goes to <out>/<date>-genes-mode/)")
    p.add_argument("--marker-panel", action="store_true", help="--genes = proving_ground.PANEL, then score it")
    p.add_argument("--cl-obo", help="CL release for the marker panel's descendant check")
    p.add_argument("--min-cells", type=int, default=MIN_CELLS)
    p.add_argument("--min-donors", type=int, default=MIN_DONORS)
    p.add_argument("--min-pseudobulk-donors", type=int, default=MIN_PB_DONORS)
    p.add_argument("--summary-include-spatial", action="store_true")
    a = p.parse_args(argv)
    if a.marker_panel:
        import proving_ground as pg
        if not a.cl_obo:
            p.error("--marker-panel needs --cl-obo")
        a.genes = a.genes or ",".join(pg.PANEL)
    return run(a)


if __name__ == "__main__":
    main()
