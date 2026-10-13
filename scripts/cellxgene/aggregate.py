#!/usr/bin/env python3
"""Per-dataset CELLxGENE worker: one schema-7 H5AD in, sparse per-group expression out.

    aggregate.py --input <local.h5ad | https URL> --out <local dir | gs://bucket/prefix>
                 --dataset-id ID --dataset-version-id ID

The chunked-CSR walk is the one `proving_ground.py` measured (single-threaded, 7.8-12.1 M
nnz/s, 5.4 GB peak on the heaviest dataset); that script stays as the record of the
measurement and this one is what runs. Pooling across datasets and ontology roll-up are not
done here: the output is one dataset's native groups, and the publish step reduces them.

Input
- read with h5py only, never `anndata.read_h5ad`, and only the obs columns used
- matrix: `raw/X` when present, else `X`. Either must hold integer counts: a float matrix
  with any non-integer, negative or NaN value is rejected (exit 3), which is what refuses a
  normalised `X` that has no `raw/X` beside it. Some datasets keep raw counts in `X` with no
  `raw/X` at all; the same check admits those
- CSR (or dense) only; rows are read in nnz-bounded chunks through `indptr`

Cells kept
- `is_primary_data` true, so a cell published in several datasets is counted once
- `assay_ontology_term_id` in `ASSAYS`. Cells of any other assay are dropped and counted per
  assay in `_SUCCESS.json`; an assay nobody has classified is dropped rather than guessed
  at, so a new one needs adding to the allow-list. Spatial assays are kept and flagged

Groups: the distinct values of `GROUP_KEYS` among kept cells, with the schema's label column
beside each ontology id. Within a group, per cell: CP10k = count * 1e4 / library size, then
log1p. Library size is the cell's total over every feature in the matrix.

Output files (Parquet, plus `_SUCCESS.json` written last as the completion marker)
- `groups.parquet`: one row per group — keys, labels, `spatial`, `n_cells`, `n_donors`
- `group_gene.parquet`: `n_expressing`, `sum_log1p_cp10k`, `mean_log1p_cp10k` (mean over all
  `n_cells`, zeros included) — **only rows with n_expressing >= 1**
- `group_donor.parquet`: per (group, donor) `n_cells` and `library_size` (summed counts)
- `group_donor_gene.parquet`: per (group, donor, gene) summed raw counts, rows > 0 only. The
  donor pseudobulk is computed from these by the publish step, which decides the minimum
  donor count; nothing here is filtered by donor count
- `genes.parquet`: every feature of the matrix's var, with `measured`

NA versus 0 — the contract consumers must keep: a gene that is `measured` here and has no
`group_gene` row for a group was measured and not seen, so its value is 0. A gene absent from
`genes.parquet`, or present with `measured` false, was not measured in this dataset, and its
value is NA. `expression_matrix()` implements exactly that and the tests hold it to it.
`measured` is false only for features `feature_is_filtered` zeroed out of `X` when `X` is
the matrix read.

Memory: three dense accumulators per gene — group x gene float64 sum and int32 count (12 bytes
per group) and group x donor float32 counts (4 bytes per pair). Together they are capped at
`--max-accum-bytes`; when they would exceed it, genes are split into blocks and the matrix is
re-read once per further block, trading I/O on the rare huge dataset for a bounded footprint.
The estimated peak is that cap plus the per-cell obs arrays and one chunk's transients. The
proving ground — not this script — measured the heaviest dataset at 5.4 GB peak with a 3 GB
accumulator, so roughly 2.4 GB of everything else; with the default 4 GB cap the estimate is
~6.5 GB, inside an 8 GB e2-standard-2. Not yet measured with this script on that dataset.

A dataset with no kept cells (all non-primary, or only excluded assays) is a valid empty
result: empty outputs and a `_SUCCESS.json` with `cells_kept` 0, so it is not re-fetched.

Exit codes: 0 done, 3 input rejected (retrying cannot help: non-count matrix, missing
mandatory column, unsupported encoding), anything else a failure that may be transient.
"""
import argparse
import concurrent.futures as cf
import json
import os
import resource
import shutil
import sys
import tempfile
import time

import numpy as np

# bump when any output file's columns or meaning change; the publish step refuses a mix
FORMAT_VERSION = 1
CHUNK_NNZ = 25_000_000
MAX_ACCUM_BYTES = 4_000_000_000

GROUP_KEYS = ["tissue_ontology_term_id", "tissue_type", "cell_type_ontology_term_id",
              "disease_ontology_term_id", "assay_ontology_term_id", "suspension_type"]
LABEL_OF = {k: k.removesuffix("_ontology_term_id") for k in GROUP_KEYS if k.endswith("_ontology_term_id")}

SPATIAL_ASSAYS = {"EFO:0022857": "Visium Spatial Gene Expression V1", "EFO:0030062": "Slide-seqV2"}
# RNA expression assays seen in human CELLxGENE datasets with primary cells (curation API,
# 2026-10-13). Deliberately absent: methylation (snmC-Seq2 EFO:0030027, snm3C-seq
# EFO:0022492, mCT-seq EFO:0030060) and chromatin (10x scATAC-seq EFO:0030007), whose
# matrices are not transcript counts. Also absent: BD Rhapsody Targeted mRNA (EFO:0700004), a
# gene panel whose CP10k is not comparable with whole-transcriptome assays (user ruling 2026-10-13)
RNA_ASSAYS = {
    "EFO:0008679": "CEL-seq", "EFO:0008722": "Drop-seq", "EFO:0008780": "inDrop",
    "EFO:0008796": "MARS-seq", "EFO:0008853": "Patch-seq", "EFO:0008877": "Quartz-seq",
    "EFO:0008919": "Seq-Well", "EFO:0008930": "Smart-seq", "EFO:0008931": "Smart-seq2",
    "EFO:0008953": "STRT-seq", "EFO:0009899": "10x 3' v2", "EFO:0009900": "10x 5' v2",
    "EFO:0009901": "10x 3' v1", "EFO:0009919": "SPLiT-seq", "EFO:0009922": "10x 3' v3",
    "EFO:0010010": "CEL-seq2", "EFO:0010058": "Fluidigm C1-based SMARTer library preparation",
    "EFO:0011025": "10x 5' v1", "EFO:0022488": "Smart-seq3",
    "EFO:0022490": "ScaleBio single cell RNA sequencing",
    "EFO:0022601": "Parse Evercode Whole Transcriptome v2", "EFO:0022604": "10x 3' v4",
    "EFO:0022605": "10x 5' v3", "EFO:0022606": "10x gene expression flex",
    "EFO:0022845": "modified STRT-seq", "EFO:0030002": "microwell-seq",
    "EFO:0030003": "10x 3' transcription profiling", "EFO:0030004": "10x 5' transcription profiling",
    "EFO:0030019": "Seq-Well S3", "EFO:0030026": "sci-Plex", "EFO:0030028": "sci-RNA-seq3",
    "EFO:0030059": "10x multiome", "EFO:0030074": "SORT-seq", "EFO:0030080": "10x transcription profiling",
    "EFO:0700003": "BD Rhapsody Whole Transcriptome Analysis",
    "EFO:0700010": "TruDrop", "EFO:0700011": "GEXSCOPE technology", "EFO:0700016": "Smart-seq v4",
    "EFO:0900000": "particle-templated instant partition sequencing", "EFO:0900001": "Asteria scRNA-seq kit",
    "EFO:0900002": "HIVE CLX Single-Cell RNAseq Solution", "EFO:0920087": "10x Next GEM Flex v1",
    "EFO:0920088": "10x GEM-X Flex v1", "EFO:0920089": "10x Flex Apex",
}
ASSAYS = RNA_ASSAYS | SPATIAL_ASSAYS

OUTPUT_FILES = ["groups.parquet", "group_gene.parquet", "group_donor.parquet",
                "group_donor_gene.parquet", "genes.parquet"]
# rows of a dense accumulator turned into Parquet per write, so the int64 index arrays of
# np.nonzero never span the whole accumulator
WRITE_CELLS = 50_000_000


class InputRejected(Exception):
    pass


def rss_mb():
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)


class Phases:
    def __init__(self):
        self.t, self.out = time.time(), {}

    def __call__(self, name):
        now = time.time()
        self.out[name] = {"s": round(now - self.t, 2), "peak_rss_mb": rss_mb()}
        print(name, self.out[name], file=sys.stderr, flush=True)
        self.t = now


# ---------------------------------------------------------------- source

def head(url):
    import requests
    r = requests.head(url, allow_redirects=True, timeout=60)
    r.raise_for_status()
    return {"url": url, "bytes": int(r.headers["Content-Length"]), "etag": r.headers.get("ETag"),
            "last_modified": r.headers.get("Last-Modified")}


def download(url, dest, parts=16, piece=64 << 20):
    """Parallel byte-range GET, pinned to the HEAD's ETag so a mid-download replacement fails.

    `sumstat_utils.fetch` is single-stream; the multi-GB files here need the parallel ranges
    to reach the measured ~275 MB/s.
    """
    import requests
    src = head(url)
    size = src["bytes"]
    pin = {"If-Match": src["etag"]} if src["etag"] else {}
    fd = os.open(dest, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    os.ftruncate(fd, size)

    def get(start):
        end = min(start + piece, size) - 1
        for attempt in range(5):
            try:
                with requests.get(url, headers={"Range": f"bytes={start}-{end}", **pin},
                                  stream=True, timeout=120) as r:
                    r.raise_for_status()
                    off = start
                    for buf in r.iter_content(4 << 20):
                        os.pwrite(fd, buf, off)
                        off += len(buf)
                    if off != end + 1:
                        raise IOError(f"short read {off} != {end + 1}")
                    return
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(2 ** attempt)

    try:
        with cf.ThreadPoolExecutor(parts) as ex:
            list(ex.map(get, range(0, size, piece)))
    finally:
        os.close(fd)
    return src


# ---------------------------------------------------------------- h5ad reading

def _decode(a):
    """h5py hands strings back as bytes; anndata writes them as UTF-8 ("Sjögren syndrome")."""
    a = np.asarray(a)
    if a.dtype.kind == "O":
        if not (a.size and isinstance(a.flat[0], bytes)):
            return a.astype(str)
        a = a.astype("S")
    return np.char.decode(a, "utf-8") if a.dtype.kind == "S" else a


def read_col(df, name):
    """(codes, categories) for any column encoding anndata has used."""
    import pandas as pd
    if name not in df:
        raise InputRejected(f"mandatory column {df.name}/{name} is missing")
    node = df[name]
    if hasattr(node, "keys"):
        codes, cats = node["codes"][()].astype(np.int64), _decode(node["categories"][()])
    elif "categories" in node.attrs:
        codes, cats = node[()].astype(np.int64), _decode(df.file[node.attrs["categories"]][()])
    else:
        codes, cats = pd.factorize(_decode(node[()]))
        codes = codes.astype(np.int64)
    if (codes < 0).any():
        # a -1 code would silently index the last category
        raise InputRejected(f"missing values in {df.name}/{name}")
    return codes, np.asarray(cats)


def col_values(df, name):
    codes, cats = read_col(df, name)
    return cats[codes]


def open_matrix(f):
    which = "raw/X" if "raw" in f and "X" in f["raw"] else "X"
    X = f[which]
    var = f["raw/var"] if which == "raw/X" else f["var"]
    if not hasattr(X, "keys"):
        return X, var, which, "dense", tuple(X.shape)
    enc = X.attrs.get("encoding-type", X.attrs.get("h5sparse_format", ""))
    enc = enc.decode() if isinstance(enc, bytes) else enc
    if enc not in ("csr_matrix", "csr"):
        raise InputRejected(f"unsupported sparse encoding {enc!r} in {which}")
    shape = tuple(int(s) for s in X.attrs.get("shape", X.attrs.get("h5sparse_shape")))
    return X, var, which, "csr", shape


def row_chunks(indptr, n_rows, layout, chunk_nnz):
    if layout == "dense":
        step = max(1, chunk_nnz // 20_000)
        for a in range(0, n_rows, step):
            yield a, min(a + step, n_rows)
        return
    a = 0
    while a < n_rows:
        b = int(np.searchsorted(indptr, indptr[a] + chunk_nnz, side="right")) - 1
        b = max(min(b, n_rows), a + 1)
        yield a, b
        a = b


def read_chunk(X, layout, indptr, a, b, n_genes):
    import scipy.sparse as sp
    if layout == "dense":
        return sp.csr_matrix(X[a:b])
    lo, hi = int(indptr[a]), int(indptr[b])
    return sp.csr_matrix((X["data"][lo:hi], X["indices"][lo:hi], indptr[a:b + 1] - lo),
                         shape=(b - a, n_genes))


def check_counts(data):
    bad = data < 0
    if data.dtype.kind == "f":
        bad = bad | (data != np.floor(data))  # NaN fails this too
    if bad.any():
        raise InputRejected("matrix holds non-integer or negative values, not raw counts "
                            "(a normalised X with no raw/X?)")


def indicator(rows, n_rows, keep_local, weights=None):
    """sparse (n_rows x chunk cells) matrix summing kept cells into their row."""
    import scipy.sparse as sp
    cols = np.flatnonzero(keep_local)
    w = np.ones(len(cols)) if weights is None else weights[cols]
    return sp.csr_matrix((w, (rows[cols], cols)), shape=(n_rows, len(keep_local)))


def add_dense(acc, m):
    m = m.tocoo()
    d = np.rint(m.data) if acc.dtype.kind == "i" else m.data
    acc[m.row, m.col] += d.astype(acc.dtype)


def gene_blocks(n_groups, n_gd, n_genes, max_bytes):
    per_gene = 12 * n_groups + 4 * n_gd  # float64 sum + int32 count per group, float32 per pair
    width = max(1, min(n_genes, int(max_bytes // max(per_gene, 1))))
    return [(lo, min(lo + width, n_genes)) for lo in range(0, n_genes, width)]


# ---------------------------------------------------------------- output

class SparseWriter:
    """Appends the nonzero cells of dense row blocks to one Parquet file."""

    def __init__(self, path, schema):
        import pyarrow.parquet as pq
        self.w = pq.ParquetWriter(path, schema)
        self.schema = schema
        self.rows = 0

    def write(self, cols):
        import pyarrow as pa
        self.w.write_table(pa.table(cols, schema=self.schema))
        self.rows += len(next(iter(cols.values())))

    def close(self):
        self.w.close()


def row_blocks(n_rows, n_cols):
    step = max(1, WRITE_CELLS // max(n_cols, 1))
    return [(a, min(a + step, n_rows)) for a in range(0, n_rows, step)]


# ---------------------------------------------------------------- aggregate

def aggregate(path, outdir, chunk_nnz=CHUNK_NNZ, max_accum_bytes=MAX_ACCUM_BYTES):
    """Write the per-dataset Parquet files into outdir; return the metrics for _SUCCESS.json."""
    import h5py
    import pandas as pd
    import pyarrow as pa

    ph = Phases()
    f = h5py.File(path, "r", rdcc_nbytes=256 << 20)
    X, var, which, layout, (n_cells, n_genes) = open_matrix(f)
    obs = f["obs"]
    n_obs = obs[obs.attrs.get("_index", "_index")].shape[0]
    if n_obs != n_cells:
        raise InputRejected(f"obs has {n_obs} rows, {which} has {n_cells}")

    primary = np.isin(col_values(obs, "is_primary_data").astype(str), ["True", "true", "1"])
    cols = {c: read_col(obs, c) for c in GROUP_KEYS + ["donor_id"]}
    assay = cols["assay_ontology_term_id"][1][cols["assay_ontology_term_id"][0]]
    assay_ok = np.isin(assay, list(ASSAYS))
    keep = primary & assay_ok
    excluded = pd.Series(assay[primary & ~assay_ok]).value_counts()
    kidx = np.flatnonzero(keep)
    grp_u, inv = np.unique(np.stack([cols[c][0][kidx] for c in GROUP_KEYS], 1).reshape(-1, len(GROUP_KEYS)),
                           axis=0, return_inverse=True)
    grp = np.full(n_cells, -1, np.int64)
    grp[kidx] = inv.ravel()
    gd_u, inv = np.unique(np.stack([grp[kidx], cols["donor_id"][0][kidx]], 1).reshape(-1, 2),
                          axis=0, return_inverse=True)
    gd = np.full(n_cells, -1, np.int64)
    gd[kidx] = inv.ravel()
    n_groups, n_gd = len(grp_u), len(gd_u)

    gene_ids = _decode(var[var.attrs.get("_index", "_index")][()])
    if len(gene_ids) != n_genes:
        raise InputRejected(f"{which} has {n_genes} columns, its var {len(gene_ids)} rows")
    genes = pd.DataFrame({"gene": np.arange(n_genes, dtype=np.int32), "gene_id": gene_ids})
    for c in ("feature_name", "feature_biotype"):
        # raw/var may lack the annotation columns; var carries them when the index matches
        src = var if c in var else f["var"]
        if c in src and np.array_equal(_decode(src[src.attrs.get("_index", "_index")][()]), gene_ids):
            genes[c] = col_values(src, c)
    measured = np.ones(n_genes, bool)
    if which == "X" and "feature_is_filtered" in var:
        measured = np.asarray(col_values(var, "feature_is_filtered")).astype(str) != "True"
    genes["measured"] = measured

    indptr = X["indptr"][()].astype(np.int64) if layout == "csr" else None
    # no kept cells: nothing to accumulate, so no pass over the matrix at all
    blocks = gene_blocks(n_groups, n_gd, n_genes, max_accum_bytes) if len(kidx) else []
    met = {"matrix": which, "layout": layout,
           "dtype": str((X["data"] if layout == "csr" else X).dtype),
           "n_cells": int(n_cells), "cells_kept": int(len(kidx)),
           "cells_not_primary": int((~primary).sum()),
           "cells_assay_excluded": {k: int(v) for k, v in excluded.items()},
           "n_groups": int(n_groups), "n_group_donor": int(n_gd), "n_genes": int(n_genes),
           "n_genes_measured": int(measured.sum()),
           "nnz": int(indptr[-1]) if layout == "csr" else None,
           "gene_blocks": len(blocks), "chunk_nnz": chunk_nnz, "max_accum_bytes": max_accum_bytes}
    ph("read_obs_var")

    lib = np.zeros(n_cells)
    n_cells_g = np.bincount(grp[kidx], minlength=n_groups)
    gdg = SparseWriter(f"{outdir}/group_donor_gene.parquet", pa.schema(
        [("gd", pa.int32()), ("gene", pa.int32()), ("counts", pa.float32())]))
    gg = SparseWriter(f"{outdir}/group_gene.parquet", pa.schema(
        [("group", pa.int32()), ("gene", pa.int32()), ("n_expressing", pa.int32()),
         ("sum_log1p_cp10k", pa.float64()), ("mean_log1p_cp10k", pa.float32())]))
    for bi, (glo, ghi) in enumerate(blocks):
        # library sizes need whole rows, so they and the count check come from the first pass
        first = bi == 0
        width = ghi - glo
        sum_log = np.zeros((n_groups, width))
        n_expr = np.zeros((n_groups, width), np.int32)
        gd_acc = np.zeros((n_gd, width), np.float32)
        for a, b in row_chunks(indptr, n_cells, layout, chunk_nnz):
            m = read_chunk(X, layout, indptr, a, b, n_genes)
            if first:
                check_counts(m.data)
                lib[a:b] = np.asarray(m.sum(1)).ravel()
            k = keep[a:b]
            if not k.any():
                continue
            if width != n_genes:
                m = m[:, glo:ghi]
            A = indicator(grp[a:b], n_groups, k)
            with np.errstate(divide="ignore"):
                scale = np.where(lib[a:b] > 0, 1e4 / lib[a:b], 0.0)
            logm = m.astype(np.float64)
            logm.data = np.log1p(logm.data * np.repeat(scale, np.diff(m.indptr)))
            add_dense(sum_log, A @ logm)
            del logm
            ex = m.copy()
            ex.data = (ex.data > 0).astype(np.float64)
            add_dense(n_expr, A @ ex)
            del ex
            add_dense(gd_acc, indicator(gd[a:b], n_gd, k) @ m)
        for a, b in row_blocks(n_gd, width):
            r, c = np.nonzero(gd_acc[a:b])
            gdg.write({"gd": (r + a).astype(np.int32), "gene": (c + glo).astype(np.int32),
                       "counts": gd_acc[a:b][r, c]})
        del gd_acc
        for a, b in row_blocks(n_groups, width):
            r, c = np.nonzero(n_expr[a:b])
            sl = sum_log[a:b][r, c]
            gg.write({"group": (r + a).astype(np.int32), "gene": (c + glo).astype(np.int32),
                      "n_expressing": n_expr[a:b][r, c], "sum_log1p_cp10k": sl,
                      "mean_log1p_cp10k": (sl / n_cells_g[r + a]).astype(np.float32)})
        del sum_log, n_expr
        ph("main_pass" if first else f"pass_{bi}")
    gdg.close()
    gg.close()

    gmeta = pd.DataFrame({"group": np.arange(n_groups, dtype=np.int32)})
    for i, c in enumerate(GROUP_KEYS):
        gmeta[c] = cols[c][1][grp_u[:, i]]
        if c in LABEL_OF:
            ids = cols[c][1][cols[c][0][kidx]]
            labels = col_values(obs, LABEL_OF[c])[kidx]
            gmeta[LABEL_OF[c]] = gmeta[c].map(pd.Series(labels).groupby(ids).first())
    gmeta["spatial"] = gmeta.assay_ontology_term_id.isin(list(SPATIAL_ASSAYS))
    gmeta["n_cells"] = n_cells_g.astype(np.int32)
    gmeta["n_donors"] = np.bincount(gd_u[:, 0], minlength=n_groups).astype(np.int32)
    gmeta.to_parquet(f"{outdir}/groups.parquet", index=False)

    pd.DataFrame({"gd": np.arange(n_gd, dtype=np.int32), "group": gd_u[:, 0].astype(np.int32),
                  "donor_id": cols["donor_id"][1][gd_u[:, 1]],
                  "n_cells": np.bincount(gd[kidx], minlength=n_gd).astype(np.int32),
                  "library_size": np.bincount(gd[kidx], weights=lib[kidx], minlength=n_gd)}
                 ).to_parquet(f"{outdir}/group_donor.parquet", index=False)
    genes.to_parquet(f"{outdir}/genes.parquet", index=False)
    met["group_gene_rows"], met["group_donor_gene_rows"] = gg.rows, gdg.rows
    ph("write")
    f.close()
    met["phases"] = ph.out
    if met["nnz"] and blocks:
        passes = sum(p["s"] for k, p in ph.out.items() if k == "main_pass" or k.startswith("pass_"))
        met["nnz_per_s"] = round(met["nnz"] * len(blocks) / max(passes, 1e-9))
    return met


def expression_matrix(outdir, gene_ids):
    """(groups, n_groups x len(gene_ids) mean log1p CP10k) with NaN for genes not measured.

    The reference reading of the NA-versus-0 contract stated in the module docstring.
    """
    import pandas as pd
    groups = pd.read_parquet(f"{outdir}/groups.parquet")
    genes = pd.read_parquet(f"{outdir}/genes.parquet")
    gg = pd.read_parquet(f"{outdir}/group_gene.parquet")
    pos = {g: i for i, g in enumerate(gene_ids)}
    out = np.full((len(groups), len(gene_ids)), np.nan)
    measured = genes[genes.measured & genes.gene_id.isin(pos)]
    out[:, [pos[g] for g in measured.gene_id]] = 0.0
    col = pd.Series(pos).reindex(genes.gene_id.values).values  # var position -> output column
    hit = gg[~np.isnan(col[gg.gene.values]) & genes.measured.values[gg.gene.values]]
    out[hit.group.values, col[hit.gene.values].astype(int)] = hit.mean_log1p_cp10k.values
    return groups, out


# ---------------------------------------------------------------- driver

def _gcs_bucket(uri):
    from google.cloud import storage
    bucket, _, prefix = uri[len("gs://"):].partition("/")
    return storage.Client().bucket(bucket), prefix.rstrip("/")


def publish(local, out, names):
    """Copy names from local to out; `_SUCCESS.json` must be last so it marks completeness."""
    if not out.startswith("gs://"):
        return
    bucket, prefix = _gcs_bucket(out)
    for n in names:
        bucket.blob(f"{prefix}/{n}").upload_from_filename(os.path.join(local, n))


def clear_marker(out):
    """A rerun must not leave the previous run's marker beside a half-written new output."""
    if out.startswith("gs://"):
        bucket, prefix = _gcs_bucket(out)
        blob = bucket.blob(f"{prefix}/_SUCCESS.json")
        if blob.exists():
            blob.delete()
    elif os.path.exists(f"{out}/_SUCCESS.json"):
        os.remove(f"{out}/_SUCCESS.json")


def run(args):
    t0 = time.time()
    remote_out = args.out.startswith("gs://")
    if not remote_out:
        os.makedirs(args.out, exist_ok=True)
    clear_marker(args.out)
    local = tempfile.mkdtemp(prefix="cxg-out-", dir=args.workdir) if remote_out else args.out
    fetched = None
    try:
        dl = {}
        if args.input.startswith(("http://", "https://")):
            fetched = os.path.join(tempfile.mkdtemp(prefix="cxg-in-", dir=args.workdir), "input.h5ad")
            t = time.time()
            source = download(args.input, fetched)
            dl = {"download": {"s": round(time.time() - t, 2), "peak_rss_mb": rss_mb()}}
            path = fetched
        else:
            path = args.input
            source = {"url": None, "bytes": os.path.getsize(path), "etag": None, "last_modified": None}
        met = aggregate(path, local, chunk_nnz=args.chunk_nnz, max_accum_bytes=args.max_accum_bytes)
        met["phases"] = dl | met["phases"]
        if fetched:
            shutil.rmtree(os.path.dirname(fetched))
            fetched = None
        t = time.time()
        publish(local, args.out, OUTPUT_FILES)
        if remote_out:
            met["phases"]["upload"] = {"s": round(time.time() - t, 2), "peak_rss_mb": rss_mb()}
        success = {"aggregation_format_version": FORMAT_VERSION,
                   "dataset_id": args.dataset_id, "dataset_version_id": args.dataset_version_id,
                   "source": source, "worker_image": os.environ.get("WORKER_IMAGE_DIGEST"),
                   "outputs": {n: os.path.getsize(os.path.join(local, n)) for n in OUTPUT_FILES},
                   **met, "total_s": round(time.time() - t0, 2), "peak_rss_mb": rss_mb()}
        tmp = os.path.join(local, "_SUCCESS.json.tmp")
        with open(tmp, "w") as fh:
            json.dump(success, fh, indent=1)
        os.replace(tmp, os.path.join(local, "_SUCCESS.json"))
        publish(local, args.out, ["_SUCCESS.json"])
        print(json.dumps(success))
        return success
    finally:
        if fetched:
            shutil.rmtree(os.path.dirname(fetched), ignore_errors=True)
        if remote_out:
            shutil.rmtree(local, ignore_errors=True)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--input", required=True, help="local .h5ad or https URL")
    p.add_argument("--out", required=True, help="local directory or gs://bucket/prefix")
    p.add_argument("--dataset-id", required=True)
    p.add_argument("--dataset-version-id", required=True)
    p.add_argument("--workdir", default=None, help="scratch for the download and staged output")
    p.add_argument("--chunk-nnz", type=int, default=CHUNK_NNZ)
    p.add_argument("--max-accum-bytes", type=int, default=MAX_ACCUM_BYTES)
    args = p.parse_args(argv)
    try:
        run(args)
    except InputRejected as e:
        print(f"input rejected: {e}", file=sys.stderr)
        sys.exit(3)


if __name__ == "__main__":
    main()
