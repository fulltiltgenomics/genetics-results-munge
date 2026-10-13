#!/usr/bin/env python3
"""Throwaway kill-criterion measurement for the CELLxGENE epic.

Not the production worker: it exists to measure throughput, memory and marker
validity on a handful of real H5ADs. Reads with h5py only (never
anndata.read_h5ad), walks CSR rows in nnz-bounded chunks via indptr, and never
densifies a sparse chunk.

  download URL DEST                 parallel byte-range GET
  process H5AD OUTDIR DATASET_ID    aggregate one dataset
  validate CL_OBO OUTDIR...         marker panel on the pooled normal-tissue summary
  gencode GENCODE_TSV OUTDIR...     CELLxGENE Ensembl ids with no v49 match
"""
import concurrent.futures as cf
import json
import os
import resource
import sys
import time

import numpy as np

CHUNK_NNZ = 25_000_000
# dense per-(group, donor) count accumulator only when it is this small;
# otherwise a first pass over X/data gives library sizes and the donor
# pseudobulk is accumulated already weighted, at group level
DENSE_GD_BYTES = 3e9

OBS_KEY = ["tissue_ontology_term_id", "tissue_type", "cell_type_ontology_term_id",
           "disease_ontology_term_id", "assay_ontology_term_id", "suspension_type"]

# marker -> expected Cell Ontology term; the expected term or any descendant counts
PANEL = {
    "CD3E": "CL:0000084", "CD3D": "CL:0000084", "CD8A": "CL:0000084",
    "MS4A1": "CL:0000236", "CD79A": "CL:0000236", "CD19": "CL:0000236",
    "JCHAIN": "CL:0000786",
    "CD14": "CL:0000576", "LYZ": "CL:0000576",
    "NKG7": "CL:0000623", "GNLY": "CL:0000623",
    "C1QA": "CL:0000235",
    "CLEC4C": "CL:0000784",
    "TPSAB1": "CL:0000097",
    "HBB": "CL:0000232", "PPBP": "CL:0000233",
    "ALB": "CL:0000182", "APOA1": "CL:0000182",
    "SFTPC": "CL:0002063", "AGER": "CL:0002062",
    "INS": "CL:0000169", "GCG": "CL:0000171",
    "PRSS1": "CL:0002064", "CPA1": "CL:0002064",
    "MBP": "CL:0000128", "PLP1": "CL:0000128", "MOG": "CL:0000128",
    "GFAP": "CL:0000127", "AQP4": "CL:0000127",
    "SNAP25": "CL:0000540", "SYT1": "CL:0000540",
    "PECAM1": "CL:0000115", "VWF": "CL:0000115", "CDH5": "CL:0000115",
    "COL1A1": "CL:0000057", "DCN": "CL:0000057", "LUM": "CL:0000057",
    "KRT5": "CL:0000646",
    "MYH11": "CL:0000192",
    "MUC2": "CL:0000160",
}
# roll-up targets beyond the panel's own terms, so unrelated cells pool into
# broad competitors rather than staying as many small leaves
EXTRA_ROLLUP = {"CL:0000451": "dendritic cell", "CL:0000066": "epithelial cell",
                "CL:0000094": "granulocyte", "CL:0000669": "pericyte",
                "CL:0000125": "glial cell"}
MIN_CELLS = 100
NORMAL = "PATO:0000461"


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


# ---------------------------------------------------------------- download

def download(url, dest, parts=16, piece=64 << 20):
    import requests
    t0 = time.time()
    size = int(requests.head(url, allow_redirects=True, timeout=60).headers["Content-Length"])
    fd = os.open(dest, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    os.ftruncate(fd, size)

    def get(start):
        end = min(start + piece, size) - 1
        for attempt in range(5):
            try:
                with requests.get(url, headers={"Range": f"bytes={start}-{end}"}, stream=True, timeout=120) as r:
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

    with cf.ThreadPoolExecutor(parts) as ex:
        list(ex.map(get, range(0, size, piece)))
    os.close(fd)
    dt = time.time() - t0
    out = {"bytes": size, "download_s": round(dt, 1), "download_MBps": round(size / 1e6 / dt, 1)}
    print(json.dumps(out))
    return out


# ---------------------------------------------------------------- h5ad reading

def _decode(a):
    a = np.asarray(a)
    return a.astype(str) if a.dtype.kind in "SO" else a


def read_obs_col(obs, name):
    """Return (codes, categories) for any obs column encoding anndata has used."""
    import pandas as pd
    node = obs[name]
    if hasattr(node, "keys"):  # anndata >= 0.8 categorical group
        return node["codes"][()].astype(np.int64), _decode(node["categories"][()])
    if "categories" in node.attrs:  # anndata < 0.8: codes plus a reference
        return node[()].astype(np.int64), _decode(obs.file[node.attrs["categories"]][()])
    codes, cats = pd.factorize(_decode(node[()]))
    return codes.astype(np.int64), np.asarray(cats)


def var_index(var):
    return _decode(var[var.attrs.get("_index", "_index")][()])


def open_matrix(f):
    """raw/X if present else X; refuse anything but CSR or dense."""
    base = "raw" if "raw" in f and "X" in f["raw"] else None
    X = f["raw/X"] if base else f["X"]
    var = f["raw/var"] if base else f["var"]
    if hasattr(X, "keys"):
        enc = X.attrs.get("encoding-type", X.attrs.get("h5sparse_format", ""))
        enc = enc.decode() if isinstance(enc, bytes) else enc
        if enc not in ("csr_matrix", "csr"):
            raise SystemExit(f"unsupported sparse encoding {enc!r}")
        return X, var, ("raw/X" if base else "X"), "csr"
    return X, var, ("raw/X" if base else "X"), "dense"


def row_chunks(indptr, n_rows, layout):
    if layout == "dense":
        step = 5000
        for a in range(0, n_rows, step):
            yield a, min(a + step, n_rows)
        return
    a = 0
    while a < n_rows:
        b = int(np.searchsorted(indptr, indptr[a] + CHUNK_NNZ, side="right")) - 1
        b = max(min(b, n_rows), a + 1)
        yield a, b
        a = b


def read_chunk(X, layout, indptr, a, b, n_genes, values=True):
    import scipy.sparse as sp
    if layout == "dense":
        return sp.csr_matrix(X[a:b])  # stored dense on disk, so a row block is already dense
    lo, hi = int(indptr[a]), int(indptr[b])
    data = X["data"][lo:hi]
    if not values:
        return data
    idx = X["indices"][lo:hi]
    return sp.csr_matrix((data, idx, (indptr[a:b + 1] - lo)), shape=(b - a, n_genes))


def check_integer(data):
    if data.dtype.kind == "f" and not np.array_equal(data, np.floor(data)):
        raise SystemExit("non-integer counts in matrix: refusing (normalised X with no raw/X)")


def indicator(rows_group, weights, n_groups, keep_local):
    """sparse (n_groups x n_chunk_cells) matrix that sums cells into groups."""
    import scipy.sparse as sp
    cols = np.flatnonzero(keep_local)
    return sp.csr_matrix((weights[cols], (rows_group[cols], cols)), shape=(n_groups, len(keep_local)))


def add_dense(acc, m):
    m = m.tocoo()
    acc[m.row, m.col] += m.data


# ---------------------------------------------------------------- process

def process(path, outdir, dataset_id):
    import h5py
    import pandas as pd
    os.makedirs(outdir, exist_ok=True)
    met = {"dataset_id": dataset_id, "file_bytes": os.path.getsize(path), "phases": {}}
    t = time.time()

    def phase(name):
        nonlocal t
        now = time.time()
        met["phases"][name] = {"s": round(now - t, 2), "peak_rss_mb": round(rss_mb())}
        print(name, met["phases"][name], flush=True)
        t = now

    f = h5py.File(path, "r", rdcc_nbytes=256 << 20)
    X, var, which, layout = open_matrix(f)
    met["matrix"], met["layout"] = which, layout
    obs = f["obs"]
    n_cells = obs[obs.attrs.get("_index", "_index")].shape[0]

    prim_c, prim_cat = read_obs_col(obs, "is_primary_data")
    keep = np.isin(np.asarray(prim_cat).astype(str), ["True", "true", "1"])[prim_c]
    cols = {}
    for c in OBS_KEY + ["donor_id"]:
        if c in obs:
            cols[c] = read_obs_col(obs, c)
        else:
            cols[c] = (np.zeros(n_cells, np.int64), np.array(["na"]))
    labels = {c: read_obs_col(obs, c.replace("_ontology_term_id", "")) for c in
              ("cell_type_ontology_term_id", "tissue_ontology_term_id") if c.replace("_ontology_term_id", "") in obs}
    gkey = np.stack([cols[c][0] for c in OBS_KEY], 1)
    grp_u, grp = np.unique(gkey, axis=0, return_inverse=True)
    grp = grp.ravel()
    n_groups = len(grp_u)
    gd_u, gd = np.unique(np.stack([grp, cols["donor_id"][0]], 1), axis=0, return_inverse=True)
    gd = gd.ravel()
    n_gd = len(gd_u)
    n_donors_g = np.bincount(gd_u[:, 0], minlength=n_groups)
    genes = var_index(var)
    n_genes = len(genes)
    names = _decode(var["feature_name"][()]) if "feature_name" in var and not hasattr(var["feature_name"], "keys") else None
    if names is None and "feature_name" in var:
        c, cats = read_obs_col(var, "feature_name")
        names = cats[c]
    if names is None:  # raw/var may lack feature_name; X's var carries it in the same order
        c, cats = read_obs_col(f["var"], "feature_name")
        names = cats[c] if len(c) == n_genes else np.array([""] * n_genes)
    marker_cols = np.flatnonzero(np.isin(names, list(PANEL)))
    indptr = X["indptr"][()].astype(np.int64) if layout == "csr" else None
    nnz = int(indptr[-1]) if layout == "csr" else None
    met.update(n_cells=int(n_cells), cells_kept=int(keep.sum()), n_groups=int(n_groups), n_group_donor=int(n_gd),
               n_genes=int(n_genes), nnz=nnz, dtype=str((X["data"] if layout == "csr" else X).dtype))
    # kept-cell library sizes are needed before the pseudobulk weights exist
    dense_gd = n_gd * n_genes * 4 <= DENSE_GD_BYTES
    met["mode"] = "single-pass dense group-donor" if dense_gd else "two-pass weighted"
    phase("obs_var")
    print(json.dumps({k: met[k] for k in ("n_cells", "cells_kept", "n_groups", "n_group_donor", "n_genes", "nnz", "mode")}), flush=True)

    lib = np.zeros(n_cells)
    if not dense_gd:
        for a, b in row_chunks(indptr, n_cells, layout):
            data = read_chunk(X, layout, indptr, a, b, n_genes, values=False)
            if layout == "dense":
                data = data.tocsr()
                check_integer(data.data)
                lib[a:b] = np.asarray(data.sum(1)).ravel()
                continue
            check_integer(data)
            cs = np.concatenate([[0.0], np.cumsum(data, dtype=np.float64)])
            p = indptr[a:b + 1] - indptr[a]
            lib[a:b] = cs[p[1:]] - cs[p[:-1]]
        phase("pass1_libsize")

    sum_log = np.zeros((n_groups, n_genes))
    n_expr = np.zeros((n_groups, n_genes), np.float32)
    pb = np.zeros((n_groups, n_genes)) if not dense_gd else None
    gd_counts = np.zeros((n_gd, n_genes), np.float32) if dense_gd else None
    gd_lib = np.zeros(n_gd)
    gd_marker = np.zeros((n_gd, len(marker_cols)))
    if not dense_gd:
        gd_lib = np.bincount(gd[keep], weights=lib[keep], minlength=n_gd)
        w_cell = np.zeros(n_cells)
        nz = keep & (gd_lib[gd] > 0)
        w_cell[nz] = 1e4 / (gd_lib[gd[nz]] * n_donors_g[grp[nz]])
    for a, b in row_chunks(indptr, n_cells, layout):
        m = read_chunk(X, layout, indptr, a, b, n_genes)
        if dense_gd:
            check_integer(m.data)
            lib[a:b] = np.asarray(m.sum(1)).ravel()
        k = keep[a:b]
        if not k.any():
            continue
        ones = np.ones(b - a)
        g_loc, gd_loc = grp[a:b], gd[a:b]
        A = indicator(g_loc, ones, n_groups, k)
        with np.errstate(divide="ignore", invalid="ignore"):
            scale = np.where(lib[a:b] > 0, 1e4 / lib[a:b], 0.0)
        logm = m.copy()
        logm.data = np.log1p(logm.data * np.repeat(scale, np.diff(m.indptr)))
        add_dense(sum_log, A @ logm)
        ex = m.copy()
        ex.data = (ex.data > 0).astype(np.float32)
        add_dense(n_expr, A @ ex)
        Agd = indicator(gd_loc, ones, n_gd, k)
        gd_marker += (Agd @ m[:, marker_cols]).toarray()
        if dense_gd:
            add_dense(gd_counts, Agd @ m)
            gd_lib += np.bincount(gd_loc[k], weights=lib[a:b][k], minlength=n_gd)
        else:
            add_dense(pb, indicator(g_loc, w_cell[a:b], n_groups, k) @ m)
    phase("main_pass")

    n_cells_g = np.bincount(grp[keep], minlength=n_groups)
    if dense_gd:
        import scipy.sparse as sp
        with np.errstate(divide="ignore"):
            w = np.where(gd_lib > 0, 1e4 / gd_lib, 0.0) / np.maximum(n_donors_g[gd_u[:, 0]], 1)
        W = sp.csc_matrix((w, (gd_u[:, 0], np.arange(n_gd))), shape=(n_groups, n_gd))
        pb = np.zeros((n_groups, n_genes))
        # blocked so the float32 accumulator is never upcast to float64 in one piece
        for s in range(0, n_gd, 2000):
            pb += W[:, s:s + 2000] @ gd_counts[s:s + 2000].astype(np.float64)
        del gd_counts
    # n_donors counts only donors with kept cells in the group
    kept_gd = np.bincount(gd[keep], minlength=n_gd) > 0
    n_donors_kept = np.bincount(gd_u[kept_gd, 0], minlength=n_groups)
    gmeta = pd.DataFrame({c: cols[c][1][grp_u[:, i]] for i, c in enumerate(OBS_KEY)})
    for c, (codes, cats) in labels.items():
        lab = pd.Series(cats[codes]).groupby(cols[c][1][cols[c][0]]).first()
        gmeta[c.replace("_ontology_term_id", "")] = gmeta[c].map(lab)
    gmeta["n_cells"], gmeta["n_donors"] = n_cells_g, n_donors_kept
    gmeta.insert(0, "group", np.arange(n_groups))
    gmeta = gmeta[gmeta.n_cells > 0]
    gmeta.to_parquet(f"{outdir}/groups.parquet")
    r, c = np.nonzero(n_expr)
    pd.DataFrame({"group": r.astype(np.int32), "gene": c.astype(np.int32),
                  "mean_log_cp10k": (sum_log[r, c] / n_cells_g[r]).astype(np.float32),
                  "frac_expressing": (n_expr[r, c] / n_cells_g[r]).astype(np.float32),
                  "pseudobulk_cp10k": pb[r, c].astype(np.float32)}).to_parquet(f"{outdir}/group_gene.parquet")
    pd.DataFrame({"gene_id": genes, "feature_name": names}).to_parquet(f"{outdir}/genes.parquet")
    donors = cols["donor_id"][1][gd_u[:, 1]]
    mk = pd.DataFrame(gd_marker, columns=names[marker_cols])
    mk.insert(0, "group", gd_u[:, 0]); mk.insert(1, "donor", donors); mk.insert(2, "lib", gd_lib)
    mk[mk.lib > 0].to_parquet(f"{outdir}/markers_gd.parquet")
    met["group_gene_rows"] = int(len(r))
    phase("write")
    tot = sum(p["s"] for p in met["phases"].values())
    met["total_s"] = round(tot, 1)
    met["peak_rss_mb"] = round(rss_mb())
    if nnz:
        met["nnz_per_s_main"] = round(nnz / met["phases"]["main_pass"]["s"])
        met["nnz_per_s_total"] = round(nnz / tot)
    json.dump(met, open(f"{outdir}/metrics.json", "w"), indent=1)
    print(json.dumps(met))


# ---------------------------------------------------------------- validity

def load_cl(path):
    parents, names, cur = {}, {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur = None
        elif line.startswith("id: "):
            cur = line[4:]
            parents.setdefault(cur, set())
        elif cur and line.startswith("name: "):
            names[cur] = line[6:]
        elif cur and line.startswith("is_a: "):
            parents[cur].add(line[6:].split()[0])
    return parents, names


def ancestors(term, parents, memo={}):
    if term not in memo:
        out = {term}
        for p in parents.get(term, ()):
            out |= ancestors(p, parents, memo)
        memo[term] = out
    return memo[term]


def validate(obo, outdirs):
    import pandas as pd
    parents, names = load_cl(obo)
    targets = set(PANEL.values()) | set(EXTRA_ROLLUP)

    # post-hoc sensitivity only: PG_ROLLUP_OVERRIDE="CL:a=CL:b,..." forces a term's roll-up target
    override = dict(kv.split("=") for kv in os.environ.get("PG_ROLLUP_OVERRIDE", "").split(",") if kv)

    def rollup(term):
        anc = ancestors(term, parents) & targets
        for o, tgt in override.items():
            if o in ancestors(term, parents):
                return tgt
        if not anc:
            return term
        # most specific: a target that is not an ancestor of another candidate
        best = [t for t in anc if not any(t != o and t in ancestors(o, parents) for o in anc)]
        return sorted(best)[0]

    gm, mk_all, measured = [], [], []
    for d in outdirs:
        did = os.path.basename(d.rstrip("/"))
        g = pd.read_parquet(f"{d}/groups.parquet")
        gg = pd.read_parquet(f"{d}/group_gene.parquet")
        genes = pd.read_parquet(f"{d}/genes.parquet")
        sel = genes.index[genes.feature_name.isin(list(PANEL))]
        gg = gg[gg.gene.isin(sel)].copy()
        gg["marker"] = genes.feature_name.values[gg.gene]
        g = g[(g.disease_ontology_term_id == NORMAL) & (g.tissue_type.isin(["tissue", "na"]))
              & (g.cell_type_ontology_term_id != "unknown")].copy()
        g["dataset"] = did
        g["rolled"] = g.cell_type_ontology_term_id.map(rollup)
        gg = gg.merge(g[["group", "n_cells", "n_donors"]], on="group")
        gg["dataset"] = did
        gm.append((g, gg))
        # a gene missing from a dataset's var is unmeasured there, not zero
        for m in set(genes.feature_name) & set(PANEL):
            measured.append(g.groupby("rolled").n_cells.sum().rename("n").reset_index().assign(marker=m))
        mk = pd.read_parquet(f"{d}/markers_gd.parquet").merge(g[["group", "rolled"]], on="group")
        mk["dataset"] = did
        mk_all.append(mk)
    groups = pd.concat([x[0] for x in gm])
    gg = pd.concat([x[1] for x in gm]).merge(groups[["dataset", "group", "rolled"]], on=["dataset", "group"])
    cells = groups.groupby("rolled").n_cells.sum()
    eligible = cells[cells >= MIN_CELLS].index
    denom = pd.concat(measured).groupby(["marker", "rolled"]).n.sum().unstack()
    denom = denom.where(denom >= MIN_CELLS).reindex(columns=eligible)
    # mean log CP10k pooled over the normal groups of the rolled type in datasets that measured the gene
    gg["w"] = gg.mean_log_cp10k * gg.n_cells
    num = gg.groupby(["marker", "rolled"]).w.sum().unstack().reindex(index=denom.index, columns=eligible).fillna(0.0)
    mean_log = num / denom
    # donor pseudobulk: counts and library summed per (rolled type, dataset, donor), then CP10k averaged over donors
    mk = pd.concat(mk_all)
    mcols = [c for c in mk.columns if c in PANEL]
    dsum = mk.groupby(["rolled", "dataset", "donor"])[mcols + ["lib"]].sum(min_count=1)
    cp = dsum[mcols].div(dsum.lib, axis=0) * 1e4
    pbk = cp.groupby(level="rolled").mean().T.reindex(columns=eligible)
    pbk = pbk.where(denom.reindex(index=pbk.index).notna())
    n_don = dsum.groupby(level="rolled").size()

    rows = []
    for marker, exp in PANEL.items():
        present = set(denom.columns[denom.loc[marker].notna()]) if marker in denom.index else set()
        covered = any(exp in ancestors(t, parents) for t in present)
        rec = {"marker": marker, "expected": f"{exp} {names.get(exp, '')}", "covered": covered}
        for label, tab in (("meanlog", mean_log), ("pseudobulk", pbk)):
            if marker not in tab.index or tab.loc[marker].notna().sum() == 0:
                rec[f"{label}_top"], rec[f"{label}_ok"] = "absent", False
                continue
            s = tab.loc[marker].dropna().sort_values(ascending=False)
            top = s.index[0]
            rec[f"{label}_top"] = f"{top} {names.get(top, '')} ({s.iloc[0]:.2f})"
            rec[f"{label}_ok"] = exp in ancestors(top, parents)
            exp_rank = [i for i, t in enumerate(s.index, 1) if exp in ancestors(t, parents)]
            rec[f"{label}_exp_rank"] = exp_rank[0] if exp_rank else None
        rows.append(rec)
    res = pd.DataFrame(rows)
    cov = res[res.covered]
    summary = {"n_markers": len(res), "n_covered": int(len(cov)),
               "pct_meanlog_covered": round(100 * cov.meanlog_ok.mean(), 1),
               "pct_pseudobulk_covered": round(100 * cov.pseudobulk_ok.mean(), 1),
               "n_rolled_types_eligible": int(len(eligible)), "normal_cells": int(groups.n_cells.sum())}
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60); pd.set_option("display.max_rows", 200)
    print(res.to_string())
    print(json.dumps(summary))
    types = pd.DataFrame({"n_cells": cells.reindex(eligible), "n_donors": n_don.reindex(eligible)})
    types["name"] = [names.get(t, "") for t in types.index]
    print(types.sort_values("n_cells", ascending=False).to_string())
    res.to_csv(os.path.join(os.path.dirname(outdirs[0].rstrip("/")), "marker_panel.tsv"), sep="\t", index=False)
    types.to_csv(os.path.join(os.path.dirname(outdirs[0].rstrip("/")), "rolled_types.tsv"), sep="\t")


def gencode(tsv, outdirs):
    import pandas as pd
    v49 = set(pd.read_csv(tsv, sep="\t", usecols=["gene_id"]).gene_id.str.split(".").str[0])
    allids = set()
    for d in outdirs:
        ids = pd.read_parquet(f"{d}/genes.parquet").gene_id
        ens = ids[ids.str.startswith("ENSG")].str.split(".").str[0]
        miss = set(ens) - v49
        allids |= set(ens)
        print(os.path.basename(d.rstrip("/")), "var", len(ids), "ENSG", len(ens), "non-ENSG", int((~ids.str.startswith("ENSG")).sum()), "missing_v49", len(miss))
    miss = sorted(allids - v49)
    print("union ENSG", len(allids), "missing in v49", len(miss), miss[:20])


if __name__ == "__main__":
    cmd, *args = sys.argv[1:]
    {"download": lambda: download(*args), "process": lambda: process(*args),
     "validate": lambda: validate(args[0], args[1:]), "gencode": lambda: gencode(args[0], args[1:])}[cmd]()
