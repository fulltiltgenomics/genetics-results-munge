"""aggregate.py against a synthetic CELLxGENE schema-7.1.0 H5AD and a brute-force reference.

Run: python -m pytest tests/cellxgene
"""
import json
import os
import sys

import h5py
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "cellxgene"))
import aggregate as agg  # noqa: E402

V3, VISIUM, SNMC, UNKNOWN = "EFO:0009922", "EFO:0022857", "EFO:0030027", "EFO:9999999"
N_CELLS, N_GENES = 400, 30
NONPRIMARY_ONLY, NEVER = 1, 2  # gene columns with special expression patterns
GENE_IDS = [f"ENSG{i:011d}" for i in range(N_GENES)]


def _strings(g, name, values):
    d = g.create_dataset(name, data=np.array(values, dtype=object), dtype=h5py.string_dtype())
    d.attrs.update({"encoding-type": "string-array", "encoding-version": "0.2.0"})


def _categorical(g, name, values):
    cats, codes = np.unique(np.asarray(values, dtype=str), return_inverse=True)
    c = g.create_group(name)
    c.attrs.update({"encoding-type": "categorical", "encoding-version": "0.2.0", "ordered": False})
    _array(c, "codes", codes.astype(np.int16))
    _strings(c, "categories", list(cats))


def _array(g, name, values):
    d = g.create_dataset(name, data=np.asarray(values))
    d.attrs.update({"encoding-type": "array", "encoding-version": "0.2.0"})


def _legacy_categorical(g, name, values):
    """anndata < 0.8: codes in the dataframe, categories behind an object reference."""
    cats, codes = np.unique(np.asarray(values, dtype=str), return_inverse=True)
    ref = g.require_group("__categories").create_dataset(
        name, data=np.array(list(cats), dtype=object), dtype=h5py.string_dtype())
    g.create_dataset(name, data=codes.astype(np.int16)).attrs["categories"] = ref.ref


def _dataframe(g, index, cols, plain=(), legacy=()):
    g.attrs.update({"encoding-type": "dataframe", "encoding-version": "0.2.0", "_index": "_index",
                    "column-order": np.array(list(cols), dtype=object)})
    _strings(g, "_index", index)
    for k, v in cols.items():
        if np.asarray(v).dtype == bool:
            _array(g, k, v)
        else:
            (_strings if k in plain else _legacy_categorical if k in legacy else _categorical)(g, k, list(v))


def _csr(g, name, dense):
    import scipy.sparse as sp
    m = sp.csr_matrix(dense)
    x = g.create_group(name)
    x.attrs.update({"encoding-type": "csr_matrix", "encoding-version": "0.1.0", "shape": m.shape})
    x.create_dataset("data", data=m.data.astype(np.float32))
    x.create_dataset("indices", data=m.indices.astype(np.int32))
    x.create_dataset("indptr", data=m.indptr.astype(np.int64))


def make_dataset(seed=0):
    rng = np.random.default_rng(seed)
    obs = {
        "tissue_ontology_term_id": rng.choice(["UBERON:0002048", "UBERON:0000178"], N_CELLS),
        "tissue_type": rng.choice(["tissue", "organoid"], N_CELLS),
        "cell_type_ontology_term_id": rng.choice(["CL:0000236", "CL:0000084", "CL:0000576"], N_CELLS),
        "disease_ontology_term_id": rng.choice(["PATO:0000461", "MONDO:0005015"], N_CELLS),
        "assay_ontology_term_id": rng.choice([V3, V3, VISIUM, SNMC, UNKNOWN], N_CELLS),
        "suspension_type": rng.choice(["cell", "nucleus"], N_CELLS),
        "donor_id": rng.choice([f"donor{i}" for i in range(5)], N_CELLS),
        "is_primary_data": rng.random(N_CELLS) < 0.8,
    }
    for k in list(obs):
        if k.endswith("_ontology_term_id"):
            obs[k.removesuffix("_ontology_term_id")] = np.char.add("label ", obs[k])
    counts = rng.poisson(0.7, (N_CELLS, N_GENES)).astype(np.float64)
    counts[:, NEVER] = 0
    counts[:, NONPRIMARY_ONLY] = np.where(obs["is_primary_data"], 0, 50)
    counts[0] = 0  # a cell with no counts at all must not divide by zero
    return obs, counts


def write_h5ad(path, obs, counts, normalised_x=False, feature_is_filtered=None, plain=(), legacy=()):
    with h5py.File(path, "w") as f:
        f.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})
        _dataframe(f.create_group("obs"), [f"cell{i}" for i in range(N_CELLS)], obs, plain, legacy)
        filtered = np.zeros(N_GENES, bool) if feature_is_filtered is None else feature_is_filtered
        var = {"feature_name": [f"GENE{i}-α" for i in range(N_GENES)], "feature_is_filtered": filtered}
        _dataframe(f.create_group("var"), GENE_IDS, var)
        if normalised_x:
            lib = np.maximum(counts.sum(1, keepdims=True), 1)
            _csr(f, "X", np.log1p(counts / lib * 1e4))
        else:
            _csr(f, "X", counts)
        if normalised_x == "with_raw":
            raw = f.create_group("raw")
            raw.attrs.update({"encoding-type": "raw", "encoding-version": "0.1.0"})
            _csr(raw, "X", counts)
            _dataframe(raw.create_group("var"), GENE_IDS, {"feature_name": var["feature_name"]})
    return path


def expected(obs, counts):
    keep = obs["is_primary_data"] & np.isin(obs["assay_ontology_term_id"], list(agg.ASSAYS))
    lib = counts.sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        logcp = np.where(lib[:, None] > 0, np.log1p(counts * 1e4 / lib[:, None]), 0.0)
    keys = pd.DataFrame({k: obs[k] for k in agg.GROUP_KEYS})[keep]
    groups, gd = {}, {}
    for key, idx in keys.groupby(agg.GROUP_KEYS).groups.items():
        idx = np.asarray(idx)
        groups[key] = {"n_cells": len(idx), "sum": logcp[idx].sum(0), "n_expr": (counts[idx] > 0).sum(0),
                       "n_donors": len(set(obs["donor_id"][idx]))}
        for d in set(obs["donor_id"][idx]):
            di = idx[obs["donor_id"][idx] == d]
            gd[key + (d,)] = {"counts": counts[di].sum(0), "lib": lib[di].sum(), "n_cells": len(di)}
    return groups, gd


def read_out(out):
    groups = pd.read_parquet(f"{out}/groups.parquet")
    keyed = {tuple(r[k] for k in agg.GROUP_KEYS): r for r in groups.to_dict("records")}
    gg = pd.read_parquet(f"{out}/group_gene.parquet")
    gdon = pd.read_parquet(f"{out}/group_donor.parquet")
    gdg = pd.read_parquet(f"{out}/group_donor_gene.parquet")
    return groups, keyed, gg, gdon, gdg, json.load(open(f"{out}/_SUCCESS.json"))


def run(tmp_path, h5ad, *extra):
    out = tmp_path / "out"
    agg.main(["--input", str(h5ad), "--out", str(out), "--dataset-id", "ds1",
              "--dataset-version-id", "dv1", *extra])
    return out


def assert_matches(out, obs, counts):
    exp_g, exp_gd = expected(obs, counts)
    groups, keyed, gg, gdon, gdg, _ = read_out(out)
    assert set(keyed) == set(exp_g)
    for key, e in exp_g.items():
        g = keyed[key]
        assert g["n_cells"] == e["n_cells"] and g["n_donors"] == e["n_donors"]
        rows = gg[gg.group == g["group"]].set_index("gene")
        # sparse: exactly the genes with at least one expressing cell
        assert set(rows.index) == set(np.flatnonzero(e["n_expr"]))
        np.testing.assert_array_equal(rows.n_expressing, e["n_expr"][rows.index])
        np.testing.assert_allclose(rows.sum_log1p_cp10k, e["sum"][rows.index], rtol=1e-12)
        np.testing.assert_allclose(rows.mean_log1p_cp10k, e["sum"][rows.index] / e["n_cells"], rtol=1e-6)
    gdon = gdon.merge(groups[["group", *agg.GROUP_KEYS]], on="group")
    assert len(gdon) == len(exp_gd)
    for r in gdon.to_dict("records"):
        e = exp_gd[tuple(r[k] for k in agg.GROUP_KEYS) + (r["donor_id"],)]
        assert r["n_cells"] == e["n_cells"] and r["library_size"] == e["lib"]
        got = np.zeros(N_GENES)
        sel = gdg[gdg.gd == r["gd"]]
        got[sel.gene] = sel.counts
        np.testing.assert_array_equal(got, e["counts"])
    assert (gdg.counts > 0).all()


@pytest.fixture
def data():
    return make_dataset()


@pytest.mark.parametrize("extra", [[], ["--chunk-nnz", "97", "--max-accum-bytes", "400"]],
                         ids=["one-chunk", "many-chunks-and-gene-blocks"])
def test_matches_bruteforce(tmp_path, data, extra):
    obs, counts = data
    out = run(tmp_path, write_h5ad(tmp_path / "d.h5ad", obs, counts), *extra)
    assert_matches(out, obs, counts)
    succ = read_out(out)[-1]
    assert succ["matrix"] == "X"
    if extra:
        assert succ["gene_blocks"] > 1


def test_raw_x_preferred_over_normalised_x(tmp_path, data):
    obs, counts = data
    out = run(tmp_path, write_h5ad(tmp_path / "d.h5ad", obs, counts, normalised_x="with_raw"))
    assert read_out(out)[-1]["matrix"] == "raw/X"
    assert_matches(out, obs, counts)


def test_non_integer_counts_rejected(tmp_path, data):
    obs, counts = data
    h5ad = write_h5ad(tmp_path / "d.h5ad", obs, counts, normalised_x=True)
    out = tmp_path / "out"
    out.mkdir()
    (out / "_SUCCESS.json").write_text("{}")  # a previous run's marker must not survive
    with pytest.raises(SystemExit) as e:
        run(tmp_path, h5ad)
    assert e.value.code == 3
    assert not (out / "_SUCCESS.json").exists()


def test_is_primary_data_filter(tmp_path, data):
    obs, counts = data
    out = run(tmp_path, write_h5ad(tmp_path / "d.h5ad", obs, counts))
    groups, _, gg, gdon, gdg, succ = read_out(out)
    allowed = np.isin(obs["assay_ontology_term_id"], list(agg.ASSAYS))
    assert succ["cells_not_primary"] == int((~obs["is_primary_data"]).sum())
    assert succ["cells_kept"] == int((obs["is_primary_data"] & allowed).sum()) == groups.n_cells.sum()
    # expressed in non-primary cells only
    assert NONPRIMARY_ONLY not in set(gg.gene) and NONPRIMARY_ONLY not in set(gdg.gene)


def test_assay_allow_list(tmp_path, data):
    obs, counts = data
    out = run(tmp_path, write_h5ad(tmp_path / "d.h5ad", obs, counts))
    groups, *_, succ = read_out(out)
    assert set(groups.assay_ontology_term_id) == {V3, VISIUM}
    assert groups.spatial.eq(groups.assay_ontology_term_id == VISIUM).all()
    prim = obs["is_primary_data"]
    assert succ["cells_assay_excluded"] == {
        a: int((prim & (obs["assay_ontology_term_id"] == a)).sum()) for a in (SNMC, UNKNOWN)}
    assert SNMC not in agg.ASSAYS and VISIUM in agg.ASSAYS and "EFO:0030062" in agg.SPATIAL_ASSAYS


def test_labels_and_keys(tmp_path, data):
    obs, counts = data
    groups = read_out(run(tmp_path, write_h5ad(tmp_path / "d.h5ad", obs, counts)))[0]
    for k, lab in agg.LABEL_OF.items():
        assert (groups[lab] == "label " + groups[k]).all()
    assert set(groups.tissue_type) == {"tissue", "organoid"}


def test_unmeasured_is_na_and_unexpressed_is_zero(tmp_path, data):
    obs, counts = data
    filtered = np.zeros(N_GENES, bool)
    filtered[3] = True  # zeroed out of X: not measured in this dataset
    counts = counts.copy()
    counts[:, 3] = 0
    out = run(tmp_path, write_h5ad(tmp_path / "d.h5ad", obs, counts, feature_is_filtered=filtered))
    query = GENE_IDS + ["ENSG99999999999"]
    groups, m = agg.expression_matrix(out, query)
    assert np.isnan(m[:, -1]).all()  # absent from var
    assert np.isnan(m[:, 3]).all()  # in var, feature_is_filtered
    assert (m[:, NEVER] == 0).all() and (m[:, NONPRIMARY_ONLY] == 0).all()  # measured, never expressed
    assert not np.isnan(np.delete(m, [3, len(query) - 1], axis=1)).any()
    gene0 = pd.read_parquet(f"{out}/group_gene.parquet").query("gene == 0").set_index("group")
    np.testing.assert_allclose(m[gene0.index, 0], gene0.mean_log1p_cp10k)


def test_success_json(tmp_path, data):
    obs, counts = data
    out = run(tmp_path, write_h5ad(tmp_path / "d.h5ad", obs, counts))
    succ = read_out(out)[-1]
    assert succ["aggregation_format_version"] == agg.FORMAT_VERSION
    assert (succ["dataset_id"], succ["dataset_version_id"]) == ("ds1", "dv1")
    assert succ["source"]["bytes"] == os.path.getsize(tmp_path / "d.h5ad")
    assert succ["nnz"] == int((counts > 0).sum())
    assert {"read_obs_var", "main_pass", "write"} <= set(succ["phases"])
    assert all(p["peak_rss_mb"] > 0 for p in succ["phases"].values())
    assert succ["outputs"] == {n: os.path.getsize(out / n) for n in agg.OUTPUT_FILES}
    assert succ["peak_rss_mb"] > 0 and succ["n_genes_measured"] == N_GENES


def test_gcs_output_marker_last(tmp_path, data, monkeypatch):
    obs, counts = data
    events = []

    class Blob:
        def __init__(self, name):
            self.name = name

        def exists(self):
            return self.name.endswith("_SUCCESS.json")

        def delete(self):
            events.append(("delete", self.name))

        def upload_from_filename(self, path):
            assert os.path.exists(path)
            events.append(("upload", self.name))

    class Bucket:
        blob = staticmethod(Blob)

    monkeypatch.setattr(agg, "_gcs_bucket", lambda uri: (Bucket(), "pre/fix"))
    agg.main(["--input", str(write_h5ad(tmp_path / "d.h5ad", obs, counts)), "--out", "gs://b/pre/fix",
              "--dataset-id", "ds1", "--dataset-version-id", "dv1", "--workdir", str(tmp_path)])
    assert events[0] == ("delete", "pre/fix/_SUCCESS.json")
    assert events[-1] == ("upload", "pre/fix/_SUCCESS.json")
    assert {n for _, n in events[1:-1]} == {f"pre/fix/{n}" for n in agg.OUTPUT_FILES}
    assert not [p for p in os.listdir(tmp_path) if p.startswith("cxg-")]  # staging cleaned up


def test_fixture_is_a_valid_h5ad(tmp_path, data):
    anndata = pytest.importorskip("anndata")
    obs, counts = data
    ad = anndata.read_h5ad(write_h5ad(tmp_path / "d.h5ad", obs, counts, normalised_x="with_raw"))
    assert ad.shape == (N_CELLS, N_GENES) and ad.raw.X.sum() == counts.sum()
    assert ad.obs.is_primary_data.dtype == bool


def _non_ascii(obs):
    obs = dict(obs)
    obs["disease"] = np.where(obs["disease_ontology_term_id"] == "MONDO:0005015", "Sjögren syndrome", "normal")
    obs["cell_type"] = np.char.add("β ", obs["cell_type_ontology_term_id"])
    obs["donor_id"] = np.char.add("dönor-", obs["donor_id"])
    return obs


@pytest.mark.parametrize("enc", [{}, {"plain": ("donor_id", "cell_type")}, {"legacy": ("donor_id", "disease")}],
                         ids=["categorical", "plain-string-array", "legacy-categories-ref"])
def test_non_ascii_strings(tmp_path, data, enc):
    obs, counts = data
    obs = _non_ascii(obs)
    out = run(tmp_path, write_h5ad(tmp_path / "d.h5ad", obs, counts, **enc))
    assert_matches(out, obs, counts)
    groups = read_out(out)[0]
    assert "Sjögren syndrome" in set(groups.disease) and groups.cell_type.str.startswith("β ").all()
    assert pd.read_parquet(f"{out}/group_donor.parquet").donor_id.str.startswith("dönor-").all()
    assert pd.read_parquet(f"{out}/genes.parquet").feature_name.str.endswith("-α").all()


@pytest.mark.parametrize("missing", ["donor_id", "is_primary_data", "cell_type"])
def test_missing_mandatory_column_rejected(tmp_path, data, missing):
    obs, counts = data
    obs = {k: v for k, v in obs.items() if k != missing}
    with pytest.raises(SystemExit) as e:
        run(tmp_path, write_h5ad(tmp_path / "d.h5ad", obs, counts))
    assert e.value.code == 3


@pytest.mark.parametrize("case", ["no-primary", "only-excluded-assays"])
def test_no_kept_cells_is_an_empty_result(tmp_path, data, case):
    obs, counts = dict(data[0]), data[1]
    if case == "no-primary":
        obs["is_primary_data"] = np.zeros(N_CELLS, bool)
    else:
        obs["assay_ontology_term_id"] = np.where(np.arange(N_CELLS) % 2, SNMC, UNKNOWN)
    out = run(tmp_path, write_h5ad(tmp_path / "d.h5ad", obs, counts))  # exits 0: no SystemExit
    groups, _, gg, gdon, gdg, succ = read_out(out)
    assert succ["cells_kept"] == 0 and succ["n_groups"] == 0 and succ["gene_blocks"] == 0
    assert len(groups) == len(gg) == len(gdon) == len(gdg) == 0
    assert set(agg.GROUP_KEYS) <= set(groups.columns) and "mean_log1p_cp10k" in gg.columns
    assert len(pd.read_parquet(f"{out}/genes.parquet")) == N_GENES


def test_group_accumulators_count_toward_the_budget():
    # a ~6k-group atlas with no donor accumulator at all must still be split
    blocks = agg.gene_blocks(6000, 0, 60000, agg.MAX_ACCUM_BYTES)
    assert len(blocks) > 1
    assert max(hi - lo for lo, hi in blocks) * 12 * 6000 <= agg.MAX_ACCUM_BYTES
    assert blocks[0][0] == 0 and blocks[-1][1] == 60000
    assert all(a[1] == b[0] for a, b in zip(blocks, blocks[1:]))


def test_gene_blocks_triggered_by_groups_end_to_end(tmp_path, data):
    obs, counts = dict(data[0]), data[1]
    obs["donor_id"] = np.full(N_CELLS, "donor0")  # one donor, so n_group_donor == n_groups
    n_groups = len(expected(obs, counts)[0])
    # the donor accumulator alone (4 bytes x pairs x genes) fits; adding the group ones does not
    budget = 16 * n_groups * 10
    assert 4 * n_groups * N_GENES <= budget < 16 * n_groups * N_GENES
    out = run(tmp_path, write_h5ad(tmp_path / "d.h5ad", obs, counts), "--max-accum-bytes", str(budget))
    assert read_out(out)[-1]["gene_blocks"] == 3
    assert_matches(out, obs, counts)
