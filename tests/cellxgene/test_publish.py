"""publish.py on synthetic worker outputs: three datasets whose expected pooled values are
worked out by hand below."""
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "cellxgene"))
import aggregate  # noqa: E402
import publish  # noqa: E402

pytestmark = pytest.mark.skipif(not (shutil.which("tabix") and shutil.which("bgzip")),
                                reason="tabix/bgzip not on PATH")

ROLLUP = os.path.join(os.path.dirname(publish.__file__), "rollup")
NORMAL = "PATO:0000461"
FIB, MYOFIB, TCELL = "CL:0000057", "CL:0000186", "CL:0000084"
G1, G2, G3 = "ENSG00000000001", "ENSG00000000002", "ENSG00000000003"


def _tsv(name):
    with open(os.path.join(ROLLUP, name)) as fh:
        return list(csv.DictReader((x for x in fh if not x.startswith("# ")), delimiter="\t"))


# a native tissue that the shipped table rolls up to a different tissue_general
TISSUE = next(r for r in _tsv("tissue_rollup.tsv")
              if r["method"] == "tissue_general" and r["tissue_ontology_term_id"] != r["tissue_general_ontology_term_id"])


def write_dataset(root, dvid, groups, donors, gene_ids, gg, gdg, measured=None):
    d = root / dvid
    d.mkdir(parents=True)
    rows = []
    for i, (ct, disease, n_cells) in enumerate(groups):
        rows.append({"group": i, "tissue_ontology_term_id": TISSUE["tissue_ontology_term_id"], "tissue_type": "tissue",
                     "cell_type_ontology_term_id": ct, "disease_ontology_term_id": disease,
                     "assay_ontology_term_id": "EFO:0009922", "suspension_type": "cell",
                     "tissue": TISSUE["tissue"], "cell_type": ct, "disease": disease, "assay": "10x 3' v3",
                     "spatial": False, "n_cells": n_cells, "n_donors": len(donors[i])})
    pd.DataFrame(rows).astype({"group": "int32", "n_cells": "int32", "n_donors": "int32"}).to_parquet(d / "groups.parquet")
    gd = [(g, don, lib) for g, ds in donors.items() for don, lib in ds]
    pd.DataFrame({"gd": range(len(gd)), "group": [x[0] for x in gd], "donor_id": [x[1] for x in gd],
                  "n_cells": 1, "library_size": [float(x[2]) for x in gd]}).astype(
        {"gd": "int32", "group": "int32", "n_cells": "int32"}).to_parquet(d / "group_donor.parquet")
    gix = {g: i for i, g in enumerate(gene_ids)}
    pd.DataFrame({"gene": range(len(gene_ids)), "gene_id": gene_ids,
                  "measured": measured or [True] * len(gene_ids)}).astype({"gene": "int32"}).to_parquet(d / "genes.parquet")
    pd.DataFrame([{"group": g, "gene": gix[gene], "n_expressing": n, "sum_log1p_cp10k": float(s),
                   "mean_log1p_cp10k": s / groups[g][2]} for g, gene, n, s in gg]).astype(
        {"group": "int32", "gene": "int32", "n_expressing": "int32", "mean_log1p_cp10k": "float32"}).to_parquet(
        d / "group_gene.parquet")
    pd.DataFrame([{"gd": gdi, "gene": gix[gene], "counts": float(c)} for gdi, gene, c in gdg],
                 columns=["gd", "gene", "counts"]).astype(
        {"gd": "int32", "gene": "int32", "counts": "float32"}).to_parquet(d / "group_donor_gene.parquet")
    (d / "_SUCCESS.json").write_text(json.dumps({"aggregation_format_version": aggregate.FORMAT_VERSION,
                                                 "dataset_id": "id-" + dvid, "dataset_version_id": dvid,
                                                 "cells_kept": sum(g[2] for g in groups), "n_groups": len(groups)}))


@pytest.fixture
def inputs(tmp_path):
    agg = tmp_path / "agg"
    a_donors = [("a1", 1000), ("a2", 1000), ("a3", 1000)]
    # A: myofibroblast (pinned to fibroblast) and fibroblast pool into one summary row; the
    # disease group is per-dataset only; the 40-cell T cell group is below min_cells
    write_dataset(agg, "dsA",
                  [(MYOFIB, NORMAL, 100), (FIB, NORMAL, 60), (FIB, "MONDO:0005002", 80), (TCELL, NORMAL, 40)],
                  {0: a_donors, 1: [(n, 2 * lib) for n, lib in a_donors],
                   2: a_donors, 3: a_donors},
                  [G1 + ".5", G2, G3, "ENSG00000099999", "ERCC-00002"],
                  [(0, G1 + ".5", 50, 100), (0, G2, 10, 20), (1, G1 + ".5", 30, 60), (2, G1 + ".5", 40, 40),
                   (3, G1 + ".5", 20, 40)],
                  [(0, G1 + ".5", 10), (1, G1 + ".5", 20), (2, G1 + ".5", 30), (0, G2, 5),
                   (3, G1 + ".5", 6), (4, G1 + ".5", 6), (5, G1 + ".5", 6)])
    # B: does not have G2 in its var at all (NA), measures G3 but never sees it (0); 2 donors
    write_dataset(agg, "dsB", [(FIB, NORMAL, 200)], {0: [("b1", 2000), ("b2", 2000)]}, [G1, G3],
                  [(0, G1, 100, 50)], [(0, G1, 4), (1, G1, 8)])
    # C: G2 is in var but measured false (filtered), its stray row must be ignored; G3 absent
    write_dataset(agg, "dsC", [(FIB, NORMAL, 100)], {0: [("c1", 1000), ("c2", 1000), ("c3", 1000)]}, [G1, G2],
                  [(0, G1, 10, 10), (0, G2, 5, 500)], [], measured=[True, False])
    gencode = tmp_path / "gencode.tsv"
    gencode.write_text("gene_id\tchrom\tgene_start\tgene_end\tgene_strand\tgene_name\tgene_type\n"
                       f"{G1}.7\t1\t1000\t1500\t+\tGENEA\tprotein_coding\n"
                       f"{G2}.2\t1\t5000\t5600\t+\tGENEB\tprotein_coding\n"
                       f"{G3}.1\t23\t2000\t2600\t-\tGENEC\tprotein_coding\n")
    return tmp_path, agg, gencode


def run(tmp_path, agg, gencode, out="rel", *extra):
    publish.main(["--agg", str(agg), "--out", str(tmp_path / out), "--release-date", "2026-01-01",
                  "--gencode", str(gencode), "--workdir", str(tmp_path / ("work-" + out)), *extra])
    return tmp_path / out / "2026-01-01"


def read(path):
    return pd.read_csv(path, sep="\t", compression="gzip", dtype={"chrom": str})


def summary_row(df, gene, ct=FIB):
    r = df[(df["#gene_id"].str.startswith(gene)) & (df.cell_type_rollup_ontology_term_id == ct)]
    assert len(r) == 1, r
    return r.iloc[0]


def test_pooling_na_versus_zero_and_rollup(inputs):
    rel = run(*inputs)
    s = read(rel / "summary.tsv.gz")
    assert set(s.tissue_general_ontology_term_id) == {TISSUE["tissue_general_ontology_term_id"]}
    # G1: measured by A (100+60 normal cells), B (200), C (100)
    r = summary_row(s, G1)
    assert (r.n_cells, r.n_donors, r.n_datasets) == (460, 8, 3)
    assert r.mean_log1p_cp10k == pytest.approx(round(220 / 460, 4))
    assert r.frac_expressing == pytest.approx(round(190 / 460, 4))
    # G2: B never had it, C filtered it -> only A's 160 cells in the denominator, not 460
    r = summary_row(s, G2)
    assert (r.n_cells, r.n_datasets, r.mean_log1p_cp10k) == (160, 1, 0.125)
    # G3: measured by A and B, never seen -> a real zero over 360 cells
    r = summary_row(s, G3)
    assert (r.n_cells, r.n_datasets, r.mean_log1p_cp10k, r.frac_expressing) == (360, 2, 0, 0)
    # the pin sends myofibroblast to fibroblast: one rolled-up cell type in the summary
    assert set(s.cell_type_rollup_ontology_term_id) == {FIB}


def test_pseudobulk(inputs):
    rel = run(*inputs)
    r = summary_row(read(rel / "summary.tsv.gz"), G1)
    # donor CP10k pooled within donor across A's two groups (library 1000 + 2000), C's donors 0
    assert r.donor_pseudobulk_n_donors == 8
    assert r.donor_pseudobulk_mean_cp10k == pytest.approx((1e4 * (16 + 26 + 36) / 3000 + 20 + 40) / 8)
    p = read(rel / "per_dataset.tsv.gz")
    a0 = p[(p.dataset_id == "id-dsA") & (p.cell_type_ontology_term_id == MYOFIB) & p["#gene_id"].str.startswith(G1)]
    assert a0.donor_pseudobulk_mean_cp10k.item() == pytest.approx(200)
    assert a0.cell_type_rollup_ontology_term_id.item() == FIB


def test_pseudobulk_needs_three_donors(inputs):
    rel = run(*inputs, "rel1", "--min-cells", "1", "--min-donors", "1")
    p = read(rel / "per_dataset.tsv.gz")
    b = p[p.dataset_id == "id-dsB"]
    assert len(b) == 1 and b.mean_log1p_cp10k.item() == 0.25 and pd.isna(b.donor_pseudobulk_mean_cp10k.item())
    rel = run(*inputs, "rel2", "--min-cells", "1", "--min-donors", "1", "--min-pseudobulk-donors", "2")
    b = read(rel / "per_dataset.tsv.gz").query("dataset_id == 'id-dsB'")
    assert b.donor_pseudobulk_mean_cp10k.item() == pytest.approx(30)


def test_suppression(inputs):
    rel = run(*inputs)
    p = read(rel / "per_dataset.tsv.gz")
    # 40-cell T cell group and B's 2-donor group are not emitted; the disease group is, per dataset
    assert TCELL not in set(p.cell_type_ontology_term_id)
    assert "id-dsB" not in set(p.dataset_id)
    assert "MONDO:0005002" in set(p.disease_ontology_term_id)
    s = read(rel / "summary.tsv.gz")
    assert TCELL not in set(s.cell_type_rollup_ontology_term_id)
    # tau needs two cell types passing suppression; fibroblast alone gives NA
    assert s.tau.isna().all()
    s1 = read(run(*inputs, "rel1", "--min-cells", "1", "--min-donors", "1") / "summary.tsv.gz")
    assert summary_row(s1, G1, TCELL).mean_log1p_cp10k == 1.0
    assert summary_row(s1, G1).tau == pytest.approx(round(1 - (220 / 460) / 1.0, 4))
    m = json.loads((rel / "manifest.json").read_text())
    assert m["thresholds"] == {"min_cells": 50, "min_donors": 3, "min_pseudobulk_donors": 3,
                               "summary_include_spatial": False}


def test_gencode_drop_counting(inputs):
    m = json.loads((run(*inputs) / "manifest.json").read_text())
    assert m["gencode"]["dropped_no_v49_match"] == 1
    assert m["gencode"]["dropped_ids"] == ["ENSG00000099999"]
    assert m["gencode"]["dropped_non_ensembl_features"] == 1
    assert m["input_dataset_version_ids"] == ["dsA", "dsB", "dsC"]
    p = read(run(*inputs, "rel1") / "per_dataset.tsv.gz")
    # v49 id and coordinates, the source's own versioned id kept beside them
    a = p[(p.dataset_id == "id-dsA") & p["#gene_id"].str.startswith(G1)]
    assert set(a["#gene_id"]) == {G1 + ".7"} and set(a.source_gene_id) == {G1 + ".5"}


def test_tabix_point_query(inputs):
    rel = run(*inputs)
    for name in ("summary", "per_dataset"):
        out = subprocess.run(["tabix", str(rel / f"{name}.tsv.gz"), "1:1000-1000"], capture_output=True,
                             text=True, check=True).stdout.splitlines()
        assert out and all(line.split("\t")[4] == "GENEA" for line in out)
    x = subprocess.run(["tabix", str(rel / "summary.tsv.gz"), "23:2000-2000"], capture_output=True,
                       text=True, check=True).stdout
    assert x.split("\t")[4] == "GENEC"


def test_deterministic_and_genes_mode(inputs):
    a, b = run(*inputs, "r1"), run(*inputs, "r2")
    for root, _, files in os.walk(a):
        for f in files:
            p = os.path.join(root, f)
            with open(p, "rb") as x, open(os.path.join(b, os.path.relpath(p, a)), "rb") as y:
                assert x.read() == y.read(), f
    rg = run(*inputs, "rg", "--genes", "GENEA")
    assert not rg.exists()  # nothing under <out>/<date>/: genes mode is its sibling
    g = rg.parent / (rg.name + "-genes-mode")
    full = read(a / "summary.tsv.gz")
    assert read(g / "summary.tsv.gz").equals(full[full.gene_name == "GENEA"].reset_index(drop=True))
    shards = json.loads((a / "manifest.json").read_text())["bq_shards"]
    for kind, cols in (("per_dataset", publish.BQ_PER_DATASET_COLUMNS), ("summary", publish.SUMMARY_COLUMNS)):
        assert shards[kind] and sum(shards[kind].values()) == json.loads((a / "manifest.json").read_text())["rows"][kind]
        for name, n in shards[kind].items():
            with gzip.open(a / name, "rt") as fh:
                assert fh.readline().rstrip("\n").split("\t") == cols
                assert sum(1 for _ in fh) == n


def test_group_keys_unique_with_noncontiguous_ids(inputs):
    tmp_path, agg, gencode = inputs
    # ids with gaps, as when groups were filtered out after numbering: dsA's 0,5,10,15 against
    # an offset of its row count (4) would give dsC's group 0 the key 5, already dsA's
    d = agg / "dsA"
    for f in ("groups.parquet", "group_donor.parquet", "group_gene.parquet"):
        x = pd.read_parquet(d / f)
        x["group"] = (x["group"] * 5).astype("int32")
        x.to_parquet(d / f)
    rel = run(tmp_path, agg, gencode, "rel", "--min-cells", "1", "--min-donors", "1")
    g = pd.read_csv(rel / "bq" / "groups.tsv.gz", sep="\t")
    assert g.group_key.is_unique and len(g) == 6
    p = pd.concat(pd.read_csv(f, sep="\t") for f in sorted((rel / "bq" / "per_dataset").glob("*.tsv.gz")))
    assert not p.duplicated(["gene_id", "group_key"]).any()
    assert (p.merge(g, on="group_key").dataset_key_x == p.merge(g, on="group_key").dataset_key_y).all()


def test_fresh_release_dir_and_checksums(inputs):
    rel = run(*inputs)
    (rel / "stale.txt").write_text("left over")
    rel = run(*inputs)
    assert not (rel / "stale.txt").exists()
    m = json.loads((rel / "manifest.json").read_text())
    assert set(m["files"]) == {os.path.relpath(os.path.join(r, f), rel) for r, _, fs in os.walk(rel)
                               for f in fs} - {"manifest.json"}
    assert m["files"]["summary.tsv.gz"]["sha256"] == publish.sha256(rel / "summary.tsv.gz")


def test_genes_mode_refuses_remote_release(inputs):
    tmp_path, agg, gencode = inputs
    with pytest.raises(SystemExit, match="genes-mode|partial release"):
        publish.main(["--agg", str(agg), "--out", "gs://nowhere/releases", "--gencode", str(gencode),
                      "--workdir", str(tmp_path / "w"), "--genes", "GENEA"])


def test_measured_genes_list(inputs):
    rel = run(*inputs)
    m = pd.read_csv(rel / "measured_genes.tsv.gz", sep="\t")
    got = {d: set(x.str.split(".").str[0]) for d, x in m.groupby("dataset_version_id").gene_id}
    # B has no G2 in var, C has it with measured false; unmapped and spike-in features never listed
    assert got == {"dsA": {G1, G2, G3}, "dsB": {G1, G3}, "dsC": {G1}}
    manifest = json.loads((rel / "manifest.json").read_text())
    assert "measured_genes.tsv.gz" in manifest["files"]
    # the BigQuery copy: same pairs, keyed by dataset_key, with the gene's chrom and symbol
    bq = pd.read_csv(rel / "bq" / "measured_genes.tsv.gz", sep="\t")
    assert list(bq.columns) == publish.BQ_MEASURED_COLUMNS
    assert len(bq) == len(m) == manifest["rows"]["measured_genes"]
    keys = pd.read_csv(rel / "bq" / "datasets.tsv", sep="\t").set_index("dataset_key").dataset_version_id
    assert set(zip(bq.dataset_key.map(keys), bq.gene_id)) == set(zip(m.dataset_version_id, m.gene_id))
    per = pd.concat(pd.read_csv(f, sep="\t") for f in sorted((rel / "bq" / "per_dataset").glob("*.tsv.gz")))
    coords = per.drop_duplicates("gene_id").set_index("gene_id")[["chrom", "gene_start", "gene_name"]]
    joined = bq.join(coords, on="gene_id", rsuffix="_pd", how="inner")
    assert len(joined) and (joined.chrom == joined.chrom_pd).all() and (joined.gene_name == joined.gene_name_pd).all()


class FakeBlob:
    def __init__(self, store, name):
        self.store, self.name = store, name

    def download_to_filename(self, path):
        shutil.copyfile(self.store.root / self.name, path)
        self.store.downloaded.append(self.name)

    def upload_from_filename(self, path):
        self.store.uploaded.append(self.name)


class FakeBucket:
    def __init__(self, root):
        self.root, self.downloaded, self.uploaded = root, [], []

    def blob(self, name):
        return FakeBlob(self, name)

    def list_blobs(self, prefix):
        return [FakeBlob(self, str(p.relative_to(self.root))) for p in sorted(self.root.rglob("*"))
                if p.is_file() and str(p.relative_to(self.root)).startswith(prefix)]


def test_gcs_without_gcloud(inputs, monkeypatch):
    tmp_path, agg, gencode = inputs
    root = tmp_path / "bucket"
    shutil.copytree(agg, root / "agg")
    shutil.copytree(agg / "dsC", root / "agg" / "dsOLD")  # a superseded version left in the prefix
    shutil.copy(gencode, root / "gencode.tsv")
    bucket = FakeBucket(root)
    google, cloud, storage = (type(sys)("google"), type(sys)("google.cloud"), type(sys)("google.cloud.storage"))
    storage.Client = lambda: type("C", (), {"bucket": lambda self, name: bucket})()
    google.cloud, cloud.storage = cloud, storage
    monkeypatch.setitem(sys.modules, "google", google)
    monkeypatch.setitem(sys.modules, "google.cloud", cloud)
    monkeypatch.setitem(sys.modules, "google.cloud.storage", storage)
    real_popen = subprocess.Popen

    def no_gcloud(cmd, *args, **kw):
        assert os.path.basename(str(cmd[0])) != "gcloud", cmd
        return real_popen(cmd, *args, **kw)
    monkeypatch.setattr(subprocess, "Popen", no_gcloud)
    monkeypatch.setattr(shutil, "which", lambda name, *a, **k: None if name == "gcloud" else "/bin/true")
    lst = tmp_path / "datasets.txt"
    lst.write_text("dsA\ndsB\n")
    publish.main(["--agg", "gs://b/agg", "--datasets", "dsC", "--datasets-file", str(lst),
                  "--out", "gs://b/releases", "--release-date", "2026-10-14-120000",
                  "--gencode", "gs://b/gencode.tsv", "--workdir", str(tmp_path / "w")])
    assert not any("dsOLD" in n for n in bucket.downloaded)
    assert bucket.uploaded[-1] == "releases/2026-10-14-120000/manifest.json"
    assert "releases/2026-10-14-120000/summary.tsv.gz" in bucket.uploaded
    with pytest.raises(SystemExit, match="--datasets"):
        publish.main(["--agg", "gs://b/agg", "--out", str(tmp_path / "o"), "--gencode", str(gencode),
                      "--workdir", str(tmp_path / "w2")])
