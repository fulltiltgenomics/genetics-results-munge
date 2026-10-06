"""Offline tests of `build_gene_sets.py`: the KEGG/BIOCARTA guard, the manifest writer, and the
rule that a refused or repeated build leaves no usable release directory. Downloads are replaced
by a fake `fetch` that writes tiny fixtures."""

import io
import json
import zipfile
from datetime import datetime, timezone

import pytest

import build_gene_sets as bgs

TODAY = datetime(2026, 10, 6, tzinfo=timezone.utc)
GOOD = "HALLMARK_X\thttp://u\t1\t2\nHALLMARK_Y\thttp://u\t3\n"


def msigdb_text(collection, n=2):
    prefix = bgs.MSIGDB_COLLECTIONS[collection]
    return "".join(f"{prefix}S{i}\thttps://www.gsea-msigdb.org/x/{prefix}S{i}\t1\t2\n"
                   for i in range(n))


def fake_download_factory(gmt_text=None, reactome_text="Metabolism\tR-HSA-1\tA\tB",
                          zip_extra=None):
    def fake_download(url, dest, **_):
        if url.endswith(".zip"):
            with zipfile.ZipFile(dest, "w") as z:
                z.writestr(bgs.REACTOME_GMT_MEMBER, reactome_text)
                for name in zip_extra or ():
                    z.writestr(name, "x")
        else:
            collection = url.rsplit("/", 1)[1].split(".v")[0]
            dest.write_text(gmt_text if gmt_text is not None else msigdb_text(collection))
        return dest
    return fake_download


def write(tmp_path, text):
    p = tmp_path / "a.gmt"
    p.write_text(text, encoding="utf-8")
    return p


def test_check_gmt_counts_sets(tmp_path):
    assert bgs.check_gmt(write(tmp_path, msigdb_text("h.all") + "\n"), "h.all") == 2


BYPASSES = ["KEGG_LEGACY_X", "KEGG_MEDICUS_Y", "BIOCARTA_Z", "kegg_lower", " KEGG_X",
            "\ufeffKEGG_X", "\xa0KEGG_X", "KEGG X", "KEGG-MEDICUS_X", "Biocarta pathway"]


@pytest.mark.parametrize("bad", BYPASSES)
def test_check_gmt_refuses_kegg_biocarta_msigdb(tmp_path, bad):
    p = write(tmp_path, msigdb_text("h.all") + f"{bad}\thttp://u\t1\n")
    with pytest.raises(SystemExit, match="refused source"):
        bgs.check_gmt(p, "h.all")


@pytest.mark.parametrize("bad", BYPASSES)
def test_check_gmt_refuses_kegg_biocarta_reactome(tmp_path, bad):
    p = write(tmp_path, f"Metabolism\tR-HSA-1\tA\n{bad}\tR-HSA-2\tA\n")
    with pytest.raises(SystemExit, match="refused source"):
        bgs.check_gmt(p, "reactome")


@pytest.mark.parametrize("collection", [*bgs.MSIGDB_COLLECTIONS, "reactome"])
def test_check_gmt_refuses_empty_name(tmp_path, collection):
    p = write(tmp_path, "\tKEGG_X\thttp://u\t1\n")
    with pytest.raises(SystemExit, match="empty set name"):
        bgs.check_gmt(p, collection)


@pytest.mark.parametrize("url", ["https://www.kegg.jp/pathway/hsa00010", "http://x/BIOCARTA_Y"])
def test_check_gmt_refuses_kegg_url_in_msigdb(tmp_path, url):
    p = write(tmp_path, f"HALLMARK_X\t{url}\t1\n")
    with pytest.raises(SystemExit, match="refused source"):
        bgs.check_gmt(p, "h.all")


@pytest.mark.parametrize("collection,name", [
    ("h.all", "GOBP_X"), ("c5.go.bp", "HALLMARK_X"), ("c5.hpo", "HP_lower"),
    ("c2.cp.pid", "PID_X Y"), ("c2.cp.wikipathways", "XWP_A"), ("c5.go.cc", " GOCC_X"),
])
def test_check_gmt_refuses_names_outside_collection_scheme(tmp_path, collection, name):
    p = write(tmp_path, f"{name}\thttp://u\t1\n")
    with pytest.raises(SystemExit, match="not a .* name"):
        bgs.check_gmt(p, collection)


@pytest.mark.parametrize("link", ["http://u", "R-HSA-", "R-MMU-1", " R-HSA-1", "R-HSA-1x"])
def test_check_gmt_reactome_requires_stable_id(tmp_path, link):
    p = write(tmp_path, f"Metabolism\t{link}\tA\n")
    with pytest.raises(SystemExit, match="not a Reactome stable id"):
        bgs.check_gmt(p, "reactome")


def test_check_gmt_refuses_malformed_and_empty(tmp_path):
    with pytest.raises(SystemExit, match="not a GMT line"):
        bgs.check_gmt(write(tmp_path, "onlyone\tfield\n"), "h.all")
    with pytest.raises(SystemExit, match="no gene sets"):
        bgs.check_gmt(write(tmp_path, "\n"), "h.all")


class FakeResponse(io.BytesIO):
    def __init__(self, body, content_length):
        super().__init__(body)
        self.headers = {} if content_length is None else {"Content-Length": str(content_length)}


@pytest.mark.parametrize("body,length,ok", [
    (b"12345", 10, False), (b"12345", 3, False), (b"12345", 5, True), (b"12345", None, True),
    (b"12345", "abc", False),
])
def test_download_checks_content_length(tmp_path, monkeypatch, body, length, ok):
    monkeypatch.setattr(bgs.urllib.request, "urlopen", lambda *a, **k: FakeResponse(body, length))
    dest = tmp_path / "f"
    if ok:
        assert bgs.download("http://x/f", dest).read_bytes() == body
    else:
        with pytest.raises(SystemExit, match="Content-Length"):
            bgs.download("http://x/f", dest)


def test_build_refuses_msigdb_gmt_without_trailing_newline(tmp_path, monkeypatch):
    monkeypatch.setattr(bgs, "download", fake_download_factory(gmt_text="HALLMARK_X\thttp://u\t1"))
    with pytest.raises(SystemExit, match="newline"):
        bgs.build(tmp_path, msigdb_release="2026.1.Hs", reactome="97", today=TODAY)
    assert list(tmp_path.iterdir()) == []


def test_build_refuses_releases_without_licence_entry(tmp_path, monkeypatch):
    def no_download(*_a, **_k):
        raise AssertionError("downloaded before the licence check")
    monkeypatch.setattr(bgs, "download", no_download)
    with pytest.raises(SystemExit, match="no committed licence entry for msigdb"):
        bgs.build(tmp_path, msigdb_release="2026.2.Hs", reactome="97", today=TODAY)
    with pytest.raises(SystemExit, match="no committed licence entry for reactome"):
        bgs.build(tmp_path, msigdb_release="2026.1.Hs", reactome="98", today=TODAY)
    monkeypatch.setattr(bgs, "reactome_release", lambda: "98")
    with pytest.raises(SystemExit, match="no committed licence entry for reactome"):
        bgs.build(tmp_path, msigdb_release="2026.1.Hs", reactome=None, today=TODAY)
    assert list(tmp_path.iterdir()) == []


def test_build_refuses_reactome_zip_with_extra_member(tmp_path, monkeypatch):
    monkeypatch.setattr(bgs, "download", fake_download_factory(zip_extra=["../evil.gmt"]))
    with pytest.raises(SystemExit, match="expected exactly"):
        bgs.build(tmp_path, msigdb_release="2026.1.Hs", reactome="97", today=TODAY)
    assert list(tmp_path.iterdir()) == []


def make_zip(tmp_path, members):
    z = tmp_path / "r.zip"
    with zipfile.ZipFile(z, "w") as archive:
        for name, text in members.items():
            archive.writestr(name, text)
    return z


def test_extract_member_refuses_missing_member(tmp_path):
    z = make_zip(tmp_path, {"Other.gmt": "x"})
    with pytest.raises(SystemExit, match="expected exactly"):
        bgs.extract_member(z, bgs.REACTOME_GMT_MEMBER, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_extract_member_size_cap(tmp_path):
    z = make_zip(tmp_path, {bgs.REACTOME_GMT_MEMBER: "0123456789"})
    with pytest.raises(SystemExit, match="over the 9 cap"):
        bgs.extract_member(z, bgs.REACTOME_GMT_MEMBER, tmp_path / "out", max_bytes=9)
    bgs.extract_member(z, bgs.REACTOME_GMT_MEMBER, tmp_path / "out", max_bytes=10)
    assert (tmp_path / "out").read_text() == "0123456789"


def test_build_refuses_leftover_partial(tmp_path, monkeypatch):
    monkeypatch.setattr(bgs, "download", fake_download_factory())
    partial = tmp_path / "gene_sets.20261006.partial"
    partial.mkdir()
    (partial / "keep").write_text("x")
    with pytest.raises(SystemExit, match="interrupted run"):
        bgs.build(tmp_path, msigdb_release="2026.1.Hs", reactome="97", today=TODAY)
    assert (partial / "keep").read_text() == "x"
    assert not (tmp_path / "gene_sets.20261006").exists()


def test_collection_list_excludes_bundles_and_kegg_sources():
    forbidden = {"c2.cp", "c2.all", "c2.cp.reactome", "c2.cp.kegg_legacy", "c2.cp.kegg_medicus",
                 "c2.cp.biocarta"}
    assert not forbidden & set(bgs.MSIGDB_COLLECTIONS)
    assert not any(c.startswith("msigdb") for c in bgs.MSIGDB_COLLECTIONS)


def test_build_writes_manifest(tmp_path, monkeypatch):
    monkeypatch.setattr(bgs, "download", fake_download_factory())
    out = bgs.build(tmp_path, msigdb_release="2026.1.Hs", reactome="97", today=TODAY)
    assert out.name == "gene_sets.20261006"
    manifest = json.loads((out / "manifest.json").read_text())
    by_name = {c["collection"]: c for c in manifest["collections"]}
    assert list(by_name) == [*bgs.MSIGDB_COLLECTIONS, "reactome"]
    hall = by_name["h.all"]
    assert hall["spdx_license"] == "CC-BY-4.0" and hall["n_sets"] == 2
    assert hall["url"].endswith("/2026.1.Hs/h.all.v2026.1.Hs.entrez.gmt")
    assert hall["sha256"] == bgs.sha256_of(out / hall["file"])
    assert "MSigDB" in hall["attribution"] and "Gene Ontology" not in hall["attribution"]
    go = by_name["c5.go.bp"]
    assert "Gene Ontology" in go["attribution"] and "2025-10-10" in go["attribution"]
    assert go["go_release"] in go["attribution"] and "no DOI" in go["attribution"]
    assert go["go_release_source"].startswith("https://docs.gsea-msigdb.org/")
    hpo = by_name["c5.hpo"]
    assert hpo["spdx_license"] == "LicenseRef-HPO" and hpo["hpo_release"] == "2026-01-08"
    assert "Human Phenotype Ontology Consortium" in hpo["attribution"]
    assert "2026-01-08" in hpo["attribution"] and "not altered" in hpo["attribution"]
    assert "Gene Ontology" not in hpo["attribution"]
    r = by_name["reactome"]
    assert r["spdx_license"] == "CC0-1.0" and r["upstream_release"] == "97"
    assert "/97/ReactomePathways.gmt.zip" in r["url"]
    assert r["gmt_sha256"] == bgs.sha256_of(out / "ReactomePathways.gmt")
    assert all(c["download_date"] == "2026-10-06" for c in manifest["collections"])
    assert not list(tmp_path.glob("*.partial"))


def test_build_refuses_existing_release(tmp_path, monkeypatch):
    monkeypatch.setattr(bgs, "download", fake_download_factory())
    existing = tmp_path / "gene_sets.20261006"
    existing.mkdir()
    (existing / "keep").write_text("x")
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        bgs.build(tmp_path, msigdb_release="2026.1.Hs", reactome="97", today=TODAY)
    assert (existing / "keep").read_text() == "x"


def test_build_guard_failure_leaves_no_directory(tmp_path, monkeypatch):
    monkeypatch.setattr(bgs, "download", fake_download_factory(reactome_text="KEGG_X\tR-HSA-1\tg\n"))
    with pytest.raises(SystemExit, match="refused source"):
        bgs.build(tmp_path, msigdb_release="2026.1.Hs", reactome="97", today=TODAY)
    assert list(tmp_path.iterdir()) == []


def test_every_licence_entry_is_complete():
    assert {source for source, _ in bgs.LICENCES} == {"msigdb", "reactome"}
    assert all(v["spdx"] and v["attribution"] for v in bgs.LICENCES.values())
    for (source, _), entry in bgs.LICENCES.items():
        if source == "msigdb":
            assert entry["go_release"] and entry["hpo_release"] and entry["release_notes"]
