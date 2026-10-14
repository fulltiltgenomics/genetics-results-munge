"""plan.py against a recorded slice of the curation API listing (listing_fixture.json).

Run: python -m pytest tests/cellxgene/test_plan.py
"""
import copy
import json
import os
import sys

import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts", "cellxgene"))
import aggregate as agg  # noqa: E402
import plan  # noqa: E402

with open(os.path.join(HERE, "listing_fixture.json")) as fh:
    LISTING = json.load(fh)
HUMAN = [d for d in LISTING if d["organism"][0]["ontology_term_id"] == plan.HUMAN_TAXON]


def _by(pred):
    return next(d for d in LISTING if pred(d))


def _ids(items):
    return [i["dataset_version_id"] for i in items]


def test_filters_keep_only_human_primary_allowed_assays():
    p = plan.build_plan(LISTING, set())
    mouse = _by(lambda d: d["organism"][0]["label"] != "Homo sapiens")
    no_primary = _by(lambda d: d["primary_cell_count"] == 0)
    methyl = _by(lambda d: any(a["ontology_term_id"] == "EFO:0030027" for a in d["assay"]))
    kept = set(_ids(p["items"]))
    assert not kept & {mouse["dataset_version_id"], no_primary["dataset_version_id"],
                       methyl["dataset_version_id"]}
    assert p["totals"]["skipped_no_primary"] == sum(d["primary_cell_count"] == 0 for d in LISTING)
    for i in p["items"]:
        assert i["assays"] and all(a["ontology_term_id"] in agg.ASSAYS for a in i["assays"])


def test_spatial_only_dataset_is_in_scope_and_allow_list_is_the_workers():
    spatial = _by(lambda d: any(a["ontology_term_id"] == "EFO:0022857" for a in d["assay"]))
    item, _ = plan.classify(spatial)
    assert item and item["assays"][0]["ontology_term_id"] == "EFO:0022857"
    assert plan.allowed_assays.__globals__["aggregate"] is agg


def test_mixed_assays_record_only_the_allowed_ones():
    d = copy.deepcopy(HUMAN[0])
    d["assay"] = d["assay"] + [{"label": "snmC-Seq2", "ontology_term_id": "EFO:0030027"}]
    item, _ = plan.classify(d)
    assert "EFO:0030027" not in {a["ontology_term_id"] for a in item["assays"]}


def test_done_markers_are_subtracted_and_counted():
    everything = plan.build_plan(LISTING, set())["items"]
    done = {everything[0]["dataset_version_id"]}
    p = plan.build_plan(LISTING, done)
    assert _ids(p["items"]) == _ids(everything[1:])
    assert p["totals"]["done"] == 1 and p["totals"]["in_scope"] == len(everything)


def test_largest_first_with_machine_class_by_estimate():
    items = plan.build_plan(LISTING, set())["items"]
    assert [i["est_nnz"] for i in items] == sorted((i["est_nnz"] for i in items), reverse=True)
    big = items[0]
    assert big["cell_count"] >= 4_000_000 and big["machine_class"] == plan.HIGHMEM
    assert all(i["machine_class"] == plan.STANDARD for i in items[1:])
    assert big["est_nnz"] == round(big["cell_count"] * _by(
        lambda d: d["dataset_version_id"] == big["dataset_version_id"])["mean_genes_per_cell"])


@pytest.mark.parametrize("nnz,cells,want", [
    (4.9e9, 2.9e6, plan.STANDARD), (5e9, 1, plan.HIGHMEM), (1, 3_000_000, plan.HIGHMEM),
    (1.56e10, 4_140_453, plan.HIGHMEM), (0, 0, plan.STANDARD)])
def test_machine_class_thresholds(nnz, cells, want):
    assert plan.machine_class(nnz, cells) == want


def test_superseded_only_when_dataset_has_a_different_current_version():
    cur = HUMAN[0]
    out = plan.build_plan(LISTING, set(), {"old-version": cur["dataset_id"], "gone": "nope",
                                           cur["dataset_version_id"]: cur["dataset_id"]})["superseded"]
    assert out == [{"dataset_id": cur["dataset_id"], "dataset_version_id": "old-version",
                    "current_dataset_version_id": cur["dataset_version_id"]}]


def test_breaker_trips_on_count_or_size_unless_approved():
    n = len(plan.build_plan(LISTING, set())["items"])
    assert not plan.build_plan(LISTING, set(), max_datasets=n, max_gb=1e9)["breaker"]["tripped"]
    assert plan.build_plan(LISTING, set(), max_datasets=n - 1, max_gb=1e9)["breaker"]["tripped"]
    assert plan.build_plan(LISTING, set(), max_datasets=n, max_gb=1.0)["breaker"]["tripped"]


def _run(tmp_path, monkeypatch, extra, done=()):
    monkeypatch.setattr(plan, "list_markers", lambda bucket, **kw: {plan.MARKER: set(done), plan.REJECTED: set()})
    monkeypatch.setattr(plan, "read_marker_dataset_id", lambda v, b, client=None, root=None: None)
    out = tmp_path / "plan.json"
    code = plan.main(["--listing", os.path.join(HERE, "listing_fixture.json"), "--out", str(out)] + extra)
    return code, json.loads(out.read_text())


def test_main_exit_codes(tmp_path, monkeypatch):
    code, p = _run(tmp_path, monkeypatch, ["--max-datasets", "1"])
    assert code == plan.BREAKER_EXIT and p["breaker"]["tripped"] and p["items"]
    code, p = _run(tmp_path, monkeypatch, ["--max-datasets", "1", "--approve-large"])
    assert code == 0 and p["breaker"]["approved"]
    code, _ = _run(tmp_path, monkeypatch, [])
    assert code == 0


def test_main_reports_done_ids_missing_from_listing(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(plan, "read_marker_dataset_id", lambda v, b, client=None, root=None: seen.append(v) or "d")
    monkeypatch.setattr(plan, "_client", lambda: object())
    monkeypatch.setattr(plan, "list_markers", lambda bucket, **kw: {
        plan.MARKER: {"orphan", HUMAN[0]["dataset_version_id"]}, plan.REJECTED: set()})
    plan.main(["--listing", os.path.join(HERE, "listing_fixture.json"), "--out", str(tmp_path / "p.json")])
    assert seen == ["orphan"]


def test_main_exits_2_on_error(tmp_path, monkeypatch):
    def boom(bucket, **kw):
        raise RuntimeError("no gcs")
    monkeypatch.setattr(plan, "list_markers", boom)
    assert plan.main(["--listing", os.path.join(HERE, "listing_fixture.json"),
                      "--out", str(tmp_path / "p.json")]) == 2


def test_empty_listing_never_reports_superseded():
    assert plan.build_plan([], {"v1", "v2"}, {"v1": "d1", "v2": "d2"})["superseded"] == []
    truncated = LISTING[:1]
    other = LISTING[1]
    p = plan.build_plan(truncated, set(), {"old": other["dataset_id"]})
    assert p["superseded"] == []


def test_rejected_is_done_not_retried_and_outside_the_breaker():
    vids = [i["dataset_version_id"] for i in plan.build_plan(LISTING, set())["items"]]
    p = plan.build_plan(LISTING, set(), rejected={vids[0]}, max_datasets=len(vids) - 1, max_gb=1e9)
    assert vids[0] not in _ids(p["items"]) and not p["breaker"]["tripped"]
    assert [r["dataset_version_id"] for r in p["rejected"]] == [vids[0]]
    assert p["totals"]["rejected"] == 1 and p["totals"]["done"] == 0
    assert p["totals"]["todo"] == len(vids) - 1


def test_list_markers_and_root_reach_the_plan():
    names = ["r/agg/v1/a/_SUCCESS.json", "r/agg/v1/b/_REJECTED.json", "r/agg/v1/c/groups.parquet",
             "r/agg/v1/d/sub/_SUCCESS.json"]

    class C:
        def list_blobs(self, bucket, prefix):
            assert (bucket, prefix) == ("bk", "r/agg/v1/")
            return [type("B", (), {"name": n}) for n in names]
    m = plan.list_markers("bk", fmt=1, root="r", client=C())
    assert m == {plan.MARKER: {"a"}, plan.REJECTED: {"b"}}
    assert plan.build_plan(LISTING, set(), bucket="bk", root="r")["agg_prefix"] == \
        f"gs://bk/r/agg/v{plan.aggregate.FORMAT_VERSION}/"
