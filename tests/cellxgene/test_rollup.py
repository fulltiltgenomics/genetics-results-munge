"""rollup.py: the checked-in tables, and the rules on a synthetic ontology.

Run: python -m pytest tests/cellxgene/test_rollup.py
Needs cellxgene-ontology-guide (tissue rule, label lookups); no network.
"""
import csv
import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "cellxgene"))
import rollup

pytest.importorskip("cellxgene_ontology_guide")

TABLES = Path(rollup.__file__).resolve().parent / "rollup"


def _rows(name):
    with open(TABLES / name) as f:
        return list(csv.DictReader(f, delimiter="\t"))


@pytest.fixture(scope="module")
def cell_types():
    return {r["cell_type_ontology_term_id"]: r for r in _rows("cell_type_rollup.tsv")}


@pytest.fixture(scope="module")
def tissues():
    return {r["tissue_ontology_term_id"]: r for r in _rows("tissue_rollup.tsv")}


@pytest.mark.parametrize("native, target, method", [
    ("CL:0000624", "CL:0000084", "ancestor"),  # CD4-positive, alpha-beta T cell -> T cell
    ("CL:0002063", "CL:0002063", "self"),      # type II pneumocyte is itself a target
    ("CL:0000186", "CL:0000057", "pin"),       # myofibroblast -> fibroblast, only via the pin
    ("CL:0000632", "CL:0000632", "self"),      # hepatic stellate cell is_a fibroblast, but a target
    ("unknown", "unknown", "unknown"),
])
def test_known_cell_type_mappings(cell_types, native, target, method):
    r = cell_types[native]
    assert (r["rollup_term_id"], r["method"]) == (target, method)


@pytest.mark.parametrize("native, general, method", [
    ("UBERON:0008952", "UBERON:0002048", "tissue_general"),     # upper lobe of left lung -> lung
    ("UBERON:0000178", "UBERON:0000178", "tissue_general"),     # blood is itself on the list
    ("CL:0002586", "CL:0002586", "cell_culture_self"),          # primary cell culture tissue
    ("CVCL_0031", "CVCL_0031", "cell_culture_self"),            # cell line
])
def test_known_tissue_mappings(tissues, native, general, method):
    r = tissues[native]
    assert (r["tissue_general_ontology_term_id"], r["method"]) == (general, method)


def test_culture_tissues_are_never_marked_tissue_only(tissues):
    for term, r in tissues.items():
        if not term.startswith("UBERON:"):
            assert r["method"] == "cell_culture_self" and "tissue" not in r["tissue_types"].split(","), term


def test_every_native_term_maps_to_exactly_one_target(cell_types):
    rows = _rows("cell_type_rollup.tsv")
    assert len(rows) == len(cell_types), "a native term appears twice"
    targets = {r["term_id"] for r in rollup._read_tsv(TABLES / "cell_type_targets.tsv")}
    pins = {r["native_term_id"]: r["target_term_id"] for r in rollup._read_tsv(TABLES / "pins.tsv")}
    for term, r in cell_types.items():
        assert r["rollup_term_id"] and ";" not in r["rollup_term_id"]
        if r["method"] in ("self", "ancestor", "pin", "pin_descendant"):
            assert r["rollup_term_id"] in targets, term
        else:
            assert r["rollup_term_id"] == term, term
    for native, target in pins.items():
        if native in cell_types:
            assert cell_types[native]["rollup_term_id"] == target


def test_metadata_records_the_releases():
    meta = json.loads((TABLES / "metadata.json").read_text())
    assert meta["ontologies"]["CL"] and meta["ontologies"]["UBERON"]
    assert meta["cl_obo"]["data_version"].endswith(meta["ontologies"]["CL"].lstrip("v") + "/cl-basic.owl")
    assert meta["n_pins"] == len(rollup._read_tsv(TABLES / "pins.tsv"))


# synthetic CL: cell <- A <- A1 <- A1x ; cell <- B ; A1 <- P (pinned to B) <- P1 ; P <- Q, Q is_a A2 too
#               C is_a A and B (tie)
OBO = """format-version: 1.2
data-version: cl/releases/{v}/cl-basic.owl

[Term]
id: CL:0000000
name: cell

[Term]
id: CL:1000001
name: a
is_a: CL:0000000

[Term]
id: CL:1000002
name: a1
is_a: CL:1000001

[Term]
id: CL:1000003
name: a1x
is_a: CL:1000002

[Term]
id: CL:1000004
name: b
is_a: CL:0000000

[Term]
id: CL:1000005
name: p
is_a: CL:1000002

[Term]
id: CL:1000006
name: p1
is_a: CL:1000005

[Term]
id: CL:1000007
name: c
is_a: CL:1000001
is_a: CL:1000004

[Term]
id: CL:1000008
name: a2
is_a: CL:1000001

[Term]
id: CL:1000009
name: q
is_a: CL:1000005
is_a: CL:1000008

[Term]
id: CL:1000010
name: orphan

[Term]
id: CL:1000011
name: old
is_obsolete: true
"""


def _synthetic(tmp_path):
    from cellxgene_ontology_guide.ontology_parser import OntologyParser
    v = OntologyParser("7.1.0").cxg_schema.supported_ontologies["CL"]["version"].lstrip("v")
    obo = tmp_path / "cl.obo"
    obo.write_text(OBO.format(v=v))
    targets = tmp_path / "targets.tsv"
    targets.write_text("# comment\nterm_id\tlabel\tlineage\torigin\n"
                       "CL:1000002\ta1\tx\ttest\nCL:1000008\ta2\tx\ttest\n"
                       "CL:1000004\tb\tx\ttest\nCL:1000001\ta\tx\ttest\n")
    pins = tmp_path / "pins.tsv"
    pins.write_text("native_term_id\tnative_label\ttarget_term_id\ttarget_label\treason\n"
                    "CL:1000005\tp\tCL:1000004\tb\tsynthetic\n")
    terms = ["CL:1000003", "CL:1000002", "CL:1000005", "CL:1000006", "CL:1000007", "CL:1000009",
             "CL:1000010", "CL:1000011", "CL:9999999", "unknown"]
    listing = tmp_path / "datasets.json"
    listing.write_text(json.dumps([
        {"organism": [{"ontology_term_id": "NCBITaxon:9606"}], "primary_cell_count": 10,
         "schema_version": "7.1.0", "assay": [{"label": "10x 3' v3"}],
         "cell_type": [{"ontology_term_id": t, "label": t} for t in terms],
         "tissue": [{"ontology_term_id": "UBERON:0008952", "label": "upper lobe of left lung",
                     "tissue_type": "tissue"},
                    {"ontology_term_id": "CL:0000084", "label": "T cell", "tissue_type": "primary cell culture"}]},
        {"organism": [{"ontology_term_id": "NCBITaxon:10090"}], "primary_cell_count": 10,
         "schema_version": "7.1.0", "assay": [], "cell_type": [{"ontology_term_id": "CL:1000004", "label": "b"}],
         "tissue": []},
    ]))
    return listing, obo, targets, pins


def test_rules_on_synthetic_ontology(tmp_path):
    listing, obo, targets, pins = _synthetic(tmp_path)
    _, ct_rows, ts_rows = rollup.build(listing, obo, tmp_path / "out", targets, pins)
    got = {r[0]: (r[2], r[4], r[5]) for r in ct_rows}
    assert got == {
        "CL:1000003": ("CL:1000002", "ancestor", ""),          # most specific target, not a
        "CL:1000002": ("CL:1000002", "self", ""),
        "CL:1000005": ("CL:1000004", "pin", ""),               # pin beats the is_a target a1
        "CL:1000006": ("CL:1000004", "pin_descendant", ""),
        "CL:1000007": ("CL:1000004", "ancestor", "CL:1000001"),  # tie: b listed before a
        "CL:1000009": ("CL:1000004", "pin_descendant", ""),    # a2 is not below the pinned term
        "CL:1000010": ("CL:1000010", "no_target", ""),
        "CL:1000011": ("CL:1000011", "not_in_ontology", ""),
        "CL:9999999": ("CL:9999999", "not_in_ontology", ""),
        "unknown": ("unknown", "unknown", ""),
    }
    assert "CL:1000004" not in got  # listed only by the mouse dataset
    assert {r[0]: (r[3], r[5]) for r in ts_rows} == {
        "UBERON:0008952": ("UBERON:0002048", "tissue_general"),
        "CL:0000084": ("CL:0000084", "cell_culture_self"),
    }


def test_rerun_is_byte_identical(tmp_path):
    listing, obo, targets, pins = _synthetic(tmp_path)
    outs = [tmp_path / "o1", tmp_path / "o2"]
    for o in outs:
        rollup.build(listing, obo, o, targets, pins)
    for name in ("cell_type_rollup.tsv", "tissue_rollup.tsv", "metadata.json"):
        assert (outs[0] / name).read_bytes() == (outs[1] / name).read_bytes(), name


@pytest.mark.parametrize("pin_row, msg", [
    ("CL:1000005\tp\tCL:1000004\tb\t \n", "empty reason"),
    ("CL:1000005\tp\tCL:1000010\torphan\tr\n", "not in the target list"),
    ("CL:1000005\twrong\tCL:1000004\tb\tr\n", "native label"),
])
def test_bad_pin_is_refused(tmp_path, pin_row, msg):
    listing, obo, targets, pins = _synthetic(tmp_path)
    pins.write_text(pins.read_text().splitlines()[0] + "\n" + pin_row)
    with pytest.raises(SystemExit, match=msg):
        rollup.build(listing, obo, tmp_path / "out", targets, pins)


def test_target_label_typo_is_refused(tmp_path):
    listing, obo, targets, pins = _synthetic(tmp_path)
    targets.write_text(targets.read_text().replace("\ta2\t", "\ta 2\t"))
    with pytest.raises(SystemExit, match="CL:1000008"):
        rollup.build(listing, obo, tmp_path / "out", targets, pins)


def _ct(parents, targets, pins=None):
    labels = {t: t for t in parents}
    return rollup.CellTypeRollup(parents, labels, set(), targets, pins or {})


def test_most_specific_beats_file_order():
    # root <- anc <- mid <- leaf; the ancestor target is listed before the descendant target
    parents = {"root": [], "anc": ["root"], "mid": ["anc"], "leaf": ["mid"]}
    ct = _ct(parents, ["anc", "mid"])
    assert ct.map("leaf") == ("mid", "ancestor", [])


def test_pin_descendant_keeps_its_own_target_below_the_pin():
    # root <- a <- p (pinned to z) <- t (a target) <- leaf ; z is a separate branch
    parents = {"root": [], "a": ["root"], "z": ["root"], "p": ["a"], "t": ["p"], "leaf": ["t"],
               "other": ["p"]}
    ct = _ct(parents, ["a", "z", "t"], {"p": "z"})
    assert ct.map("leaf") == ("t", "ancestor", [])
    assert ct.map("t") == ("t", "self", [])
    assert ct.map("other") == ("z", "pin_descendant", [])
    assert ct.map("p") == ("z", "pin", [])
