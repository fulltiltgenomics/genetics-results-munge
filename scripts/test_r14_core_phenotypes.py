"""Unit tests for `r14_core_phenotypes.py`: the row conversion, the duplicate guard, and
`--previous` accepting only an additive change against the previously served JSON while
carrying its pheweb-computed keys over."""

import pytest

from r14_core_phenotypes import build, carry_over, compare

HEADER = "phenocode\tname\tcategory\tcategory_index\tnum_cases\tnum_controls\n"


def _tsv(*rows):
    return HEADER + "".join("\t".join(r) + "\n" for r in rows)


BINARY = ("AB1_CMV_WIDE", "Cytomegaloviral disease, including avohilmo", "I Infections (AB1_)", "1", "880", "484648")
QUANT = ("BMI_IRN", "BMI, inverse-rank normalized", "Quantitative endpoints", "47", "361936", "0")


def test_rows_use_pheweb_keys_and_keep_zero_controls():
    items = build(_tsv(BINARY, QUANT))
    assert items[0] == {"phenocode": "AB1_CMV_WIDE", "phenostring": BINARY[1], "category": BINARY[2],
                        "category_index": 1, "num_cases": 880, "num_controls": 484648}
    assert items[1]["num_controls"] == 0


def test_duplicate_phenocode_aborts():
    with pytest.raises(SystemExit, match="duplicate"):
        build(_tsv(BINARY, BINARY))


def test_check_passes_on_additive_change_and_tolerates_duplicated_previous_entry():
    items = build(_tsv(BINARY, QUANT))
    previous = [dict(items[1]), dict(items[1])]
    assert compare(items, previous) == []


def test_check_reports_changed_field_and_dropped_code():
    items = build(_tsv(BINARY))
    previous = [dict(items[0], num_cases=881), {"phenocode": "GONE", "phenostring": "x"}]
    problems = compare(items, previous)
    assert any(p.startswith("AB1_CMV_WIDE.num_cases") for p in problems)
    assert any(p.startswith("GONE:") for p in problems)
    assert len(problems) == 2


def test_carry_over_copies_pheweb_keys_only_where_pheweb_had_the_code():
    items = build(_tsv(BINARY, QUANT))
    previous = [dict(items[1], gc_lambda={"0.5": 0.9}, num_gw_significant=3, num_cases_prev=100,
                     assoc_files="/mnt/disks/x.gz")]
    assert carry_over(items, previous) == 1
    assert items[1]["gc_lambda"] == {"0.5": 0.9} and items[1]["num_gw_significant"] == 3
    assert items[1]["num_cases_prev"] == 100 and "assoc_files" not in items[1]
    assert "gc_lambda" not in items[0]
