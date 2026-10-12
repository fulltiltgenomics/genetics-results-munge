"""Unit tests for the mapping-TSV row conversion in `r14_meta_phenotypes.py`: which keys are
omitted for a quantitative trait, the fixed cohort and family order, and the sumstats header
overriding the TSV on cohort presence."""

import csv

from r14_meta_phenotypes import COHORTS, header_cohorts, parse_row, restrict_to_header

HEADER = (
    "phenocode name category category_index aou_phenocode ukbb_phenocode ukbb_phenocode_type "
    "mvp_phenocode num_cases num_controls".split()
    + [f"{c}_n_cases" for c in COHORTS]
    + [f"{c}_n_controls" for c in COHORTS]
)


def _row(phenocode, category, num_controls, counts):
    row = dict.fromkeys(HEADER, "NA")
    row.update(phenocode=phenocode, name=phenocode.lower(), category=category, category_index="1",
               num_cases=str(sum(c for c, _ in counts.values())), num_controls=num_controls)
    for cohort, (n_cases, n_controls) in counts.items():
        row[f"{cohort}_n_cases"] = str(n_cases)
        row[f"{cohort}_n_controls"] = "NA" if n_controls is None else str(n_controls)
    return row


def _parse(tmp_path, rows):
    # round-trip through a real TSV so the fixture exercises the same reader as the script
    path = tmp_path / "mapping.tsv"
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with path.open() as f:
        return [parse_row(r) for r in csv.DictReader(f, delimiter="\t")]


def test_three_row_fixture(tmp_path):
    quant, every, two = _parse(tmp_path, [
        _row("LAB_X", "Labs", "NA", {"fg": (1000, None), "ukbb_eur": (500, None)}),
        _row("ALL", "Diseases", "100", {c: (10, 20) for c in reversed(COHORTS)}),
        _row("AB1_X_WIDE", "Infections", "900", {"ukbb_eur": (5, 400), "fg": (7, 500)}),
    ])

    assert "num_controls" not in quant
    assert quant["num_cases"] == 1500
    assert quant["n_by_cohort"] == {"fg": {"n_cases": 1000}, "ukbb_eur": {"n_cases": 500}}

    assert every["contributing_cohorts"] == COHORTS
    assert every["cohort_families"] == ["fg", "aou", "mvp", "ukbb"]
    assert every["n_by_cohort"]["mvp_amr"] == {"n_cases": 10, "n_controls": 20}
    assert every["num_controls"] == 100

    assert two["contributing_cohorts"] == ["fg", "ukbb_eur"]
    assert two["cohort_families"] == ["fg", "ukbb"]
    assert two["category_index"] == 1
    assert "cohorts" not in two


def test_header_cohorts_ignores_meta_and_leave_one_out_columns_but_keeps_unknown_cohorts():
    header = ["chrom", "fg_beta", "fg_sebeta", "ukbb_eur_beta", "all_inv_var_meta_beta",
              "leave_cohort_fg_inv_var_meta_beta", "new_cohort_beta"]
    assert header_cohorts(header) == {"fg", "ukbb_eur", "new_cohort"}


def test_tsv_cohort_missing_from_header_is_dropped_with_warning_and_totals_recomputed(tmp_path):
    lab, binary = _parse(tmp_path, [
        _row("LAB_X", "Labs", "NA", {"fg": (1000, None), "mvp_eur": (300, None), "mvp_amr": (40, None)}),
        _row("BIN", "Diseases", "720", {"fg": (10, 500), "mvp_afr": (3, 220)}),
    ])

    warnings, errors = restrict_to_header(lab, {"fg", "mvp_eur"})
    assert errors == []
    assert len(warnings) == 1 and "LAB_X" in warnings[0] and "mvp_amr" in warnings[0] and "40" in warnings[0]
    assert lab["contributing_cohorts"] == ["fg", "mvp_eur"]
    assert lab["cohort_families"] == ["fg", "mvp"]
    assert lab["n_by_cohort"] == {"fg": {"n_cases": 1000}, "mvp_eur": {"n_cases": 300}}
    assert lab["num_cases"] == 1300
    assert "num_controls" not in lab

    warnings, errors = restrict_to_header(binary, {"fg"})
    assert (len(warnings), errors) == (1, [])
    assert binary["cohort_families"] == ["fg"]
    assert (binary["num_cases"], binary["num_controls"]) == (10, 500)


def test_header_cohort_without_tsv_n_is_an_error(tmp_path):
    (item,) = _parse(tmp_path, [_row("X", "Diseases", "500", {"fg": (10, 500)})])
    warnings, errors = restrict_to_header(item, {"fg", "ukbb_eur"})
    assert warnings == []
    assert len(errors) == 1 and "ukbb_eur" in errors[0]


def test_matching_header_leaves_tsv_totals_alone(tmp_path):
    (item,) = _parse(tmp_path, [_row("X", "Diseases", "500", {"fg": (10, 500)})])
    before = dict(item)
    assert restrict_to_header(item, {"fg"}) == ([], [])
    assert item == before
