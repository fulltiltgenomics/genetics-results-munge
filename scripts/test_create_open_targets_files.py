"""Unit tests for the Open Targets `trait` value: a consumer cuts the `_(<studyId>)` suffix and
compares the rest by exact string equality, so the spelling is a contract."""

import polars as pl
import pytest

import create_open_targets_files as ot
from create_open_targets_files import sanitize_trait_name, trait_label


def test_spaces_become_underscores():
    assert trait_label("GCST004602", "Type 2 diabetes") == "Type_2_diabetes_(GCST004602)"


@pytest.mark.parametrize(
    "name, expected",
    [
        ("  Alcohol  withdrawal state ", "Alcohol_withdrawal_state"),
        ("herpes virus 7 IgG ", "herpes_virus_7_IgG"),
        ("line separator　wide", "line_separator_wide"),
        ("tab\tand\nnew\r\nline", "tabandnewline"),
        ("zero​width\x00null\x1f", "zerowidthnull"),
    ],
)
def test_no_whitespace_or_control_character_survives(name, expected):
    assert sanitize_trait_name(name) == expected


@pytest.mark.parametrize(
    "name",
    [
        "Disorders_of_adrenal_gland,_other_and/or_unspecified",
        "Alzheimer's_disease_(late_onset)",
        "Phosphatidylcholine_33:1_[M+H]1+/Phosphatidate_38:2_[M+NH4]1+_levels",
        "Airway_obstruction_(FEV1/FVC<70%)",
        "MetaCyc_pathway_(BIOTIN-BIOSYNTHESIS-PWY|biotin_biosynthesis_I)",
        "Sjögren's_syndrome",
        "IFNβ-1b_≥_70_years",
    ],
)
def test_punctuation_and_non_ascii_letters_are_kept(name):
    assert sanitize_trait_name(name) == name


def test_double_quote_becomes_single_quote():
    assert sanitize_trait_name('responded "yes" to "Are you a worrier?"') == (
        "responded_'yes'_to_'Are_you_a_worrier?'"
    )


def test_composed_and_decomposed_spellings_compare_equal():
    assert sanitize_trait_name("Ménière") == sanitize_trait_name("Ménière")


@pytest.mark.parametrize("name", [None, "", " \t\n", " ​"])
def test_study_without_a_usable_name_falls_back_to_the_accession(name):
    assert trait_label("GCST1", name) == "GCST1"


def test_same_name_gives_one_trait_per_study(tmp_path):
    (tmp_path / "study_metadata").mkdir()
    pl.DataFrame(
        {
            "studyId": ["GCST1", "GCST2", "GCST3", "GCST9"],
            "traitFromSource": ["Height", "Height", None, "not in the credible sets"],
        }
    ).write_parquet(tmp_path / "study_metadata" / "part-0.parquet")

    labels = ot.read_trait_labels(str(tmp_path), pl.Series(["GCST1", "GCST2", "GCST3", "GCST4"]))

    assert dict(labels.iter_rows()) == {
        "GCST1": "Height_(GCST1)",
        "GCST2": "Height_(GCST2)",
        "GCST3": "GCST3",
        "GCST4": "GCST4",
    }
