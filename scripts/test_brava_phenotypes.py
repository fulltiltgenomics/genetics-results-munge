"""Unit tests for the parts of `brava_phenotypes.py` that a wrong answer would not make obvious:
the sex tags in the file codes, and the deduplication of Table S5."""

from brava_phenotypes import build_items, parse_listing, per_ancestry, phenocode, phenotype_id


def test_all_is_a_sex_tag_and_f_is_part_of_the_name():
    assert phenotype_id("AFib_ALL") == "AFib"
    assert phenotype_id("BreastCanc_F") == "BreastCanc"
    assert phenocode("AFib_ALL", None) == "AFib"
    assert phenocode("BreastCanc_F", None) == "BreastCanc_F"
    assert phenocode("AFib_ALL", "non_EUR") == "AFib|non_EUR"
    assert phenocode("BreastCanc_F", "EUR") == "BreastCanc_F|EUR"


def test_parse_listing_reads_the_stratum_off_the_file_name():
    listing = [
        "gs://b/p/AFib_ALL_gene_meta_analysis_100_cutoff.tsv.gz",
        "AFib_ALL_gene_meta_analysis_100_cutoff.non_EUR.tsv.gz",
        "AFib_ALL_variant_meta_analysis_100_cutoff.EUR.vcf.gz",
        "",
    ]
    assert parse_listing(listing) == {"AFib_ALL": {None, "non_EUR"}}


def _tables():
    return {
        "Table S1": [{"Description": "Atrial fibrillation", "Sex": "Both", "Phenotype ID": "AFib"}],
        "Table S2": [{"Description": "Height", "Sex": "Both", "Phenotype ID": "Height"}],
        "Table S4": [
            {"Phenotype ID": "AFib", "Ancestry": a, "Biobank ID": "bb", "N cases": 10.0, "N controls": 90.0}
            for a in ("EUR", "MID")
        ],
        "Table S5": [
            {"Phenotype ID": "Height", "Ancestry": "EUR", "Biobank ID": "bb", "N": 500.0},
            {"Phenotype ID": "Height", "Ancestry": "EUR", "Biobank ID": "bb", "N": 500.0},
            {"Phenotype ID": "Height", "Ancestry": "AFR", "Biobank ID": "bb", "N": 40.0},
            {"Phenotype ID": "Height", "Ancestry": "AFR", "Biobank ID": "bb", "N": 30.0},
        ],
        "Table S6": [{"Description": "Atrial fibrillation", "Phenotype ID": "AFib", "N cases": 20.0, "N controls": 180.0}],
        "Table S7": [{"Description": "Height", "Phenotype ID": "Height", "N": 1070.0}],
    }


def test_repeated_s5_keys_count_once_and_the_smaller_n_wins():
    _, quantitative = per_ancestry(_tables())
    assert quantitative[("Height", "EUR")] == 500.0
    assert quantitative[("Height", "AFR")] == 30.0


def test_non_eur_sums_every_other_ancestry_including_mid():
    tables = _tables()
    binary, quantitative = per_ancestry(tables)
    items = {
        item["phenocode"]: item
        for item in build_items(
            tables, {"AFib_ALL": {None, "non_EUR"}, "Height_ALL": {None}}, binary, quantitative
        )
    }
    # MID has no stratum file, so non_EUR is the only place its 10/90 can reach the output
    assert items["AFib|non_EUR"] == {
        "phenocode": "AFib|non_EUR",
        "phenostring": "Atrial fibrillation (non_EUR)",
        "category": "Both",
        "num_cases": 10,
        "num_controls": 90,
    }
    assert items["AFib"]["num_cases"] == 20
    assert items["Height"] == {
        "phenocode": "Height",
        "phenostring": "Height",
        "category": "Both",
        "num_samples": 530,
    }
