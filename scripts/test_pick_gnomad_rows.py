"""pick_gnomad_rows: which gnomAD row a munge joins against, on both layouts of the sites file."""

import polars as pl

from sumstat_utils import pick_gnomad_rows

# v4.0 layout: a variant in both callsets has a genomes and an exomes row
TWO_ROW = pl.DataFrame({
    "chr": [1, 1, 1, 23],
    "pos": [100, 100, 200, 300],
    "ref": ["A", "A", "C", "G"],
    "alt": ["T", "T", "G", "A"],
    "AF": [0.1, 0.2, 0.3, 0.4],
    "genome_or_exome": ["e", "g", "e", "g"],
})

# v4.1.1 layout as build_gnomad_annotation.py writes it: one row per variant, either source
ONE_ROW = TWO_ROW.filter(~((pl.col("pos") == 100) & (pl.col("genome_or_exome") == "e")))


def keys(df: pl.DataFrame) -> list[tuple]:
    return sorted(df.select("chr", "pos", "ref", "alt", "genome_or_exome").iter_rows())


def test_any_keeps_every_variant_of_a_one_row_file():
    assert keys(pick_gnomad_rows(ONE_ROW)) == keys(ONE_ROW)


def test_any_prefers_genomes_where_both_rows_exist():
    out = pick_gnomad_rows(TWO_ROW)
    assert keys(out) == [(1, 100, "A", "T", "g"), (1, 200, "C", "G", "e"), (23, 300, "G", "A", "g")]
    assert out.filter(pl.col("pos") == 100)["AF"].to_list() == [0.2]


def test_source_restriction_drops_the_other_sources_rows():
    assert keys(pick_gnomad_rows(ONE_ROW, "g")) == [(1, 100, "A", "T", "g"), (23, 300, "G", "A", "g")]
    assert keys(pick_gnomad_rows(TWO_ROW, "e")) == [(1, 100, "A", "T", "e"), (1, 200, "C", "G", "e")]


def test_key_names_the_variant_columns():
    renamed = TWO_ROW.rename({"chr": "gnomad_chr", "pos": "gnomad_pos", "ref": "gnomad_ref", "alt": "gnomad_alt"})
    out = pick_gnomad_rows(renamed, key=("gnomad_chr", "gnomad_pos", "gnomad_ref", "gnomad_alt"))
    assert out.height == 3
