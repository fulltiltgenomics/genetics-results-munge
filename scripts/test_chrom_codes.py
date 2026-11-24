"""The suite's chromosome coding, in its polars and plain-Python forms."""

import polars as pl

from peak_utils import CANONICAL_CHROMS, numeric_chrom_expr
from sumstat_utils import AWK_CHROM_CODES, CHROM_CODES, CHROM_NAMES, chrom_code, chrom_code_expr

AUTOSOMES = [str(c) for c in range(1, 23)]
EXPECTED = {
    **{c: c for c in AUTOSOMES},
    "X": "23", "Y": "24", "M": "26", "MT": "26",
    "chr1": "1", "chrX": "23", "chrY": "24", "chrM": "26", "CHRx": "23", "chrmt": "26",
    # already numeric, 25 included: it passes through as a string and is nobody's code
    "23": "23", "24": "24", "25": "25", "26": "26",
    "GL000220.1": "GL000220.1", "chrUn_KI270742v1": "Un_KI270742v1", "": "",
}


def test_plain_python_form():
    assert {name: chrom_code(name) for name in EXPECTED} == EXPECTED


def test_polars_form_matches_and_keeps_nulls():
    names = [*EXPECTED, None]
    df = pl.DataFrame({"chr": names}, schema={"chr": pl.Utf8})
    for coded in (df.select(chrom_code_expr())["chr"], df.select(chrom_code_expr(pl.col("chr")))["chr"],
                  df.select(chrom_code_expr(df["chr"]).alias("chr"))["chr"]):
        assert coded.dtype == pl.Utf8
        assert coded.to_list() == [*EXPECTED.values(), None]


def test_integer_column_passes_through():
    assert pl.DataFrame({"chr": [1, 23, 26]}).select(chrom_code_expr())["chr"].to_list() == ["1", "23", "26"]


def test_the_cast_is_the_callers():
    df = pl.DataFrame({"chr": ["chrX", "MT", "GL000220.1"]})
    assert df.select(chrom_code_expr().cast(pl.Int32, strict=False))["chr"].to_list() == [23, 26, None]


def test_no_code_is_25_and_names_reverse_codes():
    assert "25" not in CHROM_CODES.values() and "25" not in CANONICAL_CHROMS
    assert CANONICAL_CHROMS == {*AUTOSOMES, "23", "24", "26"}
    assert CHROM_NAMES == {"23": "X", "24": "Y", "26": "M"}
    assert all(CHROM_CODES[name] == code for code, name in CHROM_NAMES.items())
    assert all(f'code["{name}"] = {code};' in AWK_CHROM_CODES for name, code in CHROM_CODES.items())


def test_peak_family_uses_the_same_function():
    df = pl.DataFrame({"chrom": [*EXPECTED, None]}, schema={"chrom": pl.Utf8})
    assert df.select(numeric_chrom_expr())["chrom"].to_list() == [*EXPECTED.values(), None]
