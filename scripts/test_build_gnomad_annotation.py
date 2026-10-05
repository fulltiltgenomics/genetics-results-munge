"""End-to-end tests of `build_gnomad_annotation.py` on small bgzipped inputs: which row of a
genome/exome pair survives, how contigs are coded and ordered, and which malformed inputs
stop the run instead of producing a file."""

import gzip
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).with_name("build_gnomad_annotation.py")
COLUMNS = ["#chr", "pos", "ref", "alt", "rsids", "filters", "AN", "AF", "most_severe",
           "gene_most_severe", "consequences", "genome_or_exome"]


def row(chrom, pos, ref, alt, an, source, most_severe="intron_variant", gene="GENE1"):
    return [str(chrom), str(pos), ref, alt, "NA", "", str(an), "1.0e-03", most_severe, gene,
            '[{"gene_symbol":"%s"}]' % gene, source]


def write_bgz(path: Path, rows, columns=COLUMNS) -> str:
    text = "".join("\t".join(r) + "\n" for r in [columns, *rows])
    with open(path, "wb") as fh:
        subprocess.run(["bgzip", "-c"], input=text.encode(), stdout=fh, check=True)
    return str(path)


def read(path: Path) -> list[list[str]]:
    with gzip.open(path, "rt") as f:
        return [line.rstrip("\n").split("\t") for line in f]


def build(tmp_path: Path, inputs: list[str], *extra: str):
    sites, cons = tmp_path / "sites.tsv.bgz", tmp_path / "consequence.tsv.bgz"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), *inputs, "--output", str(sites), "--consequence-output", str(cons),
         "--min-free-gb", "0", "--threads", "1", *extra],
        capture_output=True, text=True,
    )
    return result, sites, cons


def build_ok(tmp_path: Path, rows, *extra: str):
    result, sites, cons = build(tmp_path, ["--gnomad", write_bgz(tmp_path / "in.tsv.bgz", rows)], *extra)
    assert result.returncode == 0, result.stderr
    return read(sites), read(cons)


def build_fails(tmp_path: Path, rows) -> str:
    for stale in ("sites.tsv.bgz", "consequence.tsv.bgz"):
        (tmp_path / stale).write_bytes(b"an earlier result")
    result, sites, cons = build(tmp_path, ["--gnomad", write_bgz(tmp_path / "in.tsv.bgz", rows)])
    assert result.returncode != 0
    assert not sites.exists() and not cons.exists()
    return result.stderr


def test_tie_keeps_the_genome_row(tmp_path):
    sites, _ = build_ok(tmp_path, [row(1, 100, "A", "C", 50, "e"), row(1, 100, "A", "C", 50, "g")])
    assert sites[1:] == [row(1, 100, "A", "C", 50, "g")]


def test_exome_row_kept_only_when_its_an_is_greater(tmp_path):
    rows = [
        row(1, 100, "A", "C", 10, "g", "intron_variant", "GENOME"),
        row(1, 100, "A", "C", 11, "e", "missense_variant", "EXOME"),
        row(1, 200, "A", "C", 11, "g", "intron_variant", "GENOME"),
        row(1, 200, "A", "C", 10, "e", "missense_variant", "EXOME"),
    ]
    sites, cons = build_ok(tmp_path, rows)
    assert sites == [COLUMNS, rows[1], rows[2]]
    assert cons == [
        ["#chr", "pos", "ref", "alt", "most_severe", "gene_most_severe"],
        ["1", "100", "A", "C", "missense_variant", "EXOME"],
        ["1", "200", "A", "C", "intron_variant", "GENOME"],
    ]


def test_an_is_compared_as_a_number(tmp_path):
    sites, _ = build_ok(tmp_path, [row(1, 100, "A", "C", 9, "g"), row(1, 100, "A", "C", 10, "e")])
    assert sites[1][-1] == "e"


def test_variant_in_one_source_is_kept(tmp_path):
    rows = [row(1, 100, "A", "C", 10, "e"), row(1, 200, "A", "C", 10, "g")]
    sites, cons = build_ok(tmp_path, rows)
    assert sites[1:] == rows
    assert [r[:4] for r in cons[1:]] == [r[:4] for r in rows]


def test_sex_chromosomes_are_numeric_and_x_is_23(tmp_path):
    rows = [row(22, 5, "A", "C", 10, "g"), row("X", 100, "A", "C", 10, "g"), row("chrx", 200, "A", "C", 10, "e"),
            row("Y", 100, "A", "C", 10, "g")]
    result, _, _ = build(tmp_path, ["--gnomad", write_bgz(tmp_path / "in.tsv.bgz", rows)])
    assert result.returncode != 0 and "contig 23 resumes" in result.stderr

    sites, cons = build_ok(tmp_path, [rows[0], rows[1], rows[3]])
    assert [r[0] for r in sites[1:]] == ["22", "23", "24"]
    assert [r[0] for r in cons[1:]] == ["22", "23", "24"]
    assert sites[2][1:] == rows[1][1:]


def test_output_is_in_numeric_contig_order_whatever_the_input_order(tmp_path):
    rows = [row("X", 7, "A", "C", 10, "g"), row(1, 100, "A", "C", 10, "g"), row(10, 5, "A", "C", 10, "g"),
            row(2, 50, "A", "C", 10, "g")]
    sites, cons = build_ok(tmp_path, rows)
    assert [r[0] for r in sites[1:]] == ["1", "2", "10", "23"]
    assert [r[0] for r in cons[1:]] == ["1", "2", "10", "23"]
    for name in ("sites.tsv.bgz", "consequence.tsv.bgz"):
        hit = subprocess.run(["tabix", str(tmp_path / name), "23:7-7"], capture_output=True, text=True, check=True)
        assert hit.stdout.split("\t")[:4] == ["23", "7", "A", "C"]
    assert not list(tmp_path.glob("*.part.*"))


def test_mitochondrion_code_is_a_choice_and_scaffolds_are_dropped(tmp_path):
    rows = [row(1, 1, "A", "C", 10, "g"), row("chr1_KI270706v1_random", 5, "A", "C", 10, "g"),
            row("chrM", 9, "A", "C", 10, "g")]
    assert [r[0] for r in build_ok(tmp_path, rows)[0][1:]] == ["1", "25"]
    assert [r[0] for r in build_ok(tmp_path, rows, "--mito-code", "26")[0][1:]] == ["1", "26"]


def test_same_source_duplicate_fails(tmp_path):
    err = build_fails(tmp_path, [row(1, 100, "A", "C", 10, "g"), row(1, 100, "A", "C", 12, "g")])
    assert "same source" in err and "1:100:A:C" in err


def test_three_rows_for_one_variant_fail(tmp_path):
    err = build_fails(tmp_path, [row(3, 100, "A", "C", 10, s) for s in "geg"])
    assert "3 rows for one variant" in err and "3:100:A:C" in err


def test_unsorted_input_fails(tmp_path):
    err = build_fails(tmp_path, [row(1, 200, "A", "C", 10, "g"), row(1, 100, "A", "C", 10, "g")])
    assert "backwards" in err and "1:100:A:C" in err


def test_contig_that_comes_back_fails(tmp_path):
    err = build_fails(tmp_path, [row(1, 1, "A", "C", 10, "g"), row(2, 1, "A", "C", 10, "g"), row(1, 2, "A", "C", 10, "g")])
    assert "resumes" in err and "1:2:A:C" in err


def test_unknown_source_fails(tmp_path):
    assert "neither g nor e" in build_fails(tmp_path, [row(1, 1, "A", "C", 10, "x")])


def test_multi_allelic_position_resolves_each_alt_separately(tmp_path):
    rows = [
        row(1, 100, "A", "T", 10, "g"),
        row(1, 100, "A", "C", 10, "g"),
        row(1, 100, "A", "T", 20, "e"),
        row(1, 100, "AG", "A", 30, "e"),
        row(1, 100, "A", "C", 5, "e"),
    ]
    sites, _ = build_ok(tmp_path, rows)
    assert sites[1:] == [rows[1], rows[2], rows[3]]


MERGED = [
    row(1, 100, "A", "C", 10, "e"),
    row(1, 100, "A", "C", 10, "g"),
    row(1, 100, "A", "G", 10, "g"),
    row(1, 150, "A", "C", 10, "g"),
    row(1, 150, "A", "G", 99, "e"),
    row(1, 150, "A", "G", 10, "g"),
    row(2, 5, "T", "C", 10, "e"),
    row("X", 7, "T", "C", 10, "e"),
    row("X", 7, "T", "C", 9, "g"),
    row("X", 8, "T", "C", 9, "g"),
    row("Y", 3, "T", "C", 9, "g"),
]


@pytest.mark.parametrize("with_source_column", [True, False])
def test_two_inputs_give_the_same_output_as_the_merged_input(tmp_path, with_source_column):
    expected = build_ok(tmp_path, MERGED)
    width = len(COLUMNS) if with_source_column else len(COLUMNS) - 1
    split = []
    for flag, source in (("--genomes", "g"), ("--exomes", "e")):
        rows = [r[:width] for r in MERGED if r[-1] == source]
        split += [flag, write_bgz(tmp_path / f"{source}.tsv.bgz", rows, COLUMNS[:width])]
    two = tmp_path / "two"
    two.mkdir()
    result, sites, cons = build(two, split)
    assert result.returncode == 0, result.stderr
    assert (read(sites), read(cons)) == expected
    assert len(expected[0]) - 1 == len({tuple(r[:4]) for r in MERGED})
