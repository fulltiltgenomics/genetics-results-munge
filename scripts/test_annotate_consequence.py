"""End-to-end tests of `annotate_consequence.py` on small files: that the two annotation
columns follow the consequence file, that nothing else in a row moves, that the two modes and
the two ways of reading the consequence file agree, and which inputs stop the run instead of
producing a file."""

import gzip
import json
import os
import random
import struct
import subprocess
import sys
from pathlib import Path

import pytest

import annotate_consequence

SCRIPT = Path(annotate_consequence.__file__)
COMBINED = ["#dataset", "data_type", "trait", "trait_original", "cell_type", "chr", "pos", "ref", "alt",
            "mlog10p", "beta", "se", "pip", "cs_id", "cs_size", "cs_min_r2", "aaf", "most_severe",
            "gene_most_severe"]
QTL = COMBINED + ["trait_chr", "trait_start", "trait_end"]
CONS_COLUMNS = ["#chr", "pos", "ref", "alt", "most_severe", "gene_most_severe"]
CONS = [
    ["1", "99", "A", "C", "upstream_gene_variant", "NEIGHBOUR"],
    ["1", "100", "A", "C", "missense_variant", "GENE1"],
    ["1", "100", "A", "G", "synonymous_variant", "GENE1"],
    ["1", "100", "AT", "A", "frameshift_variant", "GENE1"],
    ["1", "200", "T", "C", "intron_variant", "GENE2"],
    ["1", "250", "T", "C", "intron_variant", "UNASKED"],
    ["1", "90000", "G", "A", "intergenic_variant", "NA"],
    ["2", "100", "A", "C", "stop_gained", "WRONGCHROM"],
    ["23", "500", "C", "T", "stop_gained", "GENEX"],
    ["24", "500", "C", "T", "stop_gained", "GENEY"],
]


def served(chrom, pos, ref, alt, most_severe="NA", gene="NA", trait="T1", gene_pos=None):
    """One served row whose other columns hold the values a careless rewrite would mangle."""
    fields = ["DS", "GWAS", trait, 'tr"ait \'x\'', "", str(chrom), str(pos), ref, alt, "11.9101",
              "7.257e-02", "1.000e-02", "0.215", "chr1_1_A_G", "7", "NA", "6.493e-01", most_severe, gene]
    return fields if gene_pos is None else fields + ["1", str(gene_pos), ""]


def variant_sorted(gene_pos=None):
    kw = {"gene_pos": gene_pos}
    return [
        served(1, 100, "A", "C", **kw),
        served(1, 100, "A", "G", "intron_variant", "OLD", **kw),
        served(1, 100, "A", "T", "missense_variant", "OLD", **kw),
        served(1, 100, "A", "C", "missense_variant", "GENE1", trait="T2", **kw),
        served(1, 200, "T", "C", "", "", **kw),
        served(1, 300, "G", "A", "intron_variant", "GONE", **kw),
        served(1, 90000, "G", "A", **kw),
        served(3, 100, "A", "C", "stop_gained", "NOCHROM", **kw),
        served(23, 500, "C", "T", **kw),
        served(23, 500, "C", "G", **kw),
    ]


def to_bytes(columns, rows) -> bytes:
    return "".join("\t".join(r) + "\n" for r in [columns, *rows]).encode()


def write_plain(path: Path, data: bytes) -> Path:
    path.write_bytes(data)
    return path


def write_bgz(path: Path, data: bytes, *tabix_args: str) -> Path:
    with open(path, "wb") as fh:
        subprocess.run(["bgzip", "-c"], input=data, stdout=fh, check=True)
    if tabix_args:
        subprocess.run(["tabix", "-f", *tabix_args, str(path)], check=True)
    return path


def write_cons(tmp_path: Path, rows=CONS, columns=CONS_COLUMNS) -> Path:
    return write_bgz(tmp_path / "cons.tsv.bgz", to_bytes(columns, rows), "-s1", "-b2", "-e2")


def read_bytes(path: Path) -> bytes:
    data = path.read_bytes()
    return gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data


def expected(columns, rows, cons=CONS) -> bytes:
    table = {tuple(c[:4]): c[4:] for c in cons}
    ms = columns.index("most_severe")
    stamped = []
    for r in rows:
        r = list(r)
        r[ms:ms + 2] = table.get((r[5], r[6], r[7], r[8]), ["NA", "NA"])
        stamped.append(r)
    return to_bytes(columns, stamped)


def run(*args) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], capture_output=True, text=True)


def stamp(src: Path, *args):
    out = src.parent / "out" / src.name
    out.parent.mkdir(exist_ok=True)
    return run("--input", src, "--output", out, *args), out


def stamp_ok(src: Path, *args):
    result, out = stamp(src, *args)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout.splitlines()[-1]), out


def stamp_fails(src: Path, *args) -> str:
    result, out = stamp(src, *args)
    assert result.returncode != 0
    assert list(out.parent.iterdir()) == []
    return result.stderr


def index_settings(path: Path) -> tuple:
    with gzip.open(str(path) + ".tbi", "rb") as f:
        return struct.unpack_from("<6i", f.read(32), 8)


@pytest.mark.parametrize("access", ["scan", "regions"])
@pytest.mark.parametrize("container", ["plain", "bgzip"])
def test_merge_stamps_only_the_annotation_columns(tmp_path, access, container):
    rows = variant_sorted()
    data = to_bytes(COMBINED, rows)
    if container == "plain":
        src = write_plain(tmp_path / "in.tsv", data)
    else:
        src = write_bgz(tmp_path / "in.tsv.gz", data, "-s6", "-b7", "-e7")
    summary, out = stamp_ok(src, "--consequence", write_cons(tmp_path), "--consequence-access", access)
    assert (out.read_bytes()[:2] == b"\x1f\x8b") == (container == "bgzip")
    assert read_bytes(out) == expected(COMBINED, rows)
    assert summary["rows"] == 10 and summary["matched"] == 6
    assert summary["most_severe"] == {"na_before": 5, "na_after": 4, "changed_non_na": 4}
    assert summary["gene_most_severe"] == {"na_before": 5, "na_after": 5, "changed_non_na": 4}


def test_merge_cases_row_by_row(tmp_path):
    src = write_plain(tmp_path / "in.tsv", to_bytes(COMBINED, variant_sorted()))
    _, out = stamp_ok(src, "--consequence", write_cons(tmp_path))
    got = [line.split("\t")[-2:] for line in out.read_text().splitlines()[1:]]
    assert got == [
        ["missense_variant", "GENE1"],      # was NA
        ["synonymous_variant", "GENE1"],    # second alt of the position, old value replaced
        ["NA", "NA"],                       # third alt, not in the consequence file: blanked
        ["missense_variant", "GENE1"],      # same variant under another trait
        ["intron_variant", "GENE2"],        # was empty
        ["NA", "NA"],                       # position absent
        ["intergenic_variant", "NA"],
        ["NA", "NA"],                       # chromosome absent
        ["stop_gained", "GENEX"],           # chromosome 23
        ["NA", "NA"],
    ]


def test_bgzip_output_is_indexed_with_the_inputs_settings(tmp_path):
    rows = sorted(variant_sorted(gene_pos=1000), key=lambda r: int(r[6]) % 7)
    for i, r in enumerate(rows):
        r[20] = str(1000 + i)
    src = write_bgz(tmp_path / "in.qtl.tsv.gz", to_bytes(QTL, rows), "-s20", "-b21", "-e21")
    _, out = stamp_ok(src, "--consequence", write_cons(tmp_path), "--mode", "lookup")
    assert index_settings(out) == index_settings(src) == (0, 20, 21, 21, ord("#"), 0)
    hit = subprocess.run(["tabix", str(out), "1:1003-1003"], capture_output=True, check=True).stdout
    assert hit == expected(QTL, [rows[3]]).split(b"\n", 1)[1]


def test_csi_index_is_rebuilt_as_csi(tmp_path):
    src = write_bgz(tmp_path / "in.tsv.gz", to_bytes(COMBINED, variant_sorted()), "-C", "-s6", "-b7", "-e7")
    _, out = stamp_ok(src, "--consequence", write_cons(tmp_path))
    assert Path(str(out) + ".csi").exists() and not Path(str(out) + ".tbi").exists()


@pytest.mark.parametrize("access", ["scan", "regions"])
@pytest.mark.parametrize("container", ["plain", "bgzip"])
def test_lookup_writes_what_merge_writes(tmp_path, access, container):
    data = to_bytes(QTL, variant_sorted(gene_pos=1000))
    if container == "plain":
        src = write_plain(tmp_path / "in.tsv", data)
    else:
        src = write_bgz(tmp_path / "in.tsv.gz", data, "-s20", "-b21", "-e21")
    cons = write_cons(tmp_path)
    _, merged = stamp_ok(src, "--consequence", cons, "--consequence-access", access)
    merged_bytes, merged_summary = merged.read_bytes(), _
    summary, looked_up = stamp_ok(src, "--consequence", cons, "--consequence-access", access, "--mode", "lookup")
    assert looked_up.read_bytes() == merged_bytes
    assert read_bytes(looked_up) == expected(QTL, variant_sorted(gene_pos=1000))
    for key in ("rows", "matched", "most_severe", "gene_most_severe"):
        assert summary[key] == merged_summary[key]


@pytest.mark.parametrize("access", ["scan", "regions"])
def test_lookup_keeps_the_row_order_of_an_unsorted_file(tmp_path, access):
    rows = variant_sorted(gene_pos=1000)
    random.Random(0).shuffle(rows)
    src = write_plain(tmp_path / "in.qtl.tsv", to_bytes(QTL, rows))
    cons = write_cons(tmp_path)
    assert "not sorted at" in stamp_fails(src, "--consequence", cons)
    _, out = stamp_ok(src, "--consequence", cons, "--mode", "lookup", "--consequence-access", access)
    assert out.read_bytes() == expected(QTL, rows)


def test_regions_and_scan_agree_across_region_boundaries(tmp_path):
    rng = random.Random(1)
    gap = annotate_consequence.REGION_GAP
    positions = sorted({rng.choice([1, 2, gap - 1, gap, gap + 1, 3 * gap]) * k + rng.randrange(3)
                        for k in range(1, 60)})
    cons = [[str(c), str(p + d), "A", alt, "intron_variant", f"G{p}"]
            for c in (1, 2, 23) for p in positions for d in (-1, 0, 1) for alt in ("C", "G") if rng.random() < 0.6]
    cons = sorted({tuple(r[:4]): r for r in cons}.values(), key=lambda r: (int(r[0]), int(r[1])))
    rows = [served(c, p, "A", alt) for c in (1, 2, 23) for p in positions for alt in ("C", "T")]
    src = write_plain(tmp_path / "in.tsv", to_bytes(COMBINED, rows))
    cons_path = write_cons(tmp_path, cons)
    outputs = []
    for access in ("scan", "regions"):
        summary, out = stamp_ok(src, "--consequence", cons_path, "--consequence-access", access)
        outputs.append(out.read_bytes())
        assert 0 < summary["matched"] < len(rows)
    assert outputs[0] == outputs[1] == expected(COMBINED, rows, cons)


def test_small_inputs_take_the_regions_path_and_large_ones_scan(tmp_path):
    src = write_plain(tmp_path / "in.tsv", to_bytes(COMBINED, variant_sorted()))
    cons = write_cons(tmp_path)
    result, _ = stamp(src, "--consequence", cons)
    assert "read by regions" in result.stderr
    result, _ = stamp(src, "--consequence", cons, "--regions-max-bytes", "10")
    assert "read by scan" in result.stderr


@pytest.mark.parametrize("access", ["scan", "regions"])
@pytest.mark.parametrize("mode", ["merge", "lookup"])
def test_a_file_with_only_a_header_comes_out_unchanged(tmp_path, access, mode):
    src = write_plain(tmp_path / "in.tsv", to_bytes(COMBINED, []))
    summary, out = stamp_ok(src, "--consequence", write_cons(tmp_path), "--consequence-access", access,
                            "--mode", mode)
    assert out.read_bytes() == src.read_bytes() and summary["rows"] == 0


def test_clear_writes_na_without_a_consequence_file(tmp_path):
    rows = variant_sorted(gene_pos=1000)
    random.Random(0).shuffle(rows)
    src = write_bgz(tmp_path / "in.qtl.tsv.gz", to_bytes(QTL, rows), "-s20", "-b21", "-e21")
    summary, out = stamp_ok(src, "--clear")
    assert read_bytes(out) == expected(QTL, rows, cons=[])
    assert summary["matched"] == 0 and summary["most_severe"]["na_after"] == len(rows)
    assert index_settings(out) == index_settings(src)


def test_a_last_line_without_a_newline_stays_that_way(tmp_path):
    rows = variant_sorted()
    for columns, rws in ((COMBINED, rows), (QTL, variant_sorted(gene_pos=1000))):
        src = write_plain(tmp_path / "in.tsv", to_bytes(columns, rws)[:-1])
        _, out = stamp_ok(src, "--consequence", write_cons(tmp_path))
        assert out.read_bytes() == expected(columns, rws)[:-1]


@pytest.mark.parametrize("column", ["chr", "pos", "ref", "alt", "most_severe", "gene_most_severe"])
@pytest.mark.parametrize("mode", ["merge", "lookup", "clear"])
def test_a_missing_column_stops_the_run(tmp_path, column, mode):
    columns = [c + "_x" if c == column else c for c in COMBINED]
    src = write_plain(tmp_path / "in.tsv", to_bytes(columns, variant_sorted()))
    args = ["--clear"] if mode == "clear" else ["--consequence", write_cons(tmp_path), "--mode", mode]
    assert column in stamp_fails(src, *args)


def test_a_consequence_file_missing_a_column_stops_the_run(tmp_path):
    src = write_plain(tmp_path / "in.tsv", to_bytes(COMBINED, variant_sorted()))
    cons = write_cons(tmp_path, columns=CONS_COLUMNS[:5] + ["gene"])
    assert "gene_most_severe" in stamp_fails(src, "--consequence", cons)


@pytest.mark.parametrize("access", ["scan", "regions"])
def test_unsorted_input_stops_merge_naming_the_position(tmp_path, access):
    rows = variant_sorted()
    rows[4], rows[5] = rows[5], rows[4]
    src = write_plain(tmp_path / "in.tsv", to_bytes(COMBINED, rows))
    stderr = stamp_fails(src, "--consequence", write_cons(tmp_path), "--consequence-access", access)
    assert "in.tsv: not sorted at 1:200" in stderr


def test_a_chromosome_that_comes_back_stops_merge(tmp_path):
    rows = variant_sorted()
    rows.append(served(2, 100, "A", "C"))
    src = write_plain(tmp_path / "in.tsv", to_bytes(COMBINED, rows))
    assert "not sorted at 2:100" in stamp_fails(src, "--consequence", write_cons(tmp_path))


@pytest.mark.parametrize("access", ["scan", "regions"])
def test_unsorted_consequence_file_stops_merge_naming_the_position(tmp_path, access):
    cons = write_cons(tmp_path)
    # tabix refuses to index an unsorted file, so the disorder goes in under an index that
    # was built while the file was sorted
    rows = list(CONS)
    rows[4], rows[5] = rows[5], rows[4]
    write_bgz(cons, to_bytes(CONS_COLUMNS, rows))
    src = write_plain(tmp_path / "in.tsv", to_bytes(COMBINED, variant_sorted()))
    stderr = stamp_fails(src, "--consequence", cons, "--consequence-access", access)
    assert "cons.tsv.bgz: not sorted at 1:200" in stderr


def test_two_consequence_rows_for_one_variant_stop_the_run(tmp_path):
    cons = write_cons(tmp_path, [CONS[1], CONS[1], *CONS[2:]])
    src = write_plain(tmp_path / "in.tsv", to_bytes(COMBINED, variant_sorted()))
    for mode in ("merge", "lookup"):
        assert "more than one row" in stamp_fails(src, "--consequence", cons, "--mode", mode)


@pytest.mark.parametrize("mode", ["merge", "lookup"])
@pytest.mark.parametrize("pos", ["0100", "1e2", "", "-5"])
def test_a_position_that_is_not_a_plain_integer_stops_the_run(tmp_path, mode, pos):
    rows = variant_sorted()
    rows[1][6] = pos
    src = write_plain(tmp_path / "in.tsv", to_bytes(COMBINED, rows))
    stamp_fails(src, "--consequence", write_cons(tmp_path), "--mode", mode)


@pytest.mark.parametrize("mode", ["merge", "lookup", "clear"])
@pytest.mark.parametrize("change", ["longer", "shorter"])
def test_a_row_with_another_column_count_stops_the_run(tmp_path, mode, change):
    rows = variant_sorted()
    rows[3] = rows[3] + [""] if change == "longer" else rows[3][:-1]
    src = write_plain(tmp_path / "in.tsv", to_bytes(COMBINED, rows))
    args = ["--clear"] if mode == "clear" else ["--consequence", write_cons(tmp_path), "--mode", mode]
    assert "row 5" in stamp_fails(src, *args)


def test_max_lookup_rows_stops_before_the_consequence_file_is_read(tmp_path):
    src = write_plain(tmp_path / "in.tsv", to_bytes(QTL, variant_sorted(gene_pos=1000)))
    cons = write_cons(tmp_path)
    # a tabix earlier on PATH that leaves a mark, to see whether the consequence file was touched
    shim, mark = tmp_path / "bin", tmp_path / "tabix-ran"
    shim.mkdir()
    (shim / "tabix").write_text(f'#!/bin/sh\ntouch "{mark}"\nPATH="{os.environ["PATH"]}" exec tabix "$@"\n')
    (shim / "tabix").chmod(0o755)
    env = {**os.environ, "PATH": f"{shim}:{os.environ['PATH']}"}
    out = tmp_path / "out.tsv"
    base = [sys.executable, str(SCRIPT), "--input", str(src), "--output", str(out), "--consequence", str(cons),
            "--mode", "lookup", "--max-lookup-rows"]
    result = subprocess.run([*base, "8"], capture_output=True, text=True, env=env)
    assert result.returncode != 0 and "--max-lookup-rows" in result.stderr
    assert not mark.exists() and not out.exists()
    result = subprocess.run([*base, "9"], capture_output=True, text=True, env=env)
    assert result.returncode == 0 and mark.exists()


def test_the_memory_ceiling_stops_the_run():
    with pytest.raises(SystemExit, match="--max-rss-gb"):
        annotate_consequence.rss_guard(1 << 20)()
    annotate_consequence.rss_guard(1 << 40)()


def test_output_must_not_be_the_input(tmp_path):
    src = write_plain(tmp_path / "in.tsv", to_bytes(COMBINED, variant_sorted()))
    before = src.read_bytes()
    assert run("--input", src, "--output", src, "--clear").returncode != 0
    assert src.read_bytes() == before


def test_gzip_that_is_not_bgzip_is_refused(tmp_path):
    src = tmp_path / "in.tsv.gz"
    src.write_bytes(gzip.compress(to_bytes(COMBINED, variant_sorted())))
    assert "not bgzip" in stamp_fails(src, "--clear")


def verify(original: Path, stamped: Path):
    result = run("--verify", original, stamped)
    return result.returncode, json.loads(result.stdout.splitlines()[-1])


def test_verify_passes_a_stamped_file_and_counts_what_changed(tmp_path):
    src = write_bgz(tmp_path / "in.tsv.gz", to_bytes(COMBINED, variant_sorted()), "-s6", "-b7", "-e7")
    summary, out = stamp_ok(src, "--consequence", write_cons(tmp_path))
    code, report = verify(src, out)
    assert code == 0 and report["ok"] and report["header_identical"]
    assert report["rows_original"] == report["rows_stamped"] == 10
    for column in ("most_severe", "gene_most_severe"):
        assert report[column] == summary[column]


@pytest.mark.parametrize("damage", ["aaf", "header", "dropped row", "empty field filled"])
def test_verify_fails_when_anything_else_differs(tmp_path, damage):
    rows = variant_sorted()
    src = write_plain(tmp_path / "in.tsv", to_bytes(COMBINED, rows))
    columns, changed = list(COMBINED), [list(r) for r in rows]
    if damage == "aaf":
        changed[2][16] = "0.6493"
    elif damage == "header":
        columns[16] = "af"
    elif damage == "dropped row":
        changed.pop()
    else:
        changed[2][4] = "NA"
    out = write_plain(tmp_path / "out.tsv", expected(columns, changed))
    code, report = verify(src, out)
    assert code != 0 and not report["ok"]
    if damage in ("aaf", "empty field filled"):
        assert report["rows_differing_outside_annotation"] == 1 and report["first_differing_row"] == 4
