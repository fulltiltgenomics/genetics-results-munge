"""End-to-end tests of `annotate_resource.py` with directories standing in for the buckets:
that a whole resource arrives under the new prefix with its names kept, that the stats are
regenerated from the stamped rows and only after the originals were reproduced, and which
situations stop the run before the destination is touched."""

import json
import os
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

from credible_set_stats import calculate_stats, get_tsv_header, stats_to_tsv_row, write_stats_json
from test_annotate_consequence import COMBINED, QTL, to_bytes, write_bgz, write_cons, write_plain

SCRIPTS = Path(__file__).parent
COMBINED_NAME, QTL_NAME = "X_credible_sets.tsv.gz", "X_credible_sets.qtl.tsv.gz"
CONS = [
    ["1", "100", "A", "C", "missense_variant", "GENE1"],
    ["1", "200", "T", "C", "stop_gained", "GENE2"],
    ["2", "300", "G", "A", "intron_variant", "GENE3"],
]


def row(trait, chrom, pos, ref, alt, cs, pip, beta, most_severe="NA", gene="NA"):
    return ["DS", "pQTL", trait, trait, "plasma", str(chrom), str(pos), ref, alt, "11.9101", beta,
            "1.000e-02", pip, cs, "2", "NA", "1.000e-01", most_severe, gene]


# T1 loses a coding credible set and gains a LoF one; T2 keeps a row gnomAD does not hold
ROWS = [
    row("T1", 1, 100, "A", "C", "cs1", "0.9", "5.0e-01"),
    row("T2", 1, 100, "A", "C", "cs1", "0.6", "-5.0e-01", "intron_variant", "OLD"),
    row("T1", 1, 200, "T", "C", "cs2", "0.8", "5.0e-01", "intron_variant", "OLD"),
    row("T1", 2, 300, "G", "A", "cs3", "0.7", "-5.0e-01", "missense_variant", "OLD"),
    row("T2", 2, 999, "G", "A", "cs2", "0.5", "5.0e-01", "stop_gained", "GONE"),
]


def tree(root: Path) -> dict:
    return {str(p.relative_to(root)): p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()}


@pytest.fixture
def resource(tmp_path):
    """A source tree whose stats come out of the code path of the munge scripts."""
    src = tmp_path / "bucket" / "resource" / "v1"
    (src / "individual").mkdir(parents=True)
    combined = write_bgz(src / COMBINED_NAME, to_bytes(COMBINED, ROWS), "-s6", "-b7", "-e7")
    qtl_rows = [r + ["1", str(5000 - i), "6000"] for i, r in enumerate(ROWS)]
    write_bgz(src / QTL_NAME, to_bytes(QTL, sorted(qtl_rows, key=lambda r: int(r[20]))), "-s20", "-b21", "-e22")
    write_plain(src / "X_credible_sets.qtl.unmapped.tsv", to_bytes(QTL, [ROWS[4] + ["NA", "NA", "NA"]]))
    write_plain(src / "X_join.log", b"1/2: T1: OK\n")
    data = pl.read_csv(combined, separator="\t", null_values=["NA"])
    stats = []
    for (trait,), part in data.partition_by("trait", as_dict=True).items():
        part.write_csv(src / "individual" / f"{trait}.SUSIE.munged.tsv", separator="\t", null_value="NA")
        stats.append(calculate_stats(part))
        write_stats_json(stats[-1], str(src / "individual" / f"{trait}.SUSIE.munged.stats.json"))
    (src / "credible_set_stats.tsv").write_text("".join(s + "\n" for s in [get_tsv_header(), *map(stats_to_tsv_row, stats)]))
    cons_dir = tmp_path / "cons"
    cons_dir.mkdir()
    return {"src": src, "dst": tmp_path / "bucket" / "resource" / "v2", "staging": tmp_path / "staging",
            "cons": write_cons(cons_dir, CONS)}


def drive(res, *extra, dest=None, wrapper=False):
    cmd = ["bash", str(SCRIPTS / "annotate_resource.sh")] if wrapper else [sys.executable, str(SCRIPTS / "annotate_resource.py")]
    return subprocess.run(
        [*cmd, "--source", str(res["src"]), "--dest", str(dest or res["dst"]), "--combined", COMBINED_NAME,
         "--qtl", QTL_NAME, "--per-phenotype-prefix", "individual/", "--per-phenotype-suffix", ".SUSIE.munged.tsv",
         "--consequence", str(res["cons"]), "--consequence-version", "test-1", "--staging", str(res["staging"]),
         "--workers", "2", *extra],
        capture_output=True, text=True, env={**os.environ, "PYTHON": sys.executable})


def report(res) -> dict:
    lines = (res["staging"] / "report" / "annotation_report.tsv").read_text().splitlines()
    cols = lines[0].split("\t")
    return {f[0]: dict(zip(cols, f)) for f in (line.split("\t") for line in lines[1:])}


def test_whole_resource(resource):
    before = tree(resource["src"])
    p = drive(resource, wrapper=True)
    assert p.returncode == 0, p.stderr
    assert tree(resource["src"]) == before
    out = tree(resource["dst"])
    assert set(out) == set(before) | {"annotation_report.tsv", "annotation_summary.txt"}
    assert out["X_join.log"] == before["X_join.log"]
    assert "RESULT: PASS" in out["annotation_summary.txt"].decode()

    stamped = pl.read_csv(resource["dst"] / COMBINED_NAME, separator="\t", infer_schema_length=0)
    assert stamped["most_severe"].to_list() == ["missense_variant", "missense_variant", "stop_gained", "intron_variant", "NA"]
    listed = subprocess.run(["tabix", str(resource["dst"] / QTL_NAME), "1:4996-4996"], capture_output=True, text=True)
    assert listed.stdout.split("\t")[17:19] == ["NA", "NA"]
    assert out["X_credible_sets.qtl.unmapped.tsv"].split(b"\t")[-5:-3] == [b"NA", b"NA"]
    assert b"\t0.9\t" in out["individual/T1.SUSIE.munged.tsv"]

    t1 = json.loads(out["individual/T1.SUSIE.munged.stats.json"])
    was = json.loads(before["individual/T1.SUSIE.munged.stats.json"])
    assert (was["n_risk_cs_with_coding"], was["n_risk_cs_with_lof"], was["n_protective_cs_with_coding"]) == (0, 0, 1)
    assert (t1["n_risk_cs_with_coding"], t1["n_risk_cs_with_lof"], t1["n_protective_cs_with_coding"]) == (2, 1, 0)
    aggregate = out["credible_set_stats.tsv"].decode().splitlines()
    assert aggregate[0] == get_tsv_header()
    assert [line.split("\t")[0] for line in aggregate[1:]] == [line.split("\t")[0] for line in before["credible_set_stats.tsv"].decode().splitlines()[1:]]
    assert aggregate[1].split("\t")[4:9] == ["2", "2", "2", "1", "1"]

    rows = report(resource)
    assert all(r["ok"] == "True" for r in rows.values())
    assert rows[COMBINED_NAME]["rows_in"] == rows[COMBINED_NAME]["rows_out"] == "5"
    assert rows[COMBINED_NAME]["most_severe_changed_non_na"] == "4"
    assert rows["credible_set_stats.tsv"]["original_reproduced"] == "True"


def test_dry_run_touches_nothing(resource):
    p = drive(resource, "--dry-run")
    assert p.returncode == 0, p.stderr
    plan = {line.split("\t")[2].rsplit("/", 1)[-1]: line.split("\t") for line in p.stdout.splitlines()[1:]}
    assert plan[COMBINED_NAME][:2] == ["combined", "stamp-merge"]
    assert plan[QTL_NAME][:2] == ["qtl", "stamp-lookup"]
    assert plan["X_credible_sets.qtl.unmapped.tsv"][:2] == ["qtl_unmapped", "stamp-lookup"]
    assert plan["T1.SUSIE.munged.tsv"][3] == f"{resource['dst']}/individual/T1.SUSIE.munged.tsv"
    assert plan[f"{COMBINED_NAME}.tbi"][1] == plan["credible_set_stats.tsv"][1] == "regenerate"
    assert plan["X_join.log"][1] == "copy"
    assert not resource["dst"].exists() and not resource["staging"].exists()


def test_rerun_resumes_and_leaves_the_destination_alone(resource):
    assert drive(resource).returncode == 0
    state = resource["staging"] / "state" / "individual"
    (state / "T2.SUSIE.munged.tsv.json").unlink()
    (resource["dst"] / "individual" / "T2.SUSIE.munged.tsv").unlink()
    kept = (state / "T1.SUSIE.munged.tsv.json").stat().st_mtime_ns
    written = {p: p.stat().st_mtime_ns for p in resource["dst"].rglob("*") if p.is_file()}
    p = drive(resource)
    assert p.returncode == 0, p.stderr
    assert "5 files to stamp, 4 already done" in p.stderr
    assert (state / "T1.SUSIE.munged.tsv.json").stat().st_mtime_ns == kept
    assert (resource["dst"] / "individual" / "T2.SUSIE.munged.tsv").exists()
    assert all(p.stat().st_mtime_ns == t for p, t in written.items())


@pytest.mark.parametrize("name", ["notes.txt", "individual/T1.SUSIE.munged.null_traits.tsv", "stray.tsv.gz.tbi"])
def test_unclassified_object_stops_the_run(resource, name):
    (resource["src"] / name).write_text("?\n")
    for extra in ([], ["--dry-run"]):
        p = drive(resource, *extra)
        assert p.returncode != 0 and name in p.stderr
    assert not resource["dst"].exists()
    assert drive(resource, "--copy", name).returncode == (1 if name.endswith(".tbi") else 0)


def test_skip_rule(resource):
    (resource["src"] / "notes.txt").write_text("?\n")
    assert drive(resource, "--skip", "*.txt").returncode == 0
    assert not (resource["dst"] / "notes.txt").exists()
    assert report(resource)["notes.txt"]["action"] == "skip"


def test_copy_refuses_a_file_with_annotation_columns(resource):
    (resource["src"] / "extra.tsv").write_bytes(to_bytes(COMBINED, ROWS[:1]))
    p = drive(resource, "--copy", "extra.tsv")
    assert p.returncode != 0 and "extra.tsv" in p.stderr
    assert not resource["dst"].exists()


def test_stats_that_cannot_be_reproduced_stop_the_run(resource):
    stats = resource["src"] / "credible_set_stats.tsv"
    stats.write_text(stats.read_text().replace("\t1\t", "\t7\t", 1))
    p = drive(resource)
    assert p.returncode == 1 and "not reproduced" in p.stderr
    assert not resource["dst"].exists()
    assert "FAIL  original stats reproduced" in (resource["staging"] / "report" / "annotation_summary.txt").read_text()


def test_destination_guards(resource):
    for dest in (resource["src"], resource["src"] / "stamped", resource["src"].parent):
        p = drive(resource, dest=dest)
        assert p.returncode != 0 and "must not contain" in p.stderr
    resource["dst"].mkdir()
    (resource["dst"] / "someone_elses.tsv").write_text("x\n")
    p = drive(resource)
    assert p.returncode != 0 and "would not write" in p.stderr
    (resource["dst"] / "someone_elses.tsv").rename(resource["dst"] / "X_join.log")
    p = drive(resource)
    assert p.returncode != 0 and "other content" in p.stderr
    assert (resource["dst"] / "X_join.log").read_text() == "x\n"
    assert len(list(resource["dst"].iterdir())) == 1


def test_staging_is_bound_to_its_arguments(resource):
    assert drive(resource).returncode == 0
    p = drive(resource, "--consequence-version", "test-2", dest=resource["dst"].parent / "v3")
    assert p.returncode != 0 and "other arguments" in p.stderr
