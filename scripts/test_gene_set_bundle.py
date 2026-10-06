"""Offline tests of `gene_set_bundle.py`: the Entrez and symbol mapping rules, per-set unmapped
accounting, the mapping-rate refusal, set id uniqueness, checksums, and the staging guard."""

import json
import subprocess

import polars as pl
import pytest

import gene_set_bundle as gsb
from build_gene_sets import sha256_of, write_manifest

HGNC_COLS = ["hgnc_id", "symbol", "locus_group", "locus_type", "status", "entrez_id",
             "ensembl_gene_id", "prev_symbol", "alias_symbol"]
HGNC_ROWS = [
    # hgnc_id, symbol, entrez, ensembl, prev, alias
    ("HGNC:1", "AAA", "1", "ENSG01", "", ""),
    ("HGNC:2", "BBB", "2", "ENSG02", "OLDB", '"SHARED|ALIASB"'),
    ("HGNC:3", "CCC", "3", "ENSG03", "", "SHARED"),
    ("HGNC:4", "NOENS", "4", "", "", ""),
    ("HGNC:5", "DUP1", "5", "ENSG05", "", ""),
    ("HGNC:6", "DUP2", "6", "ENSG05", "", ""),
    # a prev symbol that is another gene's approved symbol: the approved symbol wins
    ("HGNC:7", "DDD", "7", "ENSG07", "AAA", ""),
    # TRM1's shape: a prev symbol of two genes, neither with an Ensembl id, that is also a
    # unique alias of a third gene that has one; the prev tier decides it, ambiguously
    ("HGNC:20", "NOENSA", "20", "", "TRMX", ""),
    ("HGNC:21", "NOENSB", "21", "", "TRMX", ""),
    ("HGNC:22", "EEE", "22", "ENSG22", "", "TRMX|PNOENS"),
    # a unique prev symbol of a gene with no Ensembl id, also EEE's alias above: unmapped
    ("HGNC:23", "NOENSC", "23", "", "PNOENS", ""),
    # an alias that is another gene's approved symbol: the approved symbol wins
    ("HGNC:24", "FFF", "24", "ENSG24", "", "CCC"),
    # a lower-case alias spelling an approved symbol in another case: matching is exact
    ("HGNC:25", "GGG", "25", "ENSG25", "", "aaa"),
]
GENCODE = ("gene_id\tchrom\tgene_start\tgene_end\tgene_strand\tgene_name\tgene_type\n"
           "ENSG01.3\t1\t100\t200\t+\tAAA\tprotein_coding\n"
           "ENSG02.1\t23\t300\t400\t-\tBBB\tprotein_coding\n")


def write_inputs(tmp_path):
    lines = ["\t".join(HGNC_COLS)]
    for hid, sym, entrez, ens, prev, alias in HGNC_ROWS:
        lines.append("\t".join([hid, sym, "protein-coding gene", "gene with protein product",
                                "Approved", entrez, ens, prev, alias]))
    hgnc = tmp_path / "hgnc_complete_set.txt"
    hgnc.write_text("\n".join(lines) + "\n")
    gencode = tmp_path / "gencode.v49.annotation.genes.tsv"
    gencode.write_text(GENCODE)
    return hgnc, gencode


def write_raw(tmp_path, msigdb_text, reactome_text):
    raw = tmp_path / "gene_sets.20261006"
    raw.mkdir()
    (raw / "h.gmt").write_text(msigdb_text)
    (raw / "r.zip").write_bytes(b"zip")
    (raw / "ReactomePathways.gmt").write_text(reactome_text)
    write_manifest(raw / "manifest.json", {"release_id": "20261006", "collections": [
        {"collection": "h.all", "source": "msigdb", "file": "h.gmt",
         "sha256": sha256_of(raw / "h.gmt"), "spdx_license": "CC-BY-4.0", "attribution": "A"},
        {"collection": "reactome", "source": "reactome", "file": "r.zip",
         "sha256": sha256_of(raw / "r.zip"), "gmt_file": "ReactomePathways.gmt",
         "gmt_sha256": sha256_of(raw / "ReactomePathways.gmt"), "spdx_license": "CC0-1.0",
         "attribution": "R"},
    ]})
    return raw


MSIGDB = "HALLMARK_A\thttp://u/A\t1\t2\t2\t4\nHALLMARK_B\thttp://u/B\t8\t5\t6\n"
REACTOME = ("Path one\tR-HSA-1\tAAA\tOLDB\tALIASB\tSHARED\tgag\tNOENS\n"
            "Path two\tR-HSA-2\tBBB\tCCC\tDDD\n")


def build(tmp_path, msigdb=MSIGDB, reactome=REACTOME):
    hgnc, gencode = write_inputs(tmp_path)
    raw = write_raw(tmp_path, msigdb, reactome)
    return gsb.build_bundle(raw, hgnc, gencode, tmp_path / "out", hgnc_version="2026-06-02",
                            gencode_version="49")


@pytest.fixture
def floors(monkeypatch):
    # the fixtures map far worse than real data, on purpose; refusal is tested separately
    monkeypatch.setattr(gsb, "MIN_MAPPING_RATE", {"msigdb": 0.0, "reactome": 0.0})


def test_mapping_rules_and_accounting(tmp_path, floors):
    bundle = build(tmp_path)
    sets = {r["set_id"]: r for r in pl.read_parquet(bundle / "sets.parquet").to_dicts()}
    members = pl.read_parquet(bundle / "members.parquet")
    by_set = {s: sorted(members.filter(pl.col("set_id") == s)["ensembl_gene_id"]) for s in sets}

    # entrez: 4 has no Ensembl id; the duplicated 2 is one given member
    assert by_set["HALLMARK_A"] == ["ENSG01", "ENSG02"]
    assert (sets["HALLMARK_A"]["n_given"], sets["HALLMARK_A"]["n_mapped"],
            sets["HALLMARK_A"]["n_unmapped"]) == (3, 2, 1)
    # 5 and 6 share ENSG05 in HGNC, which is ambiguous and so unmapped
    assert by_set["HALLMARK_B"] == []
    assert sets["HALLMARK_B"]["n_unmapped"] == 3

    # symbol: approved, unique prev, unique alias map; SHARED (alias of two genes), the
    # lower-case viral gene and the gene without an Ensembl id do not; OLDB and ALIASB both
    # land on BBB, so they collapse
    assert by_set["R-HSA-1"] == ["ENSG01", "ENSG02"]
    assert (sets["R-HSA-1"]["n_given"], sets["R-HSA-1"]["n_mapped"],
            sets["R-HSA-1"]["n_unmapped"]) == (6, 2, 3)
    assert by_set["R-HSA-2"] == ["ENSG02", "ENSG03", "ENSG07"]
    assert sets["R-HSA-1"]["name"] == "Path one"
    assert sets["R-HSA-1"]["url_or_stable_id"] == "R-HSA-1"

    assert set(members["symbol"]) == {"AAA", "BBB", "CCC", "DDD"}
    genes = pl.read_parquet(bundle / "genes.parquet")
    assert "ENSG05" not in set(genes["ensembl_gene_id"])
    assert genes.filter(pl.col("ensembl_gene_id") == "ENSG02")["chr"].to_list() == [23]
    assert genes.filter(pl.col("ensembl_gene_id") == "ENSG03")["chr"].to_list() == [None]


def resolvers_for(tmp_path):
    hgnc_path, gencode_path = write_inputs(tmp_path)
    hgnc = gsb.read_hgnc(hgnc_path)
    return gsb.build_resolvers(hgnc, gsb.gene_universe(hgnc, gsb.read_gencode(gencode_path)))


def test_prev_symbol_never_overrides_an_approved_symbol(tmp_path):
    symbol = resolvers_for(tmp_path)["symbol"]
    assert symbol["AAA"] == "ENSG01"
    assert "SHARED" not in symbol


def test_first_tier_that_names_a_symbol_decides_it(tmp_path):
    symbol = resolvers_for(tmp_path)["symbol"]
    # ambiguous at the prev tier, which hides EEE's alias
    assert "TRMX" not in symbol
    # one gene at the prev tier, but it has no Ensembl id
    assert "PNOENS" not in symbol
    assert symbol["CCC"] == "ENSG03"
    assert symbol["aaa"] == "ENSG25"
    assert "Aaa" not in symbol


def test_first_tier_rule_end_to_end(tmp_path, floors):
    bundle = build(tmp_path, reactome=REACTOME + "Path three\tR-HSA-3\tTRMX\tPNOENS\tCCC\taaa\tAaa\n")
    sets = {r["set_id"]: r for r in pl.read_parquet(bundle / "sets.parquet").to_dicts()}
    members = pl.read_parquet(bundle / "members.parquet").filter(pl.col("set_id") == "R-HSA-3")
    assert sorted(members["ensembl_gene_id"]) == ["ENSG03", "ENSG25"]
    assert (sets["R-HSA-3"]["n_mapped"], sets["R-HSA-3"]["n_unmapped"]) == (2, 3)


def test_manifest_records_checksums_and_mapping(tmp_path, floors):
    bundle = build(tmp_path)
    manifest = gsb.verify_bundle(bundle)
    meta = manifest["bundle"]
    assert meta["release_id"] == "20261006.hgnc20260602.gencode49"
    assert bundle.name == "gene_sets_bundle.20261006.hgnc20260602.gencode49"
    for name, entry in meta["files"].items():
        assert entry["sha256"] == sha256_of(bundle / name)
    assert meta["hgnc"]["sha256"] == sha256_of(tmp_path / "hgnc_complete_set.txt")
    assert meta["mapping"]["h.all"]["member_pairs"] == 6
    assert meta["mapping"]["h.all"]["mapped_pairs"] == 2
    # licence and attribution are carried over from the raw manifest unchanged
    raw = json.loads((tmp_path / "gene_sets.20261006" / "manifest.json").read_text())
    assert manifest["collections"] == raw["collections"]
    assert sorted(p.name for p in bundle.iterdir()) == sorted([*gsb.BUNDLE_FILES, "manifest.json"])

    (bundle / "sets.parquet").write_bytes(b"tampered")
    with pytest.raises(SystemExit, match="sha256"):
        gsb.verify_bundle(bundle)


def test_mapping_rate_below_floor_refuses_and_leaves_nothing(tmp_path):
    with pytest.raises(SystemExit, match="h.all 0.3333 < 0.99"):
        build(tmp_path)
    assert not (tmp_path / "out").exists() or not any((tmp_path / "out").iterdir())


def test_reactome_mapping_rate_below_floor_refuses(tmp_path, monkeypatch):
    # fixture Reactome maps 6 of 9 member pairs
    monkeypatch.setattr(gsb, "MIN_MAPPING_RATE", {"msigdb": 0.0, "reactome": 0.7})
    with pytest.raises(SystemExit, match=r"reactome 0\.6667 < 0\.7"):
        build(tmp_path)


def test_mapping_rate_exactly_at_floor_passes(tmp_path, monkeypatch):
    # h.all maps 2 of 6 pairs; the floor is a minimum, so a rate equal to it passes and one
    # just above it refuses
    floor = 2 / 6
    assert floor * 6 == 2
    monkeypatch.setattr(gsb, "MIN_MAPPING_RATE", {"msigdb": floor, "reactome": 0.0})
    build(tmp_path)
    monkeypatch.setattr(gsb, "MIN_MAPPING_RATE", {"msigdb": 0.3334, "reactome": 0.0})
    (tmp_path / "again").mkdir()
    with pytest.raises(SystemExit, match="h.all"):
        build(tmp_path / "again")


def test_refused_source_with_valid_checksum_refused(tmp_path, floors):
    # the raw manifest's sha256 matches, so only this stage's own re-run of the guard stops it
    with pytest.raises(SystemExit, match="refused source"):
        build(tmp_path, reactome=REACTOME + "KEGG_X\tR-HSA-9\tAAA\n")


def test_duplicate_set_id_refused(tmp_path, floors):
    with pytest.raises(SystemExit, match="not unique"):
        build(tmp_path, reactome="One\tR-HSA-1\tAAA\nTwo\tR-HSA-1\tBBB\n")


def test_tampered_raw_file_refused(tmp_path, floors):
    hgnc, gencode = write_inputs(tmp_path)
    raw = write_raw(tmp_path, MSIGDB, REACTOME)
    (raw / "ReactomePathways.gmt").write_text(REACTOME + "KEGG_X\tR-HSA-9\tAAA\n")
    with pytest.raises(SystemExit, match="sha256"):
        gsb.build_bundle(raw, hgnc, gencode, tmp_path / "out", hgnc_version="2026-06-02",
                         gencode_version="49")


def test_existing_bundle_not_overwritten(tmp_path, floors):
    bundle = build(tmp_path)
    hgnc, gencode = tmp_path / "hgnc_complete_set.txt", tmp_path / "gencode.v49.annotation.genes.tsv"
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        gsb.build_bundle(tmp_path / "gene_sets.20261006", hgnc, gencode, bundle.parent,
                         hgnc_version="2026-06-02", gencode_version="49")


TARGET = "gs://daly-genetics-results/gene_sets/20261006.hgnc20260602.gencode49/"


def fake_run(listing_stdout="", listing_rc=1, stderr="One or more URLs matched no objects.",
             fail_copy_of=None):
    calls = []

    def run(cmd, **kwargs):
        calls.append(cmd)
        if cmd[:3] == ["gcloud", "storage", "ls"]:
            return subprocess.CompletedProcess(cmd, listing_rc, listing_stdout, stderr)
        rc = 1 if fail_copy_of and cmd[-1].endswith(fail_copy_of) else 0
        return subprocess.CompletedProcess(cmd, rc, "", "")
    return run, calls


def test_stage_dry_run_copies_nothing(tmp_path, floors, capsys):
    bundle = build(tmp_path)
    run, calls = fake_run()
    target = gsb.stage(bundle, dry_run=True, run=run)
    assert target == TARGET
    assert calls == []
    assert target in capsys.readouterr().out


def test_stage_refuses_tampered_bundle(tmp_path, floors):
    bundle = build(tmp_path)
    (bundle / "members.parquet").write_bytes(b"tampered")
    run, calls = fake_run()
    with pytest.raises(SystemExit, match="sha256"):
        gsb.stage(bundle, dry_run=False, run=run)
    assert calls == []


def test_stage_refuses_existing_remote(tmp_path, floors):
    bundle = build(tmp_path)
    listing = "".join(f"{TARGET}{n}\n" for n in (*gsb.BUNDLE_FILES, "manifest.json"))
    run, calls = fake_run(listing_stdout=listing, listing_rc=0, stderr="")
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        gsb.stage(bundle, dry_run=False, run=run)
    assert len(calls) == 1


def test_stage_refuses_partial_remote(tmp_path, floors):
    bundle = build(tmp_path)
    run, calls = fake_run(listing_stdout=f"{TARGET}sets.parquet\n", listing_rc=0, stderr="")
    with pytest.raises(SystemExit, match=f"partial upload.*gcloud storage rm -r {TARGET}"):
        gsb.stage(bundle, dry_run=False, run=run)
    assert len(calls) == 1


def test_stage_stops_at_first_failed_copy(tmp_path, floors):
    bundle = build(tmp_path)
    run, calls = fake_run(fail_copy_of="members.parquet")
    with pytest.raises(SystemExit, match="copy failed"):
        gsb.stage(bundle, dry_run=False, run=run)
    assert [c[-1] for c in calls[1:]] == [TARGET + "sets.parquet", TARGET + "members.parquet"]


def test_stage_refuses_when_listing_fails(tmp_path, floors):
    bundle = build(tmp_path)
    run, calls = fake_run(stderr="403 permission denied")
    with pytest.raises(SystemExit, match="could not list"):
        gsb.stage(bundle, dry_run=False, run=run)
    assert len(calls) == 1


def test_stage_copies_to_empty_remote(tmp_path, floors):
    bundle = build(tmp_path)
    run, calls = fake_run()
    gsb.stage(bundle, dry_run=False, run=run)
    copies = calls[1:]
    # one create-only copy per object, manifest.json last
    names = [*gsb.BUNDLE_FILES, "manifest.json"]
    assert copies == [["gcloud", "storage", "cp", "--if-generation-match=0", str(bundle / n),
                       TARGET + n] for n in names]
