"""Offline tests of `munge_protein_structure.py` on a synthetic bundle: window ownership and
is_owner, the partial-PAE refusal, PAE encoding and symmetrisation, headers and sort orders,
the variant address rule, model-less accessions, features, predictor scores, and the pointer.

Needs h5py, hdf5plugin and pyarrow (the protein-structure venv), and bgzip/tabix on PATH."""

import gzip
import hashlib
import json
import os
import subprocess
import types
import zipfile

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from pathlib import Path

import munge_protein_structure as mps
import protein_structure_contract as psc

THREE = {v: k for k, v in mps.THREE_TO_ONE.items()}
SEQ = {"Q00001": "MKTAYIAKQR", "Q00002": "MSTNP", "Q00003": "MAAAG", "Q00004": "MCCCC",
       "Q00006": "MSSSG", "Q00007": "MAUAG",
       "Q00010": "".join(np.random.default_rng(7).choice(list("ACDEFGHIKLMNPQRSTVWY"), 230))}
# accession -> (pdb files with guide range, af3, pae)
MODELS = {
    "Q00001": ([("AF-Q00001-F1-model_v6", 1, 6), ("AF-Q00001-F2-model_v6", 3, 10)], False, True),
    "Q00002": ([("AF3-Q00002-F1-model_v6", 1, 5)], True, True),
    "Q00003": ([("AF-Q00003-F1-model_v6", 1, 5)], False, False),
    "Q00006": ([("AF-Q00006-F1-model_v6", 1, 5)], False, False),
    "Q00007": ([("AF-Q00007-F1-model_v6", 1, 5)], False, False),
    # eleven windows, so string order (F1, F10, F11, F2, ...) differs from window order
    "Q00010": ([(f"AF-Q00010-F{k}-model_v6", 1 + 20 * (k - 1), 30 + 20 * (k - 1)) for k in range(1, 12)], False, False),
}


def pdb_text(seq, af3, rng):
    head = "" if af3 else ("TITLE     ALPHAFOLD MONOMER V2.0 PREDICTION FOR TEST\n"
                           "REMARK   1 ACADEMIC AND COMMERCIAL PURPOSES, UNDER CC-BY 4.0 LICENCE.\n")
    lines = []
    for i, aa in enumerate(seq, 1):
        x, y, z = rng.normal(0, 10, 3)
        for name in (" N  ", " CA "):
            lines.append(f"ATOM  {len(lines) + 1:5d} {name} {THREE[aa]} A{i:4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{50 + i:6.2f}           C  ")
    return head + "\n".join(lines) + "\nEND\n"


def pae_matrix(n, af3, rng):
    m = rng.integers(0, 320, (n, n)) / 10 if af3 else rng.integers(0, 32, (n, n)).astype(float)
    return m


def make_bundle(tmp, drop_pae=None):
    rng = np.random.default_rng(1)
    b = tmp / "bundle"
    (b / "pdb_files").mkdir(parents=True)
    (b / "pae_files").mkdir()
    guide, paes = [], {}
    for acc, (files, af3, has_pae) in MODELS.items():
        for mid, s, e in files:
            model_seq = SEQ[acc][s - 1:e].replace("U", "C")
            (b / "pdb_files" / f"{mid}.pdb.gz").write_bytes(gzip.compress(pdb_text(model_seq, af3, rng).encode()))
            pae_name = mid.replace("model", "predicted_aligned_error") + ".json.gz"
            if has_pae:
                m = pae_matrix(e - s + 1, af3, rng)
                paes[mid] = m
                if pae_name != drop_pae:
                    body = json.dumps([{"predicted_aligned_error": m.tolist() if af3 else m.astype(int).tolist(),
                                        "max_predicted_aligned_error": 31.75}])
                    (b / "pae_files" / pae_name).write_bytes(gzip.compress(body.encode()))
            guide.append((f"{mid}.pdb.gz", pae_name if has_pae else "", f"[{s}, {e}]", acc))
    guide.append(("AF-Q00004-F1-model_v6.pdb.gz", "", "[1, 5]", "Q00004"))  # no file in the bundle
    pd.DataFrame(guide, columns=["pdb_filename", "pae_filename", "pos_covered", "uniprot_id"]).to_csv(
        b / "pdb_pae_file_pos_guide.tsv", sep="\t", index=False)
    seq_rows = [
        ("Q00001", "Q00001", "ENST1", "ENSG00000000001", "GENEA", "True", "UniProt", "True", SEQ["Q00001"], "True"),
        ("Q00002", "Q00002-2", "ENST2", "ENSG00000000002", "GENEB", "False", "UniProt", "False", SEQ["Q00002"], "True"),
        ("Q00002", "Q00002-2", "ENST3", "ENSG00000000003", "GENEC", "False", "UniProt", "False", SEQ["Q00002"], "True"),
        ("Q00003", "Q00003", "ENST4", "ENSG00000000099", "OLDNAME", "True", "UniProt", "True", SEQ["Q00003"], "False"),
        ("Q00004", "Q00004", "ENST5", "ENSG00000000004", "GENED", "True", "UniProt", "True", SEQ["Q00004"], "False"),
        ("Q00006", "Q00006", "ENST6", "ENSG00000000006", "GENEF", "True", "Ensembl", "False", SEQ["Q00006"], "False"),
        ("Q00007", "Q00007", "ENST7", "ENSG00000000007", "GENEG", "True", "UniProt", "True", SEQ["Q00007"], "False"),
        ("Q00010", "Q00010", "ENST10", "ENSG00000000010", "GENEJ", "True", "UniProt", "True", SEQ["Q00010"], "False"),
    ]
    pd.DataFrame(seq_rows, columns=["uniprot_id", "uniprot_isoform", "enst", "ensg", "gene_name", "AF2DB",
                                    "sequence_source", "uniprot_canonical", "sequence", "PAE_files"]).to_csv(
        b / "protein_sequence_guide.tsv", sep="\t", index=False)

    # variants: (pos, ref, alt, aa_ref, aa_alt, aa_pos, aa_pos_file, map filename, accession)
    v = [
        (1000, "A", "G", "Y", "C", 5, 5, "AF-Q00001-F1-model_v6.pdb.gz", "Q00001"),  # map picks F1, owner is F2
        (1003, "C", "T", "Q", "R", 9, 7, "AF-Q00001-F2-model_v6.pdb.gz", "Q00001"),
        (2000, "G", "A", "S", "N", 2, 2, "AF3-Q00002-F1-model_v6.pdb.gz", "Q00002"),
        (2000, "G", "A", "T", "M", 4, 2, "AF-Q00001-F2-model_v6.pdb.gz", "Q00001"),  # same variant, 2nd protein
        (5000, "T", "C", "C", "R", 2, 2, "AF-Q00004-F1-model_v6.pdb.gz", "Q00004"),  # no model: dropped
    ]
    with h5py.File(b / "all_missense_variants_gr38.h5", "w") as h:
        for c in mps.CHROMS:
            rows = v if c == "1" else []
            h[f"chr{c}_ref_alt"] = np.array([r[1:5] for r in rows], dtype="S1").reshape(-1, 4)
            h[f"chr{c}_pos"] = np.array([(r[0], r[5], r[6]) for r in rows], dtype=np.int64).reshape(-1, 3)
            h[f"chr{c}_filename"] = np.array([r[7] for r in rows], dtype="S33").reshape(-1, 1)
            h[f"chr{c}_uniprot_id"] = np.array([r[8] for r in rows], dtype="S10").reshape(-1, 1)

    lab = tmp / "lab"
    lab.mkdir()
    pq.write_table(pa.table({"chrom": ["chr1", "chr1"], "pos": [1000, 2000], "ref": ["A", "G"],
                             "alt": ["G", "A"], "mu_snp": [1.5e-9, 2.5e-9]}), lab / "gnomad_mu_snp_by_variant.parquet")
    pq.write_table(pa.table({"uniprot_id": ["Q00001", "Q00001"], "aa_pos": pa.array([4, 5], pa.int32()),
                             "mu": [1e-8, 2e-8], "n_with_mu": [6, 7], "n_possible": [6, 8]}),
                   lab / "mu_per_residue_gnomad.parquet")
    feats = [
        ("Q00001", "Q00001-1", "Domain", "2.0", "6.0", "Kinase", "", "ECO:0000255|PROSITE-ProRule|PRU1;ECO:0000269|PubMed|1"),
        ("Q00001", "Q00001-1", "Chain", "", "10.0", "unknown start", "", ""),
        ("Q00001", "Q00001-1", "Site", "11.0", "11.0", "beyond the sequence", "", "ECO:0000250|UniProtKB|P1"),
        ("Q00001", "Q00001-1", "Region", "1.0", "3.0", "odd code", "", "ECO:0000999||"),
        ("Q00002", "Q00002-1", "Domain", "1.0", "5.0", "other isoform", "", "ECO:0000255||"),
        # same start, end and source: only feature_type, then description, orders these
        ("Q00001", "Q00001-1", "Site", "3.0", "3.0", "b-desc", "", ""),
        ("Q00001", "Q00001-1", "Region", "3.0", "3.0", "a", "", ""),
        ("Q00001", "Q00001-1", "Binding site", "4.0", "4.0", "zinc", "", ""),
        ("Q00001", "Q00001-1", "Binding site", "4.0", "4.0", "calcium", "", ""),
        # bundle sequence not UniProt's (uniprot_canonical False), bare accession
        ("Q00006", "Q00006-1", "Domain", "1.0", "3.0", "not canonical", "", ""),
        # bundle's bare accession is isoform 1, UniProt's canonical is isoform 3
        ("Q00007", "Q00007-3", "Domain", "1.0", "3.0", "other canonical", "", ""),
    ]
    pd.DataFrame(feats, columns=["uniprot_id", "canonical_isoform", "feature", "position_start", "position_end",
                                 "description", "label", "evidence"]).to_csv(
        lab / "uniprot_all_feats_evidence.tsv", sep="\t", index=False)

    gencode = tmp / "gencode.tsv"
    pd.DataFrame([
        ("ENSG00000000001.3", 1, 900, 1100, "+", "GENEA", "protein_coding"),
        ("ENSG00000000002.1", 1, 1900, 2100, "+", "GENEB", "protein_coding"),
        ("ENSG00000000003.1", 1, 1950, 2500, "+", "GENEC", "protein_coding"),
        ("ENSG00000000004.1", 1, 4900, 5100, "+", "GENED", "protein_coding"),
        ("ENSG00000000006.1", 1, 3000, 3100, "+", "GENEF", "protein_coding"),
        ("ENSG00000000007.1", 1, 3200, 3300, "+", "GENEG", "protein_coding"),
        ("ENSG00000000010.1", 1, 6000, 7000, "+", "GENEJ", "protein_coding"),
    ], columns=["gene_id", "chrom", "gene_start", "gene_end", "gene_strand", "gene_name", "gene_type"]).to_csv(
        gencode, sep="\t", index=False)

    def md5(s):
        return hashlib.md5(s.encode()).hexdigest()
    ted_rows = [
        ("AF-Q00001-F1-model_v4_TED01", md5(SEQ["Q00001"][0:3] + SEQ["Q00001"][6:9]), "high", "1-3_7-9", "6", "2",
         "90", "1", "1", "0", "1", "0", "proteome-tax_id-9606-0_v4", "3.40.50.300", "H", "foldseek", "11", "0.3",
         "Human", "Homo_sapiens", "x"),
        ("AF-Q00001-F1-model_v4_TED02", "0" * 32, "medium", "4-6", "3", "1", "90", "1", "1", "0", "1", "0",
         "proteome-tax_id-9606-0_v4", "-", "-", "-", "11", "0.3", "Human", "Homo_sapiens", "x"),
        # an AFDB v4 fragment: positions are fragment-relative, 200 residues per step
        ("AF-Q00010-F2-model_v4_TED01", md5(SEQ["Q00010"][204:210]), "high", "5-10", "6", "1", "90", "1", "1",
         "0", "1", "0", "proteome-tax_id-9606-0_v4", "-", "-", "-", "11", "0.3", "Human", "Homo_sapiens", "x"),
    ]
    ted = tmp / "ted.tsv"
    ted.write_text("\n".join("\t".join(r) for r in ted_rows) + "\n")
    make_predictors(tmp)
    return b, lab, gencode, ted, paes


# what each synthetic source scored. AlphaMissense: Q00006 on another residue at 2, Q00007 with no
# row at its selenocysteine, Q00010 one residue longer, and two accessions the release lacks.
# ESM1b: Q00002 under its isoform's name beside a canonical file with another sequence, Q00007
# with the selenocysteine read as C, Q00010 one residue short, nothing for Q00003 and Q00006
AM_SCORED = {"Q00001": SEQ["Q00001"], "Q00002": SEQ["Q00002"], "Q00003": SEQ["Q00003"], "Q00004": SEQ["Q00004"],
             "Q00006": "MTSSG", "Q00007": SEQ["Q00007"], "Q00010": SEQ["Q00010"] + "A", "Q99999": "MAAA"}
ESM_SCORED = {"Q00001": SEQ["Q00001"], "Q00002": "MSTNA", "Q00002-2": SEQ["Q00002"], "Q00007": "MACAG",
              "Q00010": SEQ["Q00010"][:-1], "Q88888": "MAAA"}
AM_FIXED = {("Q00001", 5, "C"): 0.2, ("Q00002", 2, "N"): 0.2, ("Q00001", 9, "R"): 0.9}
ESM_FIXED = {("Q00001", 5, "C"): -3.0, ("Q00001", 9, "R"): -3.0, ("Q00002", 2, "N"): -1.0}


def am_score(acc, pos, alt):
    return AM_FIXED.get((acc, pos, alt), round(0.05 + 0.04 * mps.AA.index(alt) + 0.001 * pos, 4))


def esm_score(acc, pos, alt):
    return ESM_FIXED.get((acc, pos, alt), round(-0.5 * mps.AA.index(alt) - 0.01 * pos, 3))


def write_alphamissense(tmp, extra=()):
    lines = ["# Copyright", "#", "# Licensed", "uniprot_id\tprotein_variant\tam_pathogenicity\tam_class"]
    for acc, seq in AM_SCORED.items():
        for pos, ref in enumerate(seq, 1):
            if ref != "U":
                lines += [f"{acc}\t{ref}{pos}{alt}\t{am_score(acc, pos, alt)}\tambiguous" for alt in mps.AA if alt != ref]
    lines += list(extra)
    path = tmp / "predictors" / "alphamissense" / "AlphaMissense_aa_substitutions.tsv.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode(), mtime=0))
    return {path.name: hashlib.sha256(path.read_bytes()).hexdigest()}


def write_esm1b(tmp):
    path = tmp / "predictors" / "esm1b" / "ALL_hum_isoforms_ESM1b_LLR.zip"
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w") as z:
        for iso, seq in ESM_SCORED.items():
            rows = ["," + ",".join(f"{aa} {i}" for i, aa in enumerate(seq, 1))]
            rows += [alt + "," + ",".join("0.000" if alt == aa else f"{esm_score(iso.split('-')[0], i, alt):.3f}"
                                          for i, aa in enumerate(seq, 1)) for alt in "KRHEDNQTSCGAVLIMPYFW"]
            z.writestr(zipfile.ZipInfo(f"content/ALL_hum_isoforms_ESM1b_LLR/{iso}_LLR.csv"), "\n".join(rows) + "\n")
    return {path.name: hashlib.sha256(path.read_bytes()).hexdigest()}


def make_predictors(tmp):
    # the build pins each source file's sha256; the synthetic files are pinned as they are written
    mps.SCORE_SOURCES["alphamissense"]["files"] = write_alphamissense(tmp)
    mps.SCORE_SOURCES["esm1b"]["files"] = write_esm1b(tmp)


def args_for(tmp, stages=",".join(mps.STAGES), extra=()):
    return ["--bundle", str(tmp / "bundle"), "--lab-files", str(tmp / "lab"), "--gencode", str(tmp / "gencode.tsv"),
            "--workdir", str(tmp / "build"), "--ted", str(tmp / "ted.tsv"), "--ted-meta", str(tmp / "ted_meta.json"),
            "--predictors", str(tmp / "predictors"),
            "--procs", "1", "--release", "2026-01-01-000000", "--stages", stages, *extra]


def run(tmp, stages=",".join(mps.STAGES), **kw):
    _, _, _, _, paes = make_bundle(tmp, **kw)
    (tmp / "ted_meta.json").write_text(json.dumps({"version": "test"}))
    return mps.main(args_for(tmp, stages)), paes


def read_tsv(path):
    with gzip.open(path, "rt") as f:
        header = f.readline()
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment=None, skiprows=1, header=None)
    return header, df


@pytest.fixture(scope="module")
def release(tmp_path_factory):
    return run(tmp_path_factory.mktemp("rel"))


def test_headers_equal_the_module(release):
    out, _ = release
    for name, cols in ((psc.MODELS_FILE, psc.MODEL_COLUMNS), (psc.RESIDUES_FILE, psc.RESIDUE_COLUMNS),
                       (psc.FEATURES_FILE, psc.FEATURE_COLUMNS), (psc.VARIANTS_FILE, psc.VARIANT_COLUMNS),
                       (psc.SCORES_FILE, psc.SCORE_COLUMNS)):
        header, df = read_tsv(out / name)
        assert header == psc.header_line(cols)
        assert df.shape[1] == len(cols)


def frame(out, name, cols):
    _, df = read_tsv(out / name)
    df.columns = [c for c, _ in cols]
    for c, t in cols:
        if t is int and (df[c] != "").all():
            df[c] = pd.to_numeric(df[c])
    return df


def test_models(release):
    out, _ = release
    m = frame(out, psc.MODELS_FILE, psc.MODEL_COLUMNS).set_index("model_file_id")
    assert set(m.index) == {"AF-Q00001-F1-model_v6", "AF-Q00001-F2-model_v6", "AF3-Q00002-F1-model_v6",
                            "AF-Q00003-F1-model_v6", "AF-Q00006-F1-model_v6", "AF-Q00007-F1-model_v6"} | \
        {f"AF-Q00010-F{k}-model_v6" for k in range(1, 12)}
    f1, f2 = m.loc["AF-Q00001-F1-model_v6"], m.loc["AF-Q00001-F2-model_v6"]
    assert (f1.owned_start, f1.owned_end, f2.owned_start, f2.owned_end) == (1, 4, 5, 10)
    assert m.loc["AF3-Q00002-F1-model_v6", "licence"] == psc.AF3_LICENCE
    assert m.loc["AF-Q00003-F1-model_v6", "licence"] == "CC-BY-4.0"
    # accession under two genes: both, sorted by id, names in the same order
    assert m.loc["AF3-Q00002-F1-model_v6", ["gene_id", "gene_name"]].tolist() == \
        ["ENSG00000000002,ENSG00000000003", "GENEB,GENEC"]
    assert m.loc["AF-Q00003-F1-model_v6", "gene_id"] == ""  # its gene is not in GENCODE
    q10 = m[m.accession == "Q00010"]
    assert sorted(q10.window_index.astype(int)) == list(range(1, 12))
    s = json.loads((out.parent / "_work" / "models.json").read_text())
    assert (s["models_gene_not_in_gencode"], s["models_under_several_genes"]) == (1, 1)
    assert s["models_residue_differs_from_sequence"] == {"AF-Q00007-F1-model_v6": ["3:C>U"]}
    assert m.loc["AF-Q00003-F1-model_v6", ["has_pae", "pae_path", "pae_sha256", "pae_bytes"]].tolist() == \
        ["false", "", "", ""]
    for mid, r in m.iterrows():
        p = out / r.structure_path
        assert hashlib.sha256(p.read_bytes()).hexdigest() == r.structure_sha256 == r.source_pdb_sha256


def test_pae_objects_decode_to_the_symmetrised_source(release):
    out, paes = release
    m = frame(out, psc.MODELS_FILE, psc.MODEL_COLUMNS).set_index("model_file_id")
    for mid, src in paes.items():
        n = src.shape[0]
        data = (out / psc.pae_path(mid)).read_bytes()
        assert len(data) == n * n == int(m.loc[mid, "pae_bytes"])
        got = np.frombuffer(data, np.uint8).reshape(n, n)
        sym = np.minimum(src, src.T)
        assert all(got[i, j] == psc.encode_pae(float(sym[i, j])) for i in range(n) for j in range(n))


def test_pae_parse_takes_the_matrix_not_the_header():
    m = [[0.1, 31.7], [2.5, 0.0]]
    raw = json.dumps([{"max_predicted_aligned_error": 31.75, "predicted_aligned_error": m}]).encode()
    assert np.array_equal(mps.parse_pae_json(raw), np.array(m))
    with pytest.raises(ValueError):
        mps.parse_pae_json(b'[{"predicted_aligned_error":[[1,2],[3]]}]')


def test_encode_matrix_equals_module():
    vals = np.arange(0, 400) / 10
    assert mps.encode_matrix(vals.reshape(20, 20)).ravel().tolist() == [psc.encode_pae(float(v)) for v in vals]


def test_partial_pae_refused(tmp_path):
    with pytest.raises(ValueError, match="PAE for only some"):
        run(tmp_path, drop_pae="AF-Q00001-F2-predicted_aligned_error_v6.json.gz")


def test_residues_owner_unique_and_sorted(release):
    out, _ = release
    r = frame(out, psc.RESIDUES_FILE, psc.RESIDUE_COLUMNS)
    owners = r[r.is_owner == "true"].groupby(["gene_id", "accession", "protein_position"]).size()
    allpos = r.groupby(["gene_id", "accession", "protein_position"]).size()
    assert (owners == 1).all() and len(owners) == len(allpos)
    assert r[list(psc.RESIDUE_SORT)].equals(r.sort_values(list(psc.RESIDUE_SORT))[list(psc.RESIDUE_SORT)].reset_index(drop=True))
    # Q00002 under two genes is repeated; Q00003's gene is not in GENCODE so it has no rows
    assert set(r[r.accession == "Q00002"].gene_id) == {"ENSG00000000002", "ENSG00000000003"}
    assert "Q00003" not in set(r.accession)
    q1 = r[r.accession == "Q00001"]
    assert len(q1) == 6 + 8
    q10 = r[r.accession == "Q00010"]
    assert q10.model_file_id.iloc[30] == "AF-Q00010-F10-model_v6"  # string order, not window order
    assert "".join(q1[q1.is_owner == "true"].sort_values("protein_position").residue) == SEQ["Q00001"]
    row = q1[(q1.protein_position == 5) & (q1.is_owner == "true")].iloc[0]
    assert (row.model_file_id, row.file_position, row.mu, row.n_possible) == ("AF-Q00001-F2-model_v6", 3, "2e-08", "8")
    assert q1[q1.protein_position == 6].mu.tolist() == ["", ""]


def test_variants(release):
    out, _ = release
    v = frame(out, psc.VARIANTS_FILE, psc.VARIANT_COLUMNS)
    assert v[list(psc.VARIANT_SORT)].equals(v.sort_values(list(psc.VARIANT_SORT))[list(psc.VARIANT_SORT)].reset_index(drop=True))
    assert "Q00004" not in set(v.accession)
    b = v[v.pos == 1000].iloc[0]
    assert (b.model_file_id, b.file_position, b.variant, b.gene_id, b.mu_snp) == \
        ("AF-Q00001-F2-model_v6", 3, "1:1000:A:G", "ENSG00000000001", "1.5e-09")
    two = v[v.variant == "1:2000:G:A"]
    assert sorted(two.accession) == ["Q00001", "Q00002"]
    # Q00002's two genes both contain 2000: ambiguous, so empty
    # GENEA does not reach 2000, and Q00002's two genes both do: empty either way
    assert two.gene_id.tolist() == ["", ""]
    s = json.loads((out.parent / "_work" / "variants.json").read_text())
    assert s["drops"]["rows_accession_without_model"] == 1
    assert s["accessions_dropped_without_model"] == ["Q00004"]
    # 1000 (aa 5) and 2000 (aa 4) both sit in the F1/F2 overlap on the side the map did not pick
    assert s["rows_address_differs_from_map"] == 2
    assert s["variants_with_several_accessions"] == 1


def test_features(release):
    out, _ = release
    f = frame(out, psc.FEATURES_FILE, psc.FEATURE_COLUMNS)
    assert f[list(psc.FEATURE_SORT)].equals(f.sort_values(list(psc.FEATURE_SORT))[list(psc.FEATURE_SORT)].reset_index(drop=True))
    u = f[f.source == "uniprot"]
    assert set(u.accession) == {"Q00001"}
    bs = u[u.feature_type == "Binding site"]
    assert bs.description.tolist() == ["calcium", "zinc"]
    dom = u[u.feature_type == "Domain"].iloc[0]
    assert (dom.evidence, dom.evidence_tier) == ("ECO:0000255,ECO:0000269", "experimental")
    assert u[u.description == "odd code"].evidence_tier.tolist() == ["unknown"]
    t = f[f.source == "ted"]
    assert sorted(zip(t.start, t.end))[:2] == [(1, 3), (7, 9)] and set(t.evidence_tier) == {"predicted"}
    assert not set(u.accession) & {"Q00006", "Q00007"}
    assert sorted(zip(u[u.start == 3].feature_type, u[u.start == 3].description))[:2] == [("Region", "a"), ("Site", "b-desc")]
    assert list(zip(t.accession, t.start, t.end))[-1] == ("Q00010", 205, 210)
    s = json.loads((out.parent / "_work" / "features.json").read_text())
    assert s["uniprot"]["eco_codes_unknown_tier"] == ["ECO:0000999"]
    assert s["drops"]["uniprot_rows_position_unknown"] == 1
    assert s["drops"]["uniprot_rows_position_outside_sequence"] == 1
    assert s["drops"]["uniprot_rows_bundle_sequence_not_canonical_isoform"] == 3
    assert s["drops"]["ted_domains_sequence_md5_differs"] == 1


def test_scores_join_mean_and_percentile(release):
    out, _ = release
    s = frame(out, psc.SCORES_FILE, psc.SCORE_COLUMNS)
    assert set(s["predictor"]) == set(psc.SCORE_DIRECTIONS)
    am, esm = (s[s["predictor"] == p] for p in ("alphamissense", "esm1b"))
    # one row per protein position however many windows hold it, repeated under each gene of Q00002
    assert am.groupby(["gene_name", "accession"]).size().to_dict() == {
        ("GENEA", "Q00001"): 10, ("GENEB", "Q00002"): 5, ("GENEC", "Q00002"): 5, ("GENEG", "Q00007"): 4}
    assert esm.groupby(["gene_name", "accession"]).size().to_dict() == {
        ("GENEA", "Q00001"): 10, ("GENEB", "Q00002"): 5, ("GENEC", "Q00002"): 5}
    assert 3 not in set(am[am["accession"] == "Q00007"]["protein_position"].astype(int))  # unscored: no row
    assert "".join(am[am["gene_name"] == "GENEA"]["aa_ref"]) == SEQ["Q00001"]

    def row(df, acc, pos, gene=None):
        hit = df[(df["accession"] == acc) & (df["protein_position"] == pos)]
        return (hit[hit["gene_name"] == gene] if gene else hit).iloc[0]
    for df, score in ((am, am_score), (esm, esm_score)):
        r = row(df, "Q00001", 5)
        others = [score("Q00001", 5, alt) for alt in mps.AA if alt != "Y"]
        assert (float(r["mean_snv"]), r["n_snv"], r["n_all"]) == (score("Q00001", 5, "C"), 1, 19)
        assert float(r["mean_all"]) == pytest.approx(sum(others) / 19, abs=1e-5)
        # the map's row at position 4 names T where the sequence has A, so it places nothing
        r = row(df, "Q00001", 4)
        assert (r["mean_snv"], r["n_snv"], r["percentile"], r["n_all"]) == ("", 0, "", 19)
    # ESM1b's Q00002 is the isoform file, not the canonical one of another sequence
    assert float(row(esm, "Q00002", 5, "GENEB")["mean_all"]) == pytest.approx(
        sum(esm_score("Q00002", 5, alt) for alt in mps.AA if alt != "P") / 19, abs=1e-5)

    # ties share a mid-rank, and Q00002's position counts once though it is written under two genes
    def pct(df):
        return {(r.gene_name, r.accession, int(r.protein_position)): float(r.percentile)
                for r in df[df["percentile"] != ""].itertuples()}
    low, high = psc.score_percentiles([0.2, 0.2, 0.9])[1:]
    assert (low, high) == (33.333, 83.333)
    assert pct(am) == {("GENEA", "Q00001", 5): low, ("GENEB", "Q00002", 2): low, ("GENEC", "Q00002", 2): low,
                       ("GENEA", "Q00001", 9): high}
    assert pct(esm) == {("GENEA", "Q00001", 5): low, ("GENEA", "Q00001", 9): low, ("GENEB", "Q00002", 2): high,
                        ("GENEC", "Q00002", 2): high}
    values = np.random.default_rng(3).integers(0, 50, 500) / 10
    assert np.round(mps.midrank_percentiles(values), psc.SCORE_PERCENTILE_DECIMALS).tolist() == \
        psc.score_percentiles(values.tolist())

    key = s.assign(**{c: s[c].astype(int) for c in ("chrom", "gene_start", "protein_position")})
    assert key.sort_values(list(psc.SCORE_SORT), kind="stable").index.tolist() == s.index.tolist()
    got = subprocess.run(["tabix", str(out / psc.SCORES_FILE), "1:1900-1900"], capture_output=True, text=True,
                         check=True).stdout
    assert got.count("\n") == 10


def test_scores_reference_check_and_manifest(release):
    out, _ = release
    man = json.loads((out / psc.MANIFEST_FILE).read_text())
    assert man["counts"]["scores"]["variant_map"] == {
        "variant_rows": 4, "variant_rows_not_on_the_sequence": 1, "snv_substitutions": 3}
    assert man["counts"]["scores"]["score_rows_by_predictor"] == {"alphamissense": 24, "esm1b": 20}
    am, esm = man["score_sources"]["alphamissense"], man["score_sources"]["esm1b"]
    assert am["counts"] == {"proteins": 3, "residues": 19, "residues_with_mean_snv": 3,
                            "proteins_kept_without_gencode_gene": 1, "residues_unscored_in_kept_proteins": 1,
                            "source_proteins_not_in_release": 2}
    assert esm["counts"] == {"proteins": 2, "residues": 15, "residues_with_mean_snv": 3,
                             "proteins_kept_without_gencode_gene": 0, "residues_unscored_in_kept_proteins": 0,
                             "source_proteins_not_in_release": 0, "source_files_not_used": 2}
    drops = {k.split(".", 1)[1]: v for k, v in man["drops"].items() if k.startswith("scores.")}
    assert drops == {
        "alphamissense.proteins_absent_from_source": 0,
        "alphamissense.proteins_substitution_scored_twice": 0,
        "alphamissense.proteins_source_longer_than_release_sequence": 1,
        "alphamissense.proteins_reference_residue_differs": 1,
        "alphamissense.proteins_last_residue_unscored": 0,
        "esm1b.proteins_absent_from_source": 2,
        "esm1b.proteins_substitution_scored_twice": 0,
        "esm1b.proteins_source_longer_than_release_sequence": 0,
        "esm1b.proteins_reference_residue_differs": 1,
        "esm1b.proteins_last_residue_unscored": 1,
    }
    for name, src in man["score_sources"].items():
        assert src["direction"] == psc.SCORE_DIRECTIONS[name]
        assert {"source", "url", "version", "licence", "licence_statement", "summary"} <= set(src)
        for f, o in src["files"].items():
            assert o[psc.OBJECT_SHA256] == hashlib.sha256(
                (out.parent.parent / "predictors" / name / f).read_bytes()).hexdigest()
    assert "reader" not in am


def test_failed_scores_gate_leaves_the_earlier_table_and_no_stage_json(fresh, tmp_path, monkeypatch):
    work = fresh.parent / "_work"
    assert (work / "scores.json").exists()
    before = (fresh / psc.SCORES_FILE).read_bytes()
    monkeypatch.setattr(mps, "MAX_GENE_SCORE_BYTES", 10)
    with pytest.raises(AssertionError, match="score rows"):
        mps.main(args_for(tmp_path, "scores"))
    assert not (work / "scores.json").exists()
    assert (fresh / psc.SCORES_FILE).read_bytes() == before
    with pytest.raises(RuntimeError, match="scores"):
        mps.main(args_for(tmp_path, "manifest"))


def test_scores_refuse_an_unpinned_file_and_drop_a_protein_scored_twice(tmp_path, monkeypatch):
    out, _ = run(tmp_path, stages="pae,models,residues,features,variants")
    assert not (out / psc.SCORES_FILE).exists()
    # a second reference residue at one position of Q00001, as the real file has for some proteins
    twice = write_alphamissense(tmp_path, extra=[f"Q00001\tT1{alt}\t0.5\tambiguous" for alt in mps.AA if alt != "T"])
    with pytest.raises(ValueError, match="pinned"):
        mps.main(args_for(tmp_path, "scores"))
    assert not (out / psc.SCORES_FILE).exists() and not (out / (psc.SCORES_FILE + ".tbi")).exists()
    assert not (out.parent / "_work" / "scores.json").exists()
    monkeypatch.setitem(mps.SCORE_SOURCES["alphamissense"], "files", twice)
    mps.main(args_for(tmp_path, "scores"))
    stats = json.loads((out.parent / "_work" / "scores.json").read_text())
    assert stats["drops"]["alphamissense.proteins_substitution_scored_twice"] == 1
    assert stats["drops"]["alphamissense.proteins_reference_residue_differs"] == 1  # still only Q00006
    s = frame(out, psc.SCORES_FILE, psc.SCORE_COLUMNS)
    assert set(s[s["accession"] == "Q00001"]["predictor"]) == {"esm1b"}


def score_index(seqs):
    accs = list(seqs)
    return mps.ScoreIndex(types.SimpleNamespace(
        models=pd.DataFrame({"uniprot_id": accs}), sequence=seqs, isoform={a: a for a in accs},
        pairs=pd.DataFrame({"accession": accs})))


def score_all(table, acc, seq, score):
    """every substitution of every position of `seq`, scored by score(position, alt)."""
    subs = [(pos, ref, alt) for pos, ref in enumerate(seq, 1) for alt in mps.AA if alt != ref]
    table.add([acc], np.zeros(len(subs), dtype=np.int64), np.array([p for p, _, _ in subs]),
              np.array([ord(r) for _, r, _ in subs], dtype=np.uint8),
              np.array([mps.AA.index(a) for _, _, a in subs]), np.array([score(p, a) for p, _, a in subs]))
    return table.close()


def test_mean_snv_counts_a_substitution_once_however_many_variants_give_it(tmp_path):
    rel = score_index({"P1": "MKT"})
    # two variants give K2R, one gives K2E
    rows = [(1, 10, "A", "G", "1:10:A:G", "P1", "G1", "GENE", 2, "K", "R", "AF-P1-F1-model_v6", 2, ""),
            (1, 11, "A", "C", "1:11:A:C", "P1", "G1", "GENE", 2, "K", "R", "AF-P1-F1-model_v6", 2, ""),
            (1, 12, "A", "G", "1:12:A:G", "P1", "G1", "GENE", 2, "K", "E", "AF-P1-F1-model_v6", 2, "")]
    path = tmp_path / psc.VARIANTS_FILE
    path.write_bytes(gzip.compress((psc.header_line(psc.VARIANT_COLUMNS)
                                    + "".join("\t".join(map(str, r)) + "\n" for r in rows)).encode()))
    snv, stats = mps.snv_substitutions(path, rel)
    assert stats == {"variant_rows": 3, "variant_rows_not_on_the_sequence": 0, "snv_substitutions": 2}
    _, _, mean_snv, n_snv, *_ = score_all(mps.ScoreTable(rel, snv), "P1", "MKT",
                                          lambda p, a: {(2, "R"): 0.9, (2, "E"): 0.3}.get((p, a), 0.5))
    assert (mean_snv[1], n_snv[1]) == (0.6, 2)  # over the three variant rows it would be 0.7
    assert n_snv.tolist() == [0, 2, 0]


def test_percentile_ranks_the_mean_as_written(tmp_path):
    rel = score_index({"P1": "MKTA"})
    snv = np.zeros((rel.n, len(mps.AA)), dtype=bool)
    a, c = mps.AA.index("A"), mps.AA.index("C")
    snv[0, [a, c]] = snv[1, [a, c]] = snv[2, a] = True
    # positions 1 and 2 differ only past the decimals written, so they are one value and one rank
    scores = {(1, "A"): 0.1, (1, "C"): 0.100001, (2, "A"): 0.1, (2, "C"): 0.1, (3, "A"): 0.5}
    _, _, mean_snv, _, percentile, *_ = score_all(mps.ScoreTable(rel, snv), "P1", "MKTA",
                                                  lambda p, alt: scores.get((p, alt), 0.7))
    assert mean_snv[:3].tolist() == [0.1, 0.1, 0.5]
    assert np.round(percentile[:3], psc.SCORE_PERCENTILE_DECIMALS).tolist() == psc.score_percentiles([0.1, 0.1, 0.5])
    assert percentile[0] == percentile[1] and np.isnan(percentile[3])


def test_written_percentiles_reach_the_ends_of_the_scale():
    # the mid-rank is strictly inside (0, 100); three decimals are not
    p = psc.score_percentiles(list(range(100_001)))
    assert (p[0], p[-1]) == (0.0, 100.0) and 0 < 100 * 0.5 / 100_001


def test_manifest_refuses_scores_from_another_variant_map(fresh, tmp_path):
    man = json.loads((fresh / psc.MANIFEST_FILE).read_text())
    assert mps.SCORES_VARIANTS_KEY not in man["counts"]["scores"]
    path = fresh.parent / "_work" / "scores.json"
    stats = json.loads(path.read_text())
    assert stats[mps.SCORES_VARIANTS_KEY] == man[psc.MANIFEST_OBJECTS][psc.VARIANTS_FILE][psc.OBJECT_SHA256]
    mps.main(args_for(tmp_path, "manifest"))  # up to date: rebuilt without complaint
    for change in ({mps.SCORES_VARIANTS_KEY: "0" * 64}, None):
        edited = {**stats, **change} if change else {k: v for k, v in stats.items() if k != mps.SCORES_VARIANTS_KEY}
        path.write_text(json.dumps(edited))
        with pytest.raises(RuntimeError, match="rerun --stages scores,manifest"):
            mps.main(args_for(tmp_path, "manifest"))


def test_scores_preconditions_are_checked_before_any_stage(tmp_path):
    make_bundle(tmp_path)
    (tmp_path / "ted_meta.json").write_text("{}")
    a = args_for(tmp_path)
    i = a.index("--predictors")
    with pytest.raises(ValueError, match="--predictors"):
        mps.main(a[:i] + a[i + 2:])
    (tmp_path / "predictors" / "esm1b" / "ALL_hum_isoforms_ESM1b_LLR.zip").unlink()
    with pytest.raises(ValueError, match="ALL_hum_isoforms_ESM1b_LLR.zip"):
        mps.main(a)
    out = tmp_path / "build" / REL
    assert not (out / psc.PAE_DIR).exists() and not (out / psc.MODELS_FILE).exists()
    mps.main(args_for(tmp_path, "pae,models"))  # stages that do not need the predictors still run


def test_scores_need_the_variant_stage(tmp_path):
    run(tmp_path, stages="pae,models,residues,features")
    with pytest.raises(RuntimeError, match="variants"):
        mps.main(args_for(tmp_path, "scores"))


def test_tabix_answers(release):
    out, _ = release
    got = subprocess.run(["tabix", str(out / psc.VARIANTS_FILE), "1:1000-1000"], capture_output=True, text=True,
                         check=True).stdout
    assert got.count("\n") == 1 and "1:1000:A:G" in got
    got = subprocess.run(["tabix", str(out / psc.RESIDUES_FILE), "1:900-900"], capture_output=True, text=True,
                         check=True).stdout
    assert got.count("\n") == 14


def test_manifest(release):
    out, _ = release
    man = json.loads((out / psc.MANIFEST_FILE).read_text())
    assert man["contract_version"] == psc.CONTRACT_VERSION
    assert man["no_model_accessions"] == ["Q00004"]
    assert man["genes_absent_from_gencode"] == ["ENSG00000000099"]
    assert set(man["lab_inputs"]) == set(psc.LAB_INPUTS)
    assert sorted(man[psc.MANIFEST_OBJECTS]) == sorted(psc.manifest_object_names())
    for name, o in man[psc.MANIFEST_OBJECTS].items():
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == o[psc.OBJECT_SHA256]
        assert (out / name).stat().st_size == o[psc.OBJECT_BYTES]
    assert man[psc.MANIFEST_DIRECTORIES][psc.STRUCTURE_DIR][psc.DIRECTORY_OBJECTS] == 17


class FakeBucket:
    """stands in for gcloud: objects by full gs:// uri, --no-clobber semantics, a call log."""

    def __init__(self, monkeypatch):
        self.objects, self.calls = {}, []
        monkeypatch.setattr(mps, "gs_list", self.ls)
        monkeypatch.setattr(mps, "gs_cp", self.cp)
        monkeypatch.setattr(mps, "gs_cat", lambda uri: self.objects[uri])
        monkeypatch.setattr(mps, "gs_read_if_exists", lambda uri: self.objects.get(uri))
        monkeypatch.setattr(mps, "gs_replace", self.replace)

    def ls(self, prefix):
        return {k[len(prefix) + 1:]: (len(v), mps.bytes_crc32c(v)) for k, v in self.objects.items()
                if k.startswith(prefix + "/")}

    def put(self, key, data):
        self.objects.setdefault(key, data)

    def cp(self, args, dest):
        self.calls.append(("cp", tuple(args), dest))
        for src in (a for a in args if a != "-r"):
            if src.startswith("gs://"):
                for k, v in list(self.objects.items()):
                    if k.startswith(src + "/"):
                        self.put(dest + src.rsplit("/", 1)[1] + k[len(src):], v)
            elif Path(src).is_dir():
                for f in Path(src).iterdir():
                    self.put(f"{dest}{Path(src).name}/{f.name}", f.read_bytes())
            else:
                self.put(dest + Path(src).name if dest.endswith("/") else dest, Path(src).read_bytes())

    def replace(self, local, uri):
        self.calls.append(("replace", uri))
        self.objects[uri] = Path(local).read_bytes()


REL = "2026-01-01-000000"
DEST = f"{mps.REHEARSAL_ROOT}/releases/{REL}"
PTR = f"{mps.REHEARSAL_ROOT}/{psc.POINTER}"


@pytest.fixture
def fresh(tmp_path):
    return run(tmp_path)[0]


def test_upload_and_pointer(fresh, monkeypatch):
    b = FakeBucket(monkeypatch)
    mps.upload(fresh, REL, mps.REHEARSAL_ROOT)
    assert b.calls[-1] == ("replace", PTR)
    ptr = json.loads(b.objects[PTR])
    man_sha = hashlib.sha256((fresh / psc.MANIFEST_FILE).read_bytes()).hexdigest()
    assert ptr == json.loads(psc.pointer(REL, DEST, ptr["written_at"], man_sha))
    assert hashlib.sha256(b.objects[f"{DEST}/{psc.MANIFEST_FILE}"]).hexdigest() == man_sha
    # a resumed upload of the same build is accepted
    mps.upload(fresh, REL, mps.REHEARSAL_ROOT)


def test_upload_refuses_local_file_not_matching_manifest(fresh, monkeypatch):
    b = FakeBucket(monkeypatch)
    (fresh / psc.VARIANTS_FILE).write_bytes(b"garbage")
    with pytest.raises(mps.RefuseUpload, match="variants.tsv.gz"):
        mps.upload(fresh, REL, mps.REHEARSAL_ROOT)
    assert b.objects == {}


def test_upload_refuses_differing_remote_object(fresh, monkeypatch):
    b = FakeBucket(monkeypatch)
    b.objects[f"{DEST}/{psc.RESIDUES_FILE}"] = b"older bytes"
    with pytest.raises(mps.RefuseUpload, match="residues.tsv.gz"):
        mps.upload(fresh, REL, mps.REHEARSAL_ROOT)
    assert len(b.objects) == 1 and PTR not in b.objects
    b.objects = {f"{DEST}/stray.txt": b"x"}
    with pytest.raises(mps.RefuseUpload, match="not named"):
        mps.upload(fresh, REL, mps.REHEARSAL_ROOT)


def test_rebuilt_manifest_after_upload_is_refused(fresh, monkeypatch):
    b = FakeBucket(monkeypatch)
    mps.upload(fresh, REL, mps.REHEARSAL_ROOT)
    before = b.objects[PTR]
    man = json.loads((fresh / psc.MANIFEST_FILE).read_text())
    man["built_at"] = "2099-01-01T00:00:00+00:00"
    (fresh / psc.MANIFEST_FILE).write_text(json.dumps(man))
    with pytest.raises(mps.RefuseUpload, match="new release id"):
        mps.upload(fresh, REL, mps.REHEARSAL_ROOT)
    assert b.objects[PTR] == before


def test_pointer_needs_the_manifest_read_back(fresh, monkeypatch):
    b = FakeBucket(monkeypatch)
    monkeypatch.setattr(mps, "gs_cat", lambda uri: b"not the manifest")
    with pytest.raises(mps.RefuseUpload, match="read back"):
        mps.upload(fresh, REL, mps.REHEARSAL_ROOT)
    assert PTR not in b.objects


def test_copy_from_an_earlier_release_is_checked(fresh, monkeypatch):
    b = FakeBucket(monkeypatch)
    mps.upload(fresh, REL, mps.REHEARSAL_ROOT)
    new = f"{mps.REHEARSAL_ROOT}/releases/2026-01-02-000000"
    mps.upload(fresh, "2026-01-02-000000", mps.REHEARSAL_ROOT, copy_from=DEST, move_pointer=True)
    assert b.objects[f"{new}/{psc.structure_path('AF-Q00001-F1-model_v6')}"] == \
        (fresh / psc.structure_path("AF-Q00001-F1-model_v6")).read_bytes()
    assert any(c[0] == "cp" and f"{DEST}/{psc.STRUCTURE_DIR}" in c[1] for c in b.calls)
    b.objects[f"{DEST}/{psc.pae_path('AF-Q00001-F1-model_v6')}"] = b"corrupt"
    with pytest.raises(mps.RefuseUpload, match="pae/AF-Q00001-F1"):
        mps.upload(fresh, "2026-01-03-000000", mps.REHEARSAL_ROOT, copy_from=DEST, move_pointer=True)
    # refused before copying: the new release directory stays empty
    assert not any("/2026-01-03-000000/" in k for k in b.objects)


def test_root_guards(fresh, tmp_path):
    for bad in ("gs://daly-genetics-results/elsewhere", mps.PRODUCTION_ROOT + "_x", mps.PRODUCTION_ROOT + "/sub"):
        with pytest.raises(ValueError):
            mps.check_root(bad)
        with pytest.raises(ValueError):
            mps.upload(fresh, REL, bad)
    with pytest.raises(ValueError):
        mps.main(args_for(tmp_path, "manifest", ("--root", mps.PRODUCTION_ROOT + "_x", "--workdir", str(tmp_path / "w2"))))
    assert not (tmp_path / "w2").exists()


def test_models_rebuild_is_byte_identical(fresh, tmp_path):
    before = (fresh / psc.MODELS_FILE).read_bytes()
    mps.main(args_for(tmp_path, "models"))
    assert (fresh / psc.MODELS_FILE).read_bytes() == before


def test_failed_gate_leaves_no_output(tmp_path, monkeypatch):
    out, _ = run(tmp_path, stages="pae,models")
    monkeypatch.setattr(mps, "MAX_GENE_RESIDUE_BYTES", 10)
    with pytest.raises(AssertionError, match="8 MiB"):
        mps.main(args_for(tmp_path, "residues"))
    assert not (out / psc.RESIDUES_FILE).exists() and not (out / (psc.RESIDUES_FILE + ".tbi")).exists()
    assert not (out.parent / "_work" / "residues.json").exists()
    with pytest.raises(RuntimeError, match="residues"):
        mps.main(args_for(tmp_path, "manifest"))


def test_pae_log_torn_line_and_deleted_object(fresh, tmp_path):
    log = fresh.parent / "_work" / "pae.tsv"
    lines = log.read_text().splitlines()
    log.write_text("\n".join(lines + [lines[1], "AF3-Q00002-F1-model_v6\t25\t0123456789abcdef012345678"]))
    (fresh / psc.pae_path("AF-Q00001-F2-model_v6")).unlink()
    mps.main(args_for(tmp_path, "pae,models"))
    rows = [ln.split("\t") for ln in log.read_text().splitlines()[1:]]
    assert len(rows) == len({r[0] for r in rows}) == 3
    m = frame(fresh, psc.MODELS_FILE, psc.MODEL_COLUMNS)
    assert m[m.has_pae == "true"].source_pae_sha256.str.len().eq(64).all()


def test_model_row_checks():
    good = "\t".join(["AF-X-F1-model_v6", "X", "X", "", "", "S", "AF2DB", "L", "v", "CC", "1", "1", "1", "2", "2",
                      "1", "2", "true", psc.structure_path("AF-X-F1-model_v6"), "a" * 64, "9",
                      psc.pae_path("AF-X-F1-model_v6"), "b" * 64, "4", "a" * 64, "c" * 64]) + "\n"
    mps.check_model_rows(good)
    for bad in (good.replace("c" * 64, "c" * 25), good.replace("\t4\t", "\t5\t"), good.replace("\n", "\textra\n"),
                good.replace("\t1\t1\t1\t2", "\t2\t1\t1\t2")):
        with pytest.raises(ValueError):
            mps.check_model_rows(bad)


def test_pae_size_checked_against_guide(fresh, tmp_path):
    src = tmp_path / "bundle" / "pae_files" / "AF-Q00001-F1-predicted_aligned_error_v6.json.gz"
    with pytest.raises(ValueError, match="guide says"):
        mps.convert_pae((src, tmp_path / "x.u8", 7))


def test_changed_input_is_rehashed(fresh, tmp_path):
    g = tmp_path / "gencode.tsv"
    g.write_text(g.read_text() + "ENSG00000000004.1\t1\t4900\t5100\t+\tGENED\tprotein_coding\n")
    mps.main(args_for(tmp_path, "manifest"))
    man = json.loads((fresh / psc.MANIFEST_FILE).read_text())
    assert man["gencode"]["sha256"] == hashlib.sha256(g.read_bytes()).hexdigest()


NEW = f"{mps.REHEARSAL_ROOT}/releases/2026-01-02-000000"


def test_copy_source_must_be_another_release_under_the_root(fresh, monkeypatch, tmp_path):
    b = FakeBucket(monkeypatch)
    mps.upload(fresh, REL, mps.REHEARSAL_ROOT)
    for src in ("gs://other-bucket/protein_structure_staging/releases/2026-01-01-000000",
                f"{mps.PRODUCTION_ROOT}/releases/{REL}", str(fresh), DEST, f"{DEST}/structures",
                f"{mps.REHEARSAL_ROOT}/releases/../x"):
        with pytest.raises(ValueError):
            mps.upload(fresh, REL if src == DEST else "2026-01-02-000000", mps.REHEARSAL_ROOT,
                       copy_from=src, move_pointer=True)
    with pytest.raises(SystemExit):
        mps.main(args_for(tmp_path, "", ("--copy-from", DEST)))


@pytest.mark.parametrize("change, why", [
    ("drop", "missing from the copy source"),
    ("extra", "not named by this manifest"),
    ("differ", "differs from the local object"),
])
def test_copy_source_checked_before_copying(fresh, monkeypatch, change, why):
    b = FakeBucket(monkeypatch)
    mps.upload(fresh, REL, mps.REHEARSAL_ROOT)
    key = f"{DEST}/{psc.structure_path('AF-Q00001-F1-model_v6')}"
    if change == "drop":
        del b.objects[key]
    elif change == "extra":
        b.objects[f"{DEST}/{psc.PAE_DIR}/AF-ZZZ-F1-model_v6.u8"] = b"x"
    else:
        b.objects[key] = b"y" * len(b.objects[key])
    with pytest.raises(mps.RefuseUpload, match=why):
        mps.upload(fresh, "2026-01-02-000000", mps.REHEARSAL_ROOT, copy_from=DEST, move_pointer=True)
    assert not any(k.startswith(NEW) for k in b.objects)


def test_pointer_moves_only_with_the_flag(fresh, monkeypatch):
    b = FakeBucket(monkeypatch)
    mps.upload(fresh, REL, mps.REHEARSAL_ROOT)
    with pytest.raises(mps.RefuseUpload, match="--move-pointer"):
        mps.upload(fresh, "2026-01-02-000000", mps.REHEARSAL_ROOT)
    assert json.loads(b.objects[PTR])["release"] == REL and not any(k.startswith(NEW) for k in b.objects)
    mps.upload(fresh, "2026-01-02-000000", mps.REHEARSAL_ROOT, move_pointer=True)
    assert json.loads(b.objects[PTR])["release"] == "2026-01-02-000000"
    # rerunning the older release's upload does not roll the pointer back
    with pytest.raises(mps.RefuseUpload, match="--move-pointer"):
        mps.upload(fresh, REL, mps.REHEARSAL_ROOT)
    assert json.loads(b.objects[PTR])["release"] == "2026-01-02-000000"


def test_release_id_validated(fresh, tmp_path, monkeypatch):
    FakeBucket(monkeypatch)
    for bad in ("../x", "2026-01-01-000000/..", "a/b", "2026-1-1-0", "2026-01-01-000000\n",
                "\u0662\u0660\u0662\u0666-01-01-000000", "2026-01-01-00000\u0660"):
        with pytest.raises(ValueError):
            psc.check_release(bad)
        with pytest.raises(ValueError):
            mps.upload(fresh, bad, mps.REHEARSAL_ROOT)
    a = args_for(tmp_path, "manifest")
    a[a.index("--release") + 1] = "../escape"
    a[a.index("--workdir") + 1] = str(tmp_path / "w3")
    with pytest.raises(ValueError):
        mps.main(a)
    assert not (tmp_path / "w3").exists() and not (tmp_path / "escape").exists()


def test_corrupt_local_structure_refused(fresh, monkeypatch):
    b = FakeBucket(monkeypatch)
    f = fresh / psc.structure_path("AF-Q00001-F1-model_v6")
    f.write_bytes(b"z" * f.stat().st_size)
    with pytest.raises(mps.RefuseUpload, match="models.tsv.gz row"):
        mps.upload(fresh, REL, mps.REHEARSAL_ROOT)
    assert b.objects == {}


def test_directory_totals_must_match_rows(fresh, monkeypatch):
    FakeBucket(monkeypatch)
    man = json.loads((fresh / psc.MANIFEST_FILE).read_text())
    man[psc.MANIFEST_DIRECTORIES][psc.PAE_DIR][psc.DIRECTORY_OBJECTS] += 1
    (fresh / psc.MANIFEST_FILE).write_text(json.dumps(man))
    with pytest.raises(mps.RefuseUpload, match="totals"):
        mps.upload(fresh, REL, mps.REHEARSAL_ROOT)


def test_copy_that_leaves_an_object_out_is_refused(fresh, monkeypatch):
    b = FakeBucket(monkeypatch)
    real = b.cp

    def lossy(args, dest):
        real(args, dest)
        b.objects.pop(f"{DEST}/{psc.pae_path('AF-Q00001-F1-model_v6')}", None)
    monkeypatch.setattr(mps, "gs_cp", lossy)
    with pytest.raises(mps.RefuseUpload, match="missing"):
        mps.upload(fresh, REL, mps.REHEARSAL_ROOT)
    assert PTR not in b.objects


def test_failed_models_gate_leaves_no_models_json(fresh, tmp_path, monkeypatch):
    work = fresh.parent / "_work"
    assert (work / "models.json").exists()

    def refuse(text):
        raise ValueError("gate")
    monkeypatch.setattr(mps, "check_model_rows", refuse)
    before = (fresh / psc.MODELS_FILE).read_bytes()
    with pytest.raises(ValueError, match="gate"):
        mps.main(args_for(tmp_path, "models"))
    assert not (work / "models.json").exists()
    assert (fresh / psc.MODELS_FILE).read_bytes() == before  # the gated table was not replaced


def test_failed_rerun_removes_the_earlier_stage_json(fresh, tmp_path, monkeypatch):
    work = fresh.parent / "_work"
    assert (work / "residues.json").exists()
    monkeypatch.setattr(mps, "MAX_GENE_RESIDUE_BYTES", 10)
    with pytest.raises(AssertionError):
        mps.main(args_for(tmp_path, "residues"))
    assert not (work / "residues.json").exists()


def test_pae_log_row_with_short_hash_ignored(fresh):
    log = fresh.parent / "_work" / "pae.tsv"
    mid = "AF-Q00001-F1-model_v6"
    size = (fresh / psc.pae_path(mid)).stat().st_size
    good = mps.read_pae_log(log, fresh)[mid]
    assert len(good[1]) == 64
    log.write_text(mps.PAE_LOG_HEADER + f"{mid}\t{size}\t{good[1][:25]}\t{good[2]}\n")
    assert mid not in mps.read_pae_log(log, fresh)
    log.write_text(mps.PAE_LOG_HEADER + f"{mid}\t{size}\t{good[1]}\t{good[2][:63]}\n")
    assert mid not in mps.read_pae_log(log, fresh)


def test_same_size_input_change_is_rehashed(fresh, tmp_path):
    g = tmp_path / "gencode.tsv"
    old = g.read_text()
    new = old.replace("GENEA", "GENEZ")
    assert len(new) == len(old) and new != old
    st = g.stat()
    g.write_text(new)
    os.utime(g, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))
    mps.main(args_for(tmp_path, "manifest"))
    man = json.loads((fresh / psc.MANIFEST_FILE).read_text())
    assert man["gencode"]["sha256"] == hashlib.sha256(new.encode()).hexdigest()
