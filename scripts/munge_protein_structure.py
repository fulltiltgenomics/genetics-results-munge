"""Build a protein structure release from the lab's reference bundle, per docs/protein-structure.md.

usage: munge_protein_structure.py --bundle DIR --lab-files DIR --gencode TSV --workdir DIR
           [--bundle-tar TAR] [--ted TSV --ted-meta JSON] [--release ID]
           [--stages pae,models,residues,features,variants,manifest] [--upload] [--root gs://...]
           [--copy-from gs://.../releases/<earlier release>]

The column lists, sort orders, tabix arguments, PAE encoding and the window rules are the
contract module's (`protein_structure_contract`); this script only reads the bundle and applies
them. Every stage writes into `<workdir>/<release>/`. A stage named in `--stages` runs again from
scratch, except `pae`, which converts only the files whose log row and object are both intact;
stages not named are left alone, so a failure in `variants` does not redo the PAE conversion.
A stage writes its outputs under temporary names, runs its gates, indexes, renames, and writes
`_work/<stage>.json` last; that file is what says the stage finished, it is removed when the
stage starts, and the manifest stage refuses without it. `--release` is fixed at the first run
in `<workdir>/RELEASE`; a new release needs a new workdir.

Upload (`--upload`) verifies every local object against the manifest (and each structure and PAE
object against its models.tsv.gz row) before copying anything, refuses when the release
directory in the bucket already holds an object the manifest does not name or one whose size or
crc32c differs, verifies the bucket against the manifest again after copying, and writes the
pointer only once the manifest read back from the bucket has the sha256 going into it. A manifest
rebuilt after an upload started therefore differs from the bucket's and is refused: build it
under a new release id. `--copy-from` copies structures and PAE objects bucket-to-bucket from an
earlier release instead of uploading them, under the same checks.

Inputs, as read from their bytes:
- bundle (`sir-reference-data/`): `pdb_pae_file_pos_guide.tsv` (pdb_filename, pae_filename,
  pos_covered "[a, b]", uniprot_id), `protein_sequence_guide.tsv` (uniprot_id, uniprot_isoform,
  ensg unversioned, gene_name, AF2DB, uniprot_canonical, sequence), `pdb_files/*.pdb.gz` (AF2DB
  files carry the AFDB header, AF3 files start at ATOM), `pae_files/*.json.gz` (a one-element
  list of {predicted_aligned_error: n x n, max_predicted_aligned_error}), and
  `all_missense_variants_gr38.h5` (`chr{1..22,X}_{ref_alt,pos,filename,uniprot_id}`, needs
  hdf5plugin; rows not in position order).
- lab files (`LAB_INPUTS`): `mu_per_residue_gnomad.parquet` keyed (uniprot_id, aa_pos);
  `gnomad_mu_snp_by_variant.parquet` keyed (chrom "chr1".."chrX", pos, ref, alt);
  `uniprot_all_feats_evidence.tsv` with float positions ("84.0"), evidence as
  `ECO:...|source|id` items joined by `;`, and `canonical_isoform` being UniProt's canonical
  isoform, which is not always `-1`.
- GENCODE v49 genes TSV (gene_id versioned, numeric chrom with X=23), the file
  `scripts/cellxgene/publish.py` reads; gene ids are matched with the `.N` version stripped.
- TED domain summary (optional): the human (tax 9606) rows of
  `ted_365m.domain_summary.cath.globularity.taxid.tsv.gz`, header-less, columns per the Zenodo
  record; `--ted-meta` is a JSON object recorded verbatim in the manifest (url, version, licence,
  sha256 of the bulk file).

Decisions that are not in the contract, stated where they bind:
- UniProt features are positions on UniProt's canonical isoform, so they are kept only for an
  accession whose bundle sequence is that isoform (`uniprot_canonical` true and the bundle's
  isoform, a bare accession meaning `-1`, equal to the feature file's `canonical_isoform`).
  Rows with an unknown position or one outside the bundle sequence are dropped and counted.
- `evidence` keeps the ECO ids only (the `|PubMed|id` provenance is dropped), deduplicated in
  source order and comma-joined, as the contract's column names it.
- a TED domain is kept only when the md5 of its residues read from the bundle sequence equals
  TED's `md5_domain`, which is what proves TED's AFDB v4 numbering is the bundle's. One row per
  segment of a discontinuous domain, so `start`/`end` never span residues the domain lacks.
- `mu` and `mu_snp` are written with MU_FORMAT. The lab's per-SNV rates carry 5 significant
  digits and the per-residue sums at most 8 (measured on the files), so 10 loses nothing; the
  build refuses a value that does not survive the format to 1e-12.
- structures and PAE objects are hashed per object in models.tsv.gz; the manifest lists the
  top-level files with their hashes (models.tsv.gz among them) and totals for the two
  directories rather than 41,000 entries.

The output is written by its own bgzip pipe rather than `peak_utils`: residue and variant rows
are produced per chromosome and never sit in one frame, and the header must be the contract's
`header_line()` byte for byte.
"""

import argparse
import datetime as dt
import gzip
import hashlib
import json
import multiprocessing as mp
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute  # noqa: F401  pa.compute
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import protein_structure_contract as psc  # noqa: E402

STAGES = ("pae", "models", "residues", "features", "variants", "manifest")
REHEARSAL_ROOT = "gs://daly-genetics-results/protein_structure_staging"
PRODUCTION_ROOT = "gs://daly-genetics-results/protein_structure"
CHROMS = [str(c) for c in range(1, 23)] + ["X"]
MAX_GENE_RESIDUE_BYTES = 8 * 1024 * 1024
TED_COLUMNS = ["ted_id", "md5_domain", "consensus_level", "chopping", "nres_domain", "num_segments",
               "plddt", "num_helix_strand_turn", "num_helix", "num_strand", "num_helix_strand",
               "num_turn", "proteome_id", "cath_label", "cath_assignment_level",
               "cath_assignment_method", "packing_density", "norm_rg", "tax_common_name",
               "tax_scientific_name", "tax_lineage"]
TED_ID_RE = re.compile(r"^AF-([A-Z0-9]+)-F(\d+)-model_v4_TED(\d+)$")
# AFDB v4 cut proteins over 2,700 residues into 1,400-residue fragments starting every 200
AFDB_FRAGMENT_STEP = 200
ECO_RE = re.compile(r"ECO:\d{7}")
VARIANT_BLOCK_ROWS = 300_000
MU_FORMAT = "{:.10g}"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
STAGE_JSON = ("models", "residues", "features", "variants")


class RefuseUpload(RuntimeError):
    pass


def file_digests(path, chunk=1 << 22):
    """(bytes, sha256 hex, crc32c base64) in one pass; crc32c is what a bucket listing always carries."""
    import base64

    import google_crc32c
    h, c, n = hashlib.sha256(), google_crc32c.Checksum(), 0
    with open(path, "rb") as f:
        while b := f.read(chunk):
            h.update(b)
            c.update(b)
            n += len(b)
    return n, h.hexdigest(), base64.b64encode(c.digest()).decode()


def bytes_crc32c(data):
    import base64

    import google_crc32c
    return base64.b64encode(google_crc32c.Checksum(data).digest()).decode()


def check_mu_format(values, what):
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    back = np.array([float(MU_FORMAT.format(x)) for x in v])
    bad = np.abs(back - v) > np.abs(v) * 1e-12
    if bad.any():
        raise ValueError(f"{what}: {int(bad.sum())} values do not survive {MU_FORMAT}, e.g. {v[bad][0]!r}")


def sha256_file(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


def write_json(path, obj):
    tmp = Path(str(path) + ".part")
    tmp.write_text(json.dumps(obj, indent=1, sort_keys=True))
    tmp.rename(path)


# ---------------------------------------------------------------- inputs


def load_gencode(path):
    g = pd.read_csv(path, sep="\t", dtype=str)
    g["gene_id"] = g["gene_id"].str.replace(r"\.\d+", "", regex=True)
    g = g.drop_duplicates("gene_id", keep="first")
    for c in ("chrom", "gene_start", "gene_end"):
        g[c] = g[c].astype(int)
    return g.set_index("gene_id")[["chrom", "gene_start", "gene_end", "gene_name"]]


class Bundle:
    """the bundle's guides resolved into windows: which accessions are in the release and why not."""

    def __init__(self, bundle_dir, gencode):
        self.dir = Path(bundle_dir)
        self.drops = {}
        g = pd.read_csv(self.dir / "pdb_pae_file_pos_guide.tsv", sep="\t", dtype=str)
        ab = g["pos_covered"].str.extract(r"^\[(\d+),\s*(\d+)\]$")
        if ab.isna().any().any():
            raise ValueError("unparseable pos_covered in the file guide")
        g["protein_start"], g["protein_end"] = ab[0].astype(int), ab[1].astype(int)
        g["model_file_id"] = [psc.model_file_id(f) for f in g["pdb_filename"]]
        if g["model_file_id"].duplicated().any():
            raise ValueError("model file ids are not unique")
        g["pdb_present"] = [(self.dir / "pdb_files" / f).exists() for f in g["pdb_filename"]]
        g["pae_present"] = [isinstance(f, str) and (self.dir / "pae_files" / f).exists() for f in g["pae_filename"]]

        rows, no_model = [], []
        for acc, x in g.groupby("uniprot_id", sort=False):
            if not x["pdb_present"].any():
                no_model.append(acc)
                continue
            if not x["pdb_present"].all():
                raise ValueError(f"{acc}: some but not all model files are in the bundle")
            if x["pae_present"].any() and not x["pae_present"].all():
                raise ValueError(f"{acc}: PAE for only some of its model files")
            windows = list(zip(x["protein_start"], x["protein_end"]))
            # owning_window refuses windows out of start order, which is the lab code's assumption
            owned = psc.owning_window(windows)
            x = x.copy()
            x["window_index"] = range(1, len(x) + 1)
            x["n_windows"] = len(x)
            x["owned_start"] = [o[0] for o in owned]
            x["owned_end"] = [o[1] for o in owned]
            rows.append(x)
        self.models = pd.concat(rows, ignore_index=True)
        self.models["has_pae"] = self.models["pae_present"]
        self.no_model_accessions = sorted(no_model)
        self.drops["accessions_without_model_files"] = len(no_model)
        named = set(g["pae_filename"].dropna())
        self.drops["pae_files_named_by_no_guide_row"] = sum(
            1 for p in (self.dir / "pae_files").iterdir() if p.name not in named)

        s = pd.read_csv(self.dir / "protein_sequence_guide.tsv", sep="\t", dtype=str, keep_default_na=False)
        self.seq = s
        released = set(self.models["uniprot_id"])
        first = s.drop_duplicates("uniprot_id").set_index("uniprot_id")
        if (s.groupby("uniprot_id")["sequence"].nunique() > 1).any():
            raise ValueError("an accession has two sequences in the sequence guide")
        self.sequence = first["sequence"].to_dict()
        self.isoform = first["uniprot_isoform"].to_dict()
        self.canonical = (first["uniprot_canonical"] == "True").to_dict()
        self.af2db = (first["AF2DB"] == "True").to_dict()
        missing = released - set(first.index)
        if missing:
            raise ValueError(f"accessions with models but no sequence guide row: {sorted(missing)[:5]}")
        self.source_symbols = s.groupby("uniprot_id")["gene_name"].agg(lambda v: ",".join(dict.fromkeys(x for x in v if x))).to_dict()

        # distinct (gene, accession) pairs for gene-keyed rows; a gene absent from GENCODE gets none
        pairs = s[s["uniprot_id"].isin(released)][["ensg", "uniprot_id"]].drop_duplicates()
        pairs = pairs.rename(columns={"ensg": "gene_id", "uniprot_id": "accession"})
        inv49 = pairs["gene_id"].isin(gencode.index)
        self.genes_absent_from_gencode = sorted(set(pairs.loc[~inv49, "gene_id"]))
        self.drops["gene_accession_pairs_gene_not_in_gencode"] = int((~inv49).sum())
        pairs = pairs[inv49].join(gencode, on="gene_id")
        self.pairs = pairs.sort_values(["chrom", "gene_start", "gene_id", "accession"]).reset_index(drop=True)
        ngenes = self.pairs.groupby("accession")["gene_id"].nunique()
        self.accessions_in_several_genes = sorted(ngenes[ngenes > 1].index)
        self.drops["accessions_without_gencode_gene"] = len(released - set(self.pairs["accession"]))


# ---------------------------------------------------------------- PAE


def parse_pae_json(raw: bytes) -> np.ndarray:
    """the predicted_aligned_error matrix of an AFDB/AF3 PAE json, without a Python list per value.

    Only the matrix is parsed; `max_predicted_aligned_error` is a header field, not a PAE entry,
    and must never be encoded as one.
    """
    k = raw.index(b'"predicted_aligned_error"')
    a = raw.index(b"[[", k)
    b = raw.index(b"]]", a)
    body = raw[a + 2:b].translate(None, b" \t\r\n")
    n = body.count(b"],[") + 1
    vals = np.fromstring(body.replace(b"],[", b",").decode(), sep=",", dtype=np.float64)
    if vals.size != n * n:
        raise ValueError(f"PAE matrix is not square: {n} rows, {vals.size} values")
    return vals.reshape(n, n)


def encode_matrix(m: np.ndarray) -> np.ndarray:
    """psc.encode_pae over every entry, computed once per distinct value so the module stays the definition."""
    uniq, inv = np.unique(m, return_inverse=True)
    codes = np.array([psc.encode_pae(float(v)) for v in uniq], dtype=np.uint8)
    return codes[inv].reshape(m.shape)


def convert_pae(task):
    src, dst, n_expected = task
    raw_gz = Path(src).read_bytes()
    m = parse_pae_json(gzip.decompress(raw_gz))
    if m.shape[0] != n_expected:
        raise ValueError(f"{src}: PAE is {m.shape[0]} x {m.shape[0]}, the guide says {n_expected} residues")
    m = np.minimum(m, m.T)
    enc = np.ascontiguousarray(encode_matrix(m))
    data = enc.tobytes()
    part = Path(str(dst) + ".part")
    part.write_bytes(data)
    part.rename(dst)
    return Path(dst).stem, len(data), hashlib.sha256(data).hexdigest(), hashlib.sha256(raw_gz).hexdigest()


PAE_LOG_HEADER = "model_file_id\tpae_bytes\tpae_sha256\tsource_pae_sha256\n"


def read_pae_log(log, out):
    """{model_file_id: (bytes, sha256, source sha256)} for the log rows whose object is intact.

    The log is appended as files finish, so a killed run can leave a torn last line and a rerun
    after an object was deleted appends a second row for it: a row is taken only when every
    field is well formed and the object on disk has its size, and the last such row wins.
    """
    rows = {}
    if not log.exists():
        return rows
    for line in log.read_text().split("\n")[1:]:
        f = line.split("\t")
        if len(f) != 4 or not psc.MODEL_FILE_ID_RE.match(f[0]) or not f[1].isdigit() \
                or not SHA256_RE.match(f[2]) or not SHA256_RE.match(f[3]):
            continue
        obj = out / psc.pae_path(f[0])
        if obj.exists() and obj.stat().st_size == int(f[1]):
            rows[f[0]] = (int(f[1]), f[2], f[3])
    return rows


def stage_pae(bundle, out, work, procs):
    (out / psc.PAE_DIR).mkdir(parents=True, exist_ok=True)
    m = bundle.models[bundle.models["has_pae"]]
    tasks = [(bundle.dir / "pae_files" / r.pae_filename, out / psc.pae_path(r.model_file_id),
              r.protein_end - r.protein_start + 1) for r in m.itertuples()]
    log = work / "pae.tsv"
    done = read_pae_log(log, out)
    # rewritten clean before appending, so the log never carries a torn or duplicate row forward
    part = Path(str(log) + ".part")
    part.write_text(PAE_LOG_HEADER + "".join(f"{k}\t{b}\t{h}\t{sh}\n" for k, (b, h, sh) in done.items()))
    part.rename(log)
    todo = [t for t in tasks if Path(t[1]).stem not in done]
    print(f"pae: {len(tasks)} files, {len(todo)} to convert", flush=True)
    results = list(done.items())
    with mp.Pool(procs) as pool, open(log, "a") as fh:
        for i, r in enumerate(pool.imap_unordered(convert_pae, todo, chunksize=4)):
            fh.write("\t".join(map(str, r)) + "\n")
            results.append(r)
            if i % 2000 == 0:
                fh.flush()
                print(f"  pae {i}/{len(todo)}", flush=True)
    if len(results) != len(tasks):
        raise RuntimeError("PAE conversion incomplete")


# ---------------------------------------------------------------- PDB


def parse_pdb(text: str):
    """CA residues of one model: refuses anything but one chain numbered 1..n with one CA each."""
    pos, res, xyz, bf, chains = [], [], [], [], set()
    for line in text.splitlines():
        if line.startswith("ATOM") and line[12:16] == " CA ":
            if line[16] not in " A":
                continue
            chains.add(line[21])
            pos.append(int(line[22:26]))
            res.append(line[17:20])
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            bf.append(float(line[60:66]))
    if len(chains) != 1:
        raise ValueError(f"expected one chain, found {sorted(chains)}")
    if pos != list(range(1, len(pos) + 1)):
        raise ValueError("residue numbering is not contiguous from 1 with one CA per residue")
    return pos, res, np.array(xyz), np.array(bf)


THREE_TO_ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
                "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
                "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U",
                "PYL": "O"}


def read_model(task):
    src, mid, n_expected = task
    raw = Path(src).read_bytes()
    text = gzip.decompress(raw).decode()
    try:
        pos, res, xyz, bf = parse_pdb(text)
    except ValueError as e:
        raise ValueError(f"{mid}: {e}") from None
    if len(pos) != n_expected:
        raise ValueError(f"{mid}: {len(pos)} residues, guide range has {n_expected}")
    head = text[:20000]
    licence = "CC-BY-4.0" if "CC-BY 4.0" in head else ""
    m = re.search(r"TITLE\s+ALPHAFOLD (\w+) (V[\d.]+) PREDICTION", head)
    engine = f"AlphaFold {m.group(1).title()} {m.group(2).lower()}" if m else ""
    letters = "".join(THREE_TO_ONE.get(r, "X") for r in res)
    return mid, hashlib.sha256(raw).hexdigest(), len(raw), licence, engine, letters, xyz, bf


def stage_models(bundle, out, work, procs):
    (out / psc.STRUCTURE_DIR).mkdir(parents=True, exist_ok=True)
    m = bundle.models
    tasks = [(bundle.dir / "pdb_files" / r.pdb_filename, r.model_file_id, r.protein_end - r.protein_start + 1)
             for r in m.itertuples()]
    # compact columns (model index, one-byte letters, float32 would lose the PDB's 3 decimals so
    # coordinates stay float64): the stage ran into the VM's memory limit holding object arrays
    meta, mids, cols = {}, [], {"model_idx": [], "file_position": [], "residue": [], "ca_x": [], "ca_y": [],
                                "ca_z": [], "plddt": []}
    with mp.Pool(procs) as pool:
        for i, (mid, sha, size, licence, engine, letters, xyz, bf) in enumerate(
                pool.imap(read_model, tasks, chunksize=16)):
            meta[mid] = (sha, size, licence, engine)
            n = len(letters)
            cols["model_idx"].append(np.full(n, len(mids), dtype=np.int32))
            mids.append(mid)
            cols["file_position"].append(np.arange(1, n + 1, dtype=np.int32))
            cols["residue"].append(np.frombuffer(letters.encode(), dtype="S1"))
            for j, c in enumerate(("ca_x", "ca_y", "ca_z")):
                cols[c].append(xyz[:, j])
            cols["plddt"].append(bf)
            if i % 5000 == 0:
                print(f"  models {i}/{len(tasks)}", flush=True)
    idx = np.concatenate(cols.pop("model_idx"))
    arrays = {"model_file_id": pa.DictionaryArray.from_arrays(pa.array(idx), pa.array(mids))}
    arrays["residue"] = pa.array(np.concatenate(cols.pop("residue")).astype(str))
    arrays.update({k: pa.array(np.concatenate(v)) for k, v in cols.items()})
    del cols, idx
    core = pa.table(arrays)
    del arrays
    pq.write_table(core, work / "residue_core.parquet")

    pae = read_pae_log(work / "pae.tsv", out)
    gene_of = bundle.pairs.groupby("accession")["gene_id"].agg(lambda v: sorted(set(v))).to_dict()
    gname = bundle.pairs.drop_duplicates("gene_id").set_index("gene_id")["gene_name"].to_dict()
    rows, mismatched_source = [], 0
    for r in m.itertuples():
        sha, size, licence, engine = meta[r.model_file_id]
        dst = out / psc.structure_path(r.model_file_id)
        if not (dst.exists() and dst.stat().st_size == size):
            shutil.copyfile(bundle.dir / "pdb_files" / r.pdb_filename, str(dst) + ".part")
            Path(str(dst) + ".part").rename(dst)
        source = "AF3" if r.model_file_id.startswith("AF3-") else "AF2DB"
        if (source == "AF2DB") != bundle.af2db[r.uniprot_id]:
            mismatched_source += 1
        if source == "AF2DB":
            if not licence or not engine:
                raise ValueError(f"{r.model_file_id}: AF2DB file without its licence REMARK or TITLE")
            afdb = re.search(r"_(v\d+)$", r.model_file_id).group(1)
            label, version = "AlphaFold Protein Structure Database", f"AFDB {afdb}, {engine}"
        else:
            label, version, licence = psc.AF3_SOURCE_LABEL, psc.AF3_MODEL_VERSION, psc.AF3_LICENCE
        genes = gene_of.get(r.uniprot_id, [])
        has = bool(r.has_pae)
        if has and r.model_file_id not in pae:
            raise ValueError(f"{r.model_file_id}: no intact PAE object and log row; rerun the pae stage")
        p = dict(zip(("pae_bytes", "pae_sha256", "source_pae_sha256"), pae[r.model_file_id])) if has else None
        rows.append({
            "model_file_id": r.model_file_id, "accession": r.uniprot_id,
            "uniprot_isoform": bundle.isoform[r.uniprot_id], "gene_id": ",".join(genes),
            "gene_name": ",".join(gname[g] for g in genes), "source_gene_symbol": bundle.source_symbols[r.uniprot_id],
            "source": source, "source_label": label, "model_version": version, "licence": licence,
            "window_index": r.window_index, "n_windows": r.n_windows,
            "protein_start": r.protein_start, "protein_end": r.protein_end,
            "n_residues": r.protein_end - r.protein_start + 1,
            "owned_start": r.owned_start, "owned_end": r.owned_end, "has_pae": has,
            "structure_path": psc.structure_path(r.model_file_id), "structure_sha256": sha,
            "structure_bytes": size,
            "pae_path": psc.pae_path(r.model_file_id) if has else "",
            "pae_sha256": p["pae_sha256"] if has else "", "pae_bytes": p["pae_bytes"] if has else "",
            "source_pdb_sha256": sha, "source_pae_sha256": p["source_pae_sha256"] if has else "",
        })
    if mismatched_source:
        raise ValueError(f"{mismatched_source} model files whose prefix disagrees with the sequence guide's AF2DB flag")
    df = pd.DataFrame(rows)[[c for c, _ in psc.MODEL_COLUMNS]]
    text = df.assign(has_pae=tf(df["has_pae"])).to_csv(sep="\t", header=False, index=False, lineterminator="\n")
    check_model_rows(text)
    part = Path(str(out / psc.MODELS_FILE) + ".part")
    # mtime 0 so a rebuild from the same inputs is byte-identical
    with open(part, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as fh:
        fh.write((psc.header_line(psc.MODEL_COLUMNS) + text).encode())
    part.rename(out / psc.MODELS_FILE)
    write_json(work / "models.json", {
        "models": len(df), "proteins": int(df["accession"].nunique()),
        "proteins_with_pae": int(df[df["has_pae"]]["accession"].nunique()),
        "proteins_without_pae": int(df[~df["has_pae"]]["accession"].nunique()),
        "models_by_source": df["source"].value_counts().to_dict(),
        "models_gene_not_in_gencode": int((df["gene_id"] == "").sum()),
        "models_under_several_genes": int(df["gene_id"].str.contains(",").sum()),
        "structure_bytes_total": int(df["structure_bytes"].sum()),
        "pae_bytes_total": int(pd.to_numeric(df["pae_bytes"], errors="coerce").fillna(0).sum()),
        "models_residue_differs_from_sequence": sequence_mismatches(bundle, core),
    })


def sequence_mismatches(bundle, core):
    """{model_file_id: ["<file position>:<model letter>><sequence letter>", ...]} where the model and the sequence guide differ."""
    m = bundle.models.set_index("model_file_id")
    t = core.select(["model_file_id", "residue"]).to_pandas()
    out = {}
    for mid, letters in t.groupby("model_file_id", sort=False)["residue"]:
        r = m.loc[mid]
        seq = bundle.sequence[r.uniprot_id][r.protein_start - 1:r.protein_end]
        diff = [f"{i}:{a}>{b}" for i, (a, b) in enumerate(zip("".join(letters), seq), 1) if a != b]
        if diff:
            out[mid] = diff
    return out


_INT_COLUMNS = {c for c, t in psc.MODEL_COLUMNS if t is int}
_PAE_COLUMNS = ("pae_path", "pae_sha256", "pae_bytes", "source_pae_sha256")


def check_model_rows(text):
    """refuses a models.tsv body that breaks the contract's column rules, naming the first bad row."""
    names = [c for c, _ in psc.MODEL_COLUMNS]
    windows = {}
    for i, line in enumerate(text.rstrip("\n").split("\n"), 2):
        f = line.split("\t")

        def bad(why, f=f, i=i):
            raise ValueError(f"models.tsv line {i} ({f[0]}): {why}")
        if len(f) != len(names):
            bad(f"{len(f)} fields, the contract has {len(names)}")
        r = dict(zip(names, f))
        if not psc.MODEL_FILE_ID_RE.match(r["model_file_id"]):
            bad("model_file_id outside the allowed charset")
        if r["has_pae"] not in ("true", "false"):
            bad(f"has_pae {r['has_pae']!r}")
        has = r["has_pae"] == "true"
        for c in _INT_COLUMNS:
            if not (r[c].isdigit() or (c == "pae_bytes" and not has and r[c] == "")):
                bad(f"{c} {r[c]!r} is not an integer")
        for c in ("structure_sha256", "source_pdb_sha256") + (("pae_sha256", "source_pae_sha256") if has else ()):
            if not SHA256_RE.match(r[c]):
                bad(f"{c} {r[c]!r} is not a sha256")
        if not has and any(r[c] for c in _PAE_COLUMNS):
            bad("PAE fields set without has_pae")
        if r["source"] not in psc.SOURCES:
            bad(f"source {r['source']!r}")
        if r["structure_path"] != psc.structure_path(r["model_file_id"]):
            bad("structure_path is not the derived path")
        if has and r["pae_path"] != psc.pae_path(r["model_file_id"]):
            bad("pae_path is not the derived path")
        if has and int(r["pae_bytes"]) != int(r["n_residues"]) ** 2:
            bad("pae_bytes is not n_residues squared")
        if int(r["n_residues"]) != int(r["protein_end"]) - int(r["protein_start"]) + 1:
            bad("n_residues does not match the protein range")
        psc.split_model_genes(r["gene_id"], r["gene_name"])
        windows.setdefault(r["accession"], []).append((int(r["window_index"]), int(r["n_windows"])))
    for acc, w in windows.items():
        if sorted(k for k, _ in w) != list(range(1, len(w) + 1)) or {n for _, n in w} != {len(w)}:
            raise ValueError(f"models.tsv: {acc}'s window_index is not 1..n_windows")


# ---------------------------------------------------------------- tabix writing


class TabixWriter:
    """header + rows piped through bgzip, then indexed with the contract's tabix arguments."""

    def __init__(self, path, columns, tabix_args):
        self.path, self.columns, self.tabix_args = Path(path), columns, tabix_args
        self.part = self.path.with_name(self.path.name + ".part")
        self.fh = open(self.part, "wb")
        self.proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=self.fh)
        self.proc.stdin.write(psc.header_line(columns).encode())
        self.rows = 0

    def write(self, df, float_format=None):
        self.proc.stdin.write(df[[c for c, _ in self.columns]].to_csv(
            sep="\t", header=False, index=False, na_rep="", lineterminator="\n",
            float_format=float_format).encode())
        self.rows += len(df)

    def close(self):
        """finish and index under the temporary name; nothing looks complete until commit()."""
        self.proc.stdin.close()
        if self.proc.wait() != 0:
            raise RuntimeError("bgzip failed")
        self.fh.close()
        subprocess.run(["tabix", "-f", *self.tabix_args, str(self.part)], check=True)
        return self.rows

    def commit(self):
        Path(str(self.part) + ".tbi").rename(str(self.path) + ".tbi")
        self.part.rename(self.path)


def tf(s):
    return np.where(s, "true", "false")


# ---------------------------------------------------------------- residues


def stage_residues(bundle, out, work, lab):
    """residue rows, one accession at a time from numpy arrays: a whole-table frame did not fit on the VM."""
    # columns read one at a time and mu kept as slices of flat arrays: this stage runs on a VM
    # shared with other jobs, and a frame or a dict of copies did not fit in the memory it had
    core_path = work / "residue_core.parquet"
    mids = pq.read_table(core_path, columns=["model_file_id"]).column(0)
    # the core is in guide order, one contiguous block per model file
    m = bundle.models.assign(n=bundle.models["protein_end"] - bundle.models["protein_start"] + 1)
    m["offset"] = np.concatenate([[0], np.cumsum(m["n"].to_numpy())[:-1]])
    if int(m["n"].sum()) != len(mids) or mids[int(m["offset"].iat[-1])].as_py() != m["model_file_id"].iat[-1]:
        raise ValueError("residue core does not follow the guide's model order")
    del mids
    letters = pq.read_table(core_path, columns=["residue"]).column(0).combine_chunks()
    if letters.buffers()[2].size < len(letters) or pa.compute.max(pa.compute.utf8_length(letters)).as_py() != 1:
        raise ValueError("residue core letters are not one character each")
    residue = np.frombuffer(letters.buffers()[2], dtype="S1", count=len(letters)).astype("U1")
    del letters
    xyz = np.column_stack([pq.read_table(core_path, columns=[c]).column(0).to_numpy() for c in ("ca_x", "ca_y", "ca_z")])
    plddt = pq.read_table(core_path, columns=["plddt"]).column(0).to_numpy()
    models_of = {a: x.sort_values("model_file_id") for a, x in m.groupby("uniprot_id", sort=False)}

    mt = pq.read_table(lab / "mu_per_residue_gnomad.parquet")
    mt = mt.filter(pa.compute.is_in(mt.column("uniprot_id"), value_set=pa.array(list(models_of))))
    mt = mt.sort_by([("uniprot_id", "ascending"), ("aa_pos", "ascending")])
    acc = mt.column("uniprot_id").combine_chunks().dictionary_encode()
    codes = acc.indices.to_numpy()
    mu_pos = mt.column("aa_pos").to_numpy()
    mu_val, mu_nw, mu_np = (mt.column(c).to_numpy() for c in ("mu", "n_with_mu", "n_possible"))
    del mt
    starts = np.flatnonzero(np.r_[True, codes[1:] != codes[:-1]]) if len(codes) else np.array([], int)
    ends = np.r_[starts[1:], len(codes)]
    names = acc.dictionary.to_pylist()
    mu_of = {names[codes[a]]: (a, b) for a, b in zip(starts, ends)}
    del codes, acc
    dup = np.flatnonzero(np.diff(mu_pos) == 0)
    if any(np.searchsorted(ends, i, side="right") == np.searchsorted(ends, i + 1, side="right") for i in dup):
        raise ValueError("mu_per_residue has duplicate (uniprot_id, aa_pos)")
    check_mu_format(mu_val, "mu")

    def body(acc):
        """the residue lines of one accession after the gene anchor, and its owner count per position."""
        lines, owners = [], {}
        a, b = mu_of.get(acc, (0, 0))
        pos_mu, val_mu, nw, npo = mu_pos[a:b], mu_val[a:b], mu_nw[a:b], mu_np[a:b]
        for r in models_of[acc].itertuples():
            fp = np.arange(1, r.n + 1)
            pp = r.protein_start + fp - 1
            own = (pp >= r.owned_start) & (pp <= r.owned_end)
            for p in pp[own]:
                owners[p] = owners.get(p, 0) + 1
            k = np.searchsorted(pos_mu, pp)
            hit = (k < len(pos_mu)) & (pos_mu[np.minimum(k, len(pos_mu) - 1)] == pp) if len(pos_mu) else np.zeros(r.n, bool)
            sl = slice(r.offset, r.offset + r.n)
            for i, (f, q, o, aa, (x, y, z), b) in enumerate(zip(fp, pp, own, residue[sl], xyz[sl], plddt[sl])):
                mus = (MU_FORMAT.format(val_mu[k[i]]), str(nw[k[i]]), str(npo[k[i]])) if hit[i] else ("", "", "")
                lines.append(f"{acc}\t{r.model_file_id}\t{f}\t{q}\t{aa}\t{'true' if o else 'false'}\t"
                             f"{x:.3f}\t{y:.3f}\t{z:.3f}\t{b:.2f}\t{mus[0]}\t{mus[1]}\t{mus[2]}\n")
        return lines, owners

    w = TabixWriter(out / psc.RESIDUES_FILE, psc.RESIDUE_COLUMNS, psc.GENE_TABIX_ARGS)
    gene_bytes, owner_rows, bad_owner = {}, 0, 0
    for pair in bundle.pairs.itertuples():
        lines, owners = body(pair.accession)
        bad_owner += sum(1 for v in owners.values() if v != 1)
        bad_owner += max(models_of[pair.accession]["protein_end"]) - len(owners)
        owner_rows += len(owners)
        prefix = f"{pair.gene_id}\t{pair.chrom}\t{pair.gene_start}\t{pair.gene_end}\t{pair.gene_name}\t"
        data = "".join(prefix + ln for ln in lines).encode()
        gene_bytes[pair.gene_id] = gene_bytes.get(pair.gene_id, 0) + len(data)
        w.proc.stdin.write(data)
        w.rows += len(lines)
    rows = w.close()
    if bad_owner:
        raise AssertionError(f"is_owner is not exactly one row for {bad_owner} (gene, accession, position)s")
    largest = max(gene_bytes, key=gene_bytes.get)
    if gene_bytes[largest] > MAX_GENE_RESIDUE_BYTES:
        raise AssertionError(f"gene {largest}'s residue rows are {gene_bytes[largest]} bytes, over the 8 MiB gate")
    w.commit()
    write_json(work / "residues.json", {
        "residue_rows": rows, "residue_owner_rows": owner_rows,
        "residue_gene_largest": {"gene_id": largest, "bytes": gene_bytes[largest], "limit": MAX_GENE_RESIDUE_BYTES}})


# ---------------------------------------------------------------- features


def eco_codes(evidence):
    return ",".join(dict.fromkeys(ECO_RE.findall(evidence or "")))


def uniprot_features(bundle, path, drops):
    f = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    n0 = len(f)
    f = f[f["uniprot_id"].isin(bundle.pairs["accession"])]
    drops["uniprot_rows_accession_not_in_release_or_gencode"] = n0 - len(f)

    def same_isoform(acc, iso):
        g = bundle.isoform[acc]
        return bundle.canonical[acc] and (g if "-" in g else g + "-1") == iso

    keep = np.array([same_isoform(a, i) for a, i in zip(f["uniprot_id"], f["canonical_isoform"])], dtype=bool)
    drops["uniprot_rows_bundle_sequence_not_canonical_isoform"] = int((~keep).sum())
    drops["uniprot_accessions_bundle_sequence_not_canonical_isoform"] = int(f.loc[~keep, "uniprot_id"].nunique())
    f = f[keep]
    known = (f["position_start"] != "") & (f["position_end"] != "")
    drops["uniprot_rows_position_unknown"] = int((~known).sum())
    f = f[known].copy()
    for c in ("position_start", "position_end"):
        v = f[c].astype(float)
        if (v != v.round()).any():
            raise ValueError(f"non-integer {c} in the UniProt feature file")
        f[c] = v.astype(int)
    length = f["uniprot_id"].map(lambda a: len(bundle.sequence[a]))
    ok = (f["position_start"] >= 1) & (f["position_start"] <= f["position_end"]) & (f["position_end"] <= length)
    drops["uniprot_rows_position_outside_sequence"] = int((~ok).sum())
    f = f[ok]
    codes = set(ECO_RE.findall(" ".join(f["evidence"])))
    out = pd.DataFrame({
        "accession": f["uniprot_id"], "source": "uniprot", "feature_type": f["feature"],
        "start": f["position_start"], "end": f["position_end"], "description": f["description"],
        "evidence": f["evidence"].map(eco_codes)})
    out["evidence_tier"] = out["evidence"].map(lambda e: psc.evidence_tier("uniprot", e))
    unknown = sorted(c for c in codes if psc.evidence_tier("uniprot", c) == "unknown")
    return out, {"eco_codes_present": sorted(codes), "eco_codes_unknown_tier": unknown,
                 "rows_without_evidence": int((out["evidence"] == "").sum())}


def ted_features(bundle, path, drops):
    t = pd.read_csv(path, sep="\t", header=None, names=TED_COLUMNS, dtype=str, keep_default_na=False)
    ids = t["ted_id"].str.extract(TED_ID_RE)
    if ids.isna().any().any():
        raise ValueError("unrecognised TED domain id")
    t["accession"], t["fragment"], t["domain"] = ids[0], ids[1].astype(int), ids[2]
    stats = {"ted_rows_human": len(t), "ted_fragments_seen": sorted(int(v) for v in t["fragment"].unique())}
    t = t[t["accession"].isin(bundle.pairs["accession"])]
    drops["ted_domains_accession_not_in_release_or_gencode"] = stats["ted_rows_human"] - len(t)
    rows, bad_md5 = [], 0
    for r in t.itertuples():
        seq = bundle.sequence[r.accession]
        off = (r.fragment - 1) * AFDB_FRAGMENT_STEP
        segs = [tuple(int(v) + off for v in s.split("-")) for s in r.chopping.split("_")]
        if any(a < 1 or b > len(seq) or a > b for a, b in segs) or \
                hashlib.md5("".join(seq[a - 1:b] for a, b in segs).encode()).hexdigest() != r.md5_domain:
            bad_md5 += 1
            continue
        cath = "" if r.cath_label == "-" else f" CATH {r.cath_label} ({r.cath_assignment_level}, {r.cath_assignment_method})"
        for i, (a, b) in enumerate(segs, 1):
            seg = f" segment {i}/{len(segs)}" if len(segs) > 1 else ""
            rows.append((r.accession, "ted", "TED domain", a, b,
                         f"TED{r.domain} {r.consensus_level} consensus{cath}{seg}", ""))
    drops["ted_domains_sequence_md5_differs"] = bad_md5
    stats["ted_domains_kept"] = int(t.shape[0] - bad_md5)
    stats["ted_accessions_kept"] = len({r[0] for r in rows})
    out = pd.DataFrame(rows, columns=["accession", "source", "feature_type", "start", "end", "description", "evidence"])
    out["evidence_tier"] = psc.evidence_tier("ted", "")
    return out, stats


def stage_features(bundle, out, work, lab, ted):
    drops = {}
    u, ustats = uniprot_features(bundle, lab / "uniprot_all_feats_evidence.tsv", drops)
    parts, stats = [u], {"uniprot": ustats}
    if ted:
        tdf, stats["ted"] = ted_features(bundle, ted, drops)
        parts.append(tdf)
    f = pd.concat(parts, ignore_index=True)
    n = len(f)
    f = f.drop_duplicates()
    drops["feature_rows_exact_duplicate"] = n - len(f)
    f = bundle.pairs.merge(f, on="accession", how="inner")
    f = f.sort_values(list(psc.FEATURE_SORT), kind="stable")
    w = TabixWriter(out / psc.FEATURES_FILE, psc.FEATURE_COLUMNS, psc.GENE_TABIX_ARGS)
    w.write(f)
    w.close()
    w.commit()
    stats["feature_rows_by_source"] = f["source"].value_counts().to_dict()
    stats["feature_rows_by_tier"] = {f"{a}/{b}": int(n) for (a, b), n in
                                     f.groupby(["source", "evidence_tier"]).size().items()}
    stats["feature_rows"] = len(f)
    stats["drops"] = drops
    write_json(work / "features.json", stats)


# ---------------------------------------------------------------- variants


def load_mu_snp(path, chrom):
    import pyarrow.dataset as ds
    t = ds.dataset(path).to_table(filter=ds.field("chrom") == f"chr{chrom}", columns=["pos", "ref", "alt", "mu_snp"])
    # kept as Arrow and sliced per block: a whole chromosome as a frame is several hundred MB
    return t.sort_by("pos")


def stage_variants(bundle, out, work, lab):
    import h5py
    import hdf5plugin  # noqa: F401  registers the HDF5 filters the map is compressed with

    m = bundle.models
    ranges = m[["uniprot_id", "model_file_id", "protein_start", "protein_end", "owned_start", "owned_end"]]
    proteins = set(m["uniprot_id"])
    plen = {a: len(s) for a, s in bundle.sequence.items()}
    pairs = bundle.pairs
    stats = {"rows_in_map": 0, "rows_address_differs_from_map": 0, "rows_map_offset_inconsistent": 0,
             "rows_aa_ref_differs_from_sequence": 0, "rows_no_mu_snp": 0, "rows_gene_ambiguous": 0,
             "rows_no_gene": 0, "variants_with_several_accessions": 0}
    drops = {"rows_accession_without_model": 0, "rows_aa_pos_outside_models": 0,
             "rows_duplicate_variant_accession": 0}
    dropped_acc = set()
    w = TabixWriter(out / psc.VARIANTS_FILE, psc.VARIANT_COLUMNS, psc.VARIANT_TABIX_ARGS)
    with h5py.File(bundle.dir / "all_missense_variants_gr38.h5", "r") as h:
        for c in CHROMS:
            t0 = time.time()
            # kept as fixed-width numpy bytes and cut into position-ordered blocks: a whole
            # chromosome as a frame, merged and sorted, does not fit beside everything else on the VM
            ra_all = h[f"chr{c}_ref_alt"][:]
            p_all = h[f"chr{c}_pos"][:]
            acc_all = h[f"chr{c}_uniprot_id"][:, 0]
            file_all = h[f"chr{c}_filename"][:, 0]
            order = np.argsort(p_all[:, 0], kind="stable")
            sorted_pos = p_all[order, 0]
            mu_all = load_mu_snp(lab / "gnomad_mu_snp_by_variant.parquet", c)
            mu_pos = mu_all.column("pos").to_numpy()
            bounds = [0]
            while bounds[-1] < len(order):
                nxt = min(bounds[-1] + VARIANT_BLOCK_ROWS, len(order))
                if nxt < len(order):
                    nxt = int(np.searchsorted(sorted_pos, sorted_pos[nxt], side="left")) or nxt
                    if nxt <= bounds[-1]:
                        nxt = int(np.searchsorted(sorted_pos, sorted_pos[bounds[-1]], side="right"))
                bounds.append(nxt)
            for lo, hi in zip(bounds[:-1], bounds[1:]):
                idx = order[lo:hi]
                ra = ra_all[idx].astype(str)
                p = p_all[idx]
                v = pd.DataFrame({"pos": p[:, 0], "aa_pos": p[:, 1], "aa_pos_file": p[:, 2],
                                  "ref": ra[:, 0], "alt": ra[:, 1], "aa_ref": ra[:, 2], "aa_alt": ra[:, 3],
                                  "accession": acc_all[idx].astype(str), "map_file": file_all[idx].astype(str)})
                del ra, p
                stats["rows_in_map"] += len(v)
                inrel = v["accession"].isin(proteins)
                drops["rows_accession_without_model"] += int((~inrel).sum())
                dropped_acc |= set(v.loc[~inrel, "accession"])
                v = v[inrel]
                if v.empty:
                    continue
                n = len(v)
                v = v.drop_duplicates(["pos", "ref", "alt", "accession"])
                drops["rows_duplicate_variant_accession"] += n - len(v)
                # the map's own window choice, and its offset, which must be protein_start - 1
                mf = v["map_file"].str.replace(r"\.pdb\.gz$", "", regex=True)
                mstart = mf.map(m.set_index("model_file_id")["protein_start"])
                stats["rows_map_offset_inconsistent"] += int(((v["aa_pos"] - v["aa_pos_file"]) != (mstart - 1)).sum())
                n = len(v)
                v = v.assign(map_mid=mf).merge(ranges, left_on="accession", right_on="uniprot_id")
                v = v[(v["aa_pos"] >= v["owned_start"]) & (v["aa_pos"] <= v["owned_end"])].reset_index(drop=True)
                drops["rows_aa_pos_outside_models"] += n - len(v)
                stats["rows_address_differs_from_map"] += int((v["model_file_id"] != v["map_mid"]).sum())
                v["file_position"] = v["aa_pos"] - v["protein_start"] + 1
                seq_ref = [bundle.sequence[a][i - 1] if i <= plen[a] else "" for a, i in zip(v["accession"], v["aa_pos"])]
                stats["rows_aa_ref_differs_from_sequence"] += int((v["aa_ref"].to_numpy() != np.array(seq_ref)).sum())
                chrom = 23 if c == "X" else int(c)
                v["chrom"] = chrom
                v["variant"] = f"{chrom}:" + v["pos"].astype(str) + ":" + v["ref"] + ":" + v["alt"]
                # the accession's gene whose span contains pos; empty when none or several do
                g = v[["pos", "accession"]].reset_index().merge(pairs[pairs["chrom"] == chrom], on="accession")
                g = g[(g["pos"] >= g["gene_start"]) & (g["pos"] <= g["gene_end"])]
                hits = g.groupby("index")["gene_id"].agg(["count", "first"])
                name = pairs.drop_duplicates("gene_id").set_index("gene_id")["gene_name"]
                one = hits[hits["count"] == 1]["first"]
                v["gene_id"] = one.reindex(v.index).fillna("")
                v["gene_name"] = v["gene_id"].map(name).fillna("")
                stats["rows_gene_ambiguous"] += int((hits["count"] > 1).sum())
                stats["rows_no_gene"] += int(len(v) - len(hits))
                a, b = np.searchsorted(mu_pos, [v["pos"].min(), v["pos"].max() + 1])
                mu = mu_all.slice(a, b - a).to_pandas()
                if mu.duplicated(["pos", "ref", "alt"]).any():
                    raise ValueError(f"chr{c}: duplicate variants in the mu_snp file")
                check_mu_format(mu["mu_snp"].to_numpy(), f"chr{c} mu_snp")
                v = v.merge(mu, on=["pos", "ref", "alt"], how="left")
                stats["rows_no_mu_snp"] += int(v["mu_snp"].isna().sum())
                del mu
                stats["variants_with_several_accessions"] += int(
                    (v.groupby(["pos", "ref", "alt"]).size() > 1).sum())
                v = v.sort_values(list(psc.VARIANT_SORT))
                w.write(v, float_format="%.10g")
            print(f"  variants chr{c}: {hi} map rows, {time.time() - t0:.0f}s", flush=True)
            del ra_all, p_all, acc_all, file_all, order, sorted_pos, mu_all, mu_pos
    stats["variant_rows"] = w.close()
    w.commit()
    stats["accessions_dropped_without_model"] = sorted(dropped_acc)
    stats["drops"] = drops
    write_json(work / "variants.json", stats)


# ---------------------------------------------------------------- manifest and upload


def stage_manifest(bundle, out, work, args, release):
    cached = work / "inputs_sha256.json"
    sums = json.loads(cached.read_text()) if cached.exists() else {}

    def h(_name, path):
        st = Path(path).stat()
        key = f"{Path(path).resolve()}|{st.st_size}|{st.st_mtime_ns}"
        if key not in sums:
            sums[key] = sha256_file(path)
            write_json(cached, sums)
        return sums[key]

    lab = Path(args.lab_files)
    guide = bundle.dir / "pdb_pae_file_pos_guide.tsv"
    manifest = {
        "contract_version": psc.CONTRACT_VERSION, "release": release,
        "built_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "bundle": {
            "tar": Path(args.bundle_tar).name if args.bundle_tar else None,
            "tar_sha256": h("bundle_tar", args.bundle_tar) if args.bundle_tar else None,
            # the tar's mtime is its download time; the guide's is when the lab wrote the bundle
            "date": dt.datetime.fromtimestamp(guide.stat().st_mtime, dt.timezone.utc).date().isoformat(),
        },
        "gencode": {"version": psc.GENCODE_VERSION, "file": Path(args.gencode).name,
                    "sha256": h("gencode", args.gencode)},
        "lab_inputs": {n: {"sha256": h(n, lab / n)} for n in psc.LAB_INPUTS},
        "feature_sources": {
            "uniprot": {"file": "uniprot_all_feats_evidence.tsv", "version": "not stated by source",
                        "licence": "CC-BY-4.0"},
            "ted": json.loads(Path(args.ted_meta).read_text()) if args.ted_meta else None,
        },
        "no_model_accessions": bundle.no_model_accessions,
        "accessions_in_several_genes": bundle.accessions_in_several_genes,
        "genes_absent_from_gencode": bundle.genes_absent_from_gencode,
        "drops": dict(bundle.drops),
        "counts": {},
        "objects": {},
        "directories": {},
    }
    for stage in STAGE_JSON:
        if not (work / f"{stage}.json").exists():
            raise RuntimeError(f"stage {stage} has not finished: no _work/{stage}.json")
        s = json.loads((work / f"{stage}.json").read_text())
        manifest["drops"].update({f"{stage}.{k}": v for k, v in s.pop("drops", {}).items()})
        manifest["counts"][stage] = s
    for f in psc.manifest_object_names():
        p = out / f
        manifest[psc.MANIFEST_OBJECTS][f] = {psc.OBJECT_BYTES: p.stat().st_size, psc.OBJECT_SHA256: sha256_file(p)}
    for d in (psc.STRUCTURE_DIR, psc.PAE_DIR):
        files = list((out / d).glob("*"))
        manifest[psc.MANIFEST_DIRECTORIES][d] = {psc.DIRECTORY_OBJECTS: len(files),
                                                 psc.OBJECT_BYTES: sum(f.stat().st_size for f in files)}
    if set(manifest) - {"no_model_accessions", "accessions_in_several_genes", "genes_absent_from_gencode"} \
            != set(psc.MANIFEST_KEYS):
        raise AssertionError("manifest keys differ from the contract's MANIFEST_KEYS")
    write_json(out / psc.MANIFEST_FILE, manifest)
    if (out / (psc.MANIFEST_FILE + ".part")).exists():
        raise RuntimeError("manifest not written")


def check_root(root):
    root = root.rstrip("/")
    if root not in (REHEARSAL_ROOT, PRODUCTION_ROOT):
        raise ValueError(f"{root} is neither the rehearsal nor the production root")
    return root


# the bucket side of the upload, kept to these four calls so tests can stand in for the bucket
def gs_list(prefix):
    """{path relative to prefix: (bytes, crc32c base64)} for every object under the prefix."""
    r = subprocess.run(["gcloud", "storage", "objects", "list", f"{prefix}/**",
                        "--format=json(name,size,crc32c_hash)"], capture_output=True, text=True)
    if r.returncode:
        if "matched no objects" in r.stderr:
            return {}
        raise RuntimeError(f"listing {prefix} failed: {r.stderr.strip()}")
    base = prefix.split("/", 3)[3] + "/"
    return {o["name"][len(base):]: (int(o["size"]), o["crc32c_hash"]) for o in json.loads(r.stdout or "[]")}


def gs_cp(args, dest):
    subprocess.run(["gcloud", "storage", "cp", "--no-clobber", *args, dest], check=True)


def gs_cat(uri):
    return subprocess.run(["gcloud", "storage", "cat", uri], check=True, capture_output=True).stdout


def gs_replace(local, uri):
    subprocess.run(["gcloud", "storage", "cp", str(local), uri], check=True)


def expected_objects(out):
    """{release-relative path: (bytes, crc32c)} for everything the manifest names, verified locally.

    Top-level files against the manifest's bytes and sha256, structures and PAE objects against
    their models.tsv.gz rows, and the directory totals against those rows.
    """
    man_bytes = (out / psc.MANIFEST_FILE).read_bytes()
    man = json.loads(man_bytes)
    exp = {psc.MANIFEST_FILE: (len(man_bytes), bytes_crc32c(man_bytes))}
    for name, o in man[psc.MANIFEST_OBJECTS].items():
        if not (out / name).exists():
            raise RefuseUpload(f"{name}: named by the manifest but missing locally")
        n, sha, crc = file_digests(out / name)
        if (n, sha) != (o[psc.OBJECT_BYTES], o[psc.OBJECT_SHA256]):
            raise RefuseUpload(f"{name}: local file does not match the manifest")
        exp[name] = (n, crc)
    models = pd.read_csv(out / psc.MODELS_FILE, sep="\t", dtype=str, keep_default_na=False)
    models.columns = [c for c, _ in psc.MODEL_COLUMNS]
    totals = {psc.STRUCTURE_DIR: [0, 0], psc.PAE_DIR: [0, 0]}
    for r in models.itertuples():
        objs = [(r.structure_path, r.structure_bytes, r.structure_sha256, psc.STRUCTURE_DIR)]
        if r.has_pae == "true":
            objs.append((r.pae_path, r.pae_bytes, r.pae_sha256, psc.PAE_DIR))
        for path, size, sha, d in objs:
            if not (out / path).exists():
                raise RefuseUpload(f"{path}: named by models.tsv.gz but missing locally")
            n, s256, crc = file_digests(out / path)
            if (n, s256) != (int(size), sha):
                raise RefuseUpload(f"{path}: local file does not match its models.tsv.gz row")
            exp[path] = (n, crc)
            totals[d][0] += 1
            totals[d][1] += n
    for d, (count, size) in totals.items():
        o = man[psc.MANIFEST_DIRECTORIES][d]
        if (o[psc.DIRECTORY_OBJECTS], o[psc.OBJECT_BYTES]) != (count, size):
            raise RefuseUpload(f"{d}/: the manifest's totals do not match models.tsv.gz")
    return exp, man_bytes


def check_bucket(dest, exp, complete):
    remote = gs_list(dest)
    for rel, got in sorted(remote.items()):
        if rel not in exp:
            raise RefuseUpload(f"{dest}/{rel} is in the bucket but not named by this manifest")
        if got != exp[rel]:
            if rel == psc.MANIFEST_FILE:
                raise RefuseUpload(f"{dest}/{rel} differs from the local manifest: it was rebuilt after an "
                                   "upload started; build under a new release id")
            raise RefuseUpload(f"{dest}/{rel} differs from the local object (size or crc32c)")
    if complete:
        missing = sorted(set(exp) - set(remote))
        if missing:
            raise RefuseUpload(f"{len(missing)} objects missing from {dest}, e.g. {missing[0]}")


def upload(out, release, root, copy_from=None):
    root = check_root(root)
    dest = f"{root}/releases/{release}"
    if not (out / psc.MANIFEST_FILE).exists():
        raise RuntimeError("no manifest: run the manifest stage first")
    exp, man_bytes = expected_objects(out)
    check_bucket(dest, exp, complete=False)
    dirs = [psc.STRUCTURE_DIR, psc.PAE_DIR]
    if copy_from:
        gs_cp(["-r", *(f"{copy_from.rstrip('/')}/{d}" for d in dirs)], f"{dest}/")
    else:
        gs_cp(["-r", *(str(out / d) for d in dirs)], f"{dest}/")
    top = [n for n in exp if "/" not in n and n != psc.MANIFEST_FILE]
    gs_cp([str(out / n) for n in top], f"{dest}/")
    gs_cp([str(out / psc.MANIFEST_FILE)], f"{dest}/{psc.MANIFEST_FILE}")
    check_bucket(dest, exp, complete=True)
    manifest_sha = hashlib.sha256(man_bytes).hexdigest()
    if hashlib.sha256(gs_cat(f"{dest}/{psc.MANIFEST_FILE}")).hexdigest() != manifest_sha:
        raise RefuseUpload(f"{dest}/{psc.MANIFEST_FILE} read back differs from the local manifest")
    ptr = out.parent / f"{release}.CURRENT"
    ptr.write_text(psc.pointer(release, dest, dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                               manifest_sha))
    # the pointer is the one object that is replaced, and it goes last
    gs_replace(ptr, f"{root}/{psc.POINTER}")
    print(ptr.read_text())


# ---------------------------------------------------------------- main


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--lab-files", required=True)
    ap.add_argument("--gencode", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--bundle-tar")
    ap.add_argument("--ted")
    ap.add_argument("--ted-meta")
    ap.add_argument("--release")
    ap.add_argument("--stages", default=",".join(STAGES))
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--upload", action="store_true")
    ap.add_argument("--root", default=REHEARSAL_ROOT)
    ap.add_argument("--copy-from")
    args = ap.parse_args(argv)
    root = check_root(args.root)

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    rel_file = workdir / "RELEASE"
    if rel_file.exists():
        release = rel_file.read_text().strip()
        if args.release and args.release != release:
            raise ValueError(f"workdir holds release {release}")
    else:
        release = args.release or dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d-%H%M%S")
        rel_file.write_text(release + "\n")
    out = workdir / release
    work = workdir / "_work"
    out.mkdir(exist_ok=True)
    work.mkdir(exist_ok=True)

    stages = [s for s in args.stages.split(",") if s]
    unknown = set(stages) - set(STAGES)
    if unknown:
        raise ValueError(f"unknown stages {sorted(unknown)}")
    gencode = load_gencode(args.gencode)
    bundle = Bundle(args.bundle, gencode)
    lab = Path(args.lab_files)
    timings = json.loads((work / "timings.json").read_text()) if (work / "timings.json").exists() else {}
    for s in STAGES:
        if s not in stages:
            continue
        t0 = time.time()
        print(f"stage {s}", flush=True)
        if s in STAGE_JSON:
            (work / f"{s}.json").unlink(missing_ok=True)
        if s == "pae":
            stage_pae(bundle, out, work, args.procs)
        elif s == "models":
            stage_models(bundle, out, work, args.procs)
        elif s == "residues":
            stage_residues(bundle, out, work, lab)
        elif s == "features":
            stage_features(bundle, out, work, lab, args.ted)
        elif s == "variants":
            stage_variants(bundle, out, work, lab)
        elif s == "manifest":
            stage_manifest(bundle, out, work, args, release)
        timings[s] = round(time.time() - t0, 1)
        write_json(work / "timings.json", timings)
        print(f"stage {s} done in {timings[s]}s", flush=True)
    if args.upload:
        upload(out, release, root, args.copy_from)
    return out


if __name__ == "__main__":
    main()
