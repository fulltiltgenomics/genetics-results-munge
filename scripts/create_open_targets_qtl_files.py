"""
Convert the QTL credible sets of an Open Targets release to the credible set TSV format used by
the rest of this repository, laid out like the eQTL Catalogue munge.

Input
  <data_dir>/credible_set/*.parquet    credible set parquet files of an Open Targets release
  <data_dir>/study_metadata/*.parquet  its study table
  <data_dir>/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz
                                       eQTL Catalogue gene metadata, for gene_id -> gene_name
  metadata/eqtl_catalogue_studies.tsv  eQTL Catalogue study table, for the GTEx tissue labels

Output, under <data_dir>/opentargets_qtl_per_study/
  <sub-study>.SUSIE.munged.tsv               one file per project x tissue/cell type x
                                             quantification method, sorted by position

No credible set stats are written, unlike the eQTL Catalogue munge: per gene they come to more
than a million stats.json files and hours of compute here and twice more in annotate_resource.sh,
and nothing reads them for this dataset (results-api's optional stats_file is not set for it).

Only the projects in PROJECTS are kept: the other QTL projects of the release are eQTL
Catalogue studies, which are munged from eQTL Catalogue itself.

A study of the release is one molecular trait in one tissue, so the release has millions of them;
the files are split per sub-study the way eQTL Catalogue's QTD ids split them instead. The
columns follow create_eqtl_catalogue_files.py: `trait` is the gene name, `trait_original` is
`<molecular_trait_id>|<quant_method>`, `cell_type` is `<tissue>|<condition>` and `data_type` is
eQTL or sQTL by quantification method. GTEx tissues take the eQTL Catalogue tissue label of the
same sample group, so that GTEx v10 and the v8 studies of eQTL Catalogue read alike.

Only the variants flagged as belonging to the 95 % credible set are written, as for the GWAS
credible sets in create_open_targets_files.py. `aaf`, `most_severe` and `gene_most_severe` are
NA until annotate_resource.sh stamps them.

Each parquet file is converted on its own and appended to the per-sub-study files, so only one
file's worth of exploded loci is held at a time.
"""

import glob
import os
import re
import sys

import polars as pl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from create_eqtl_catalogue_files import QUANT_DATA_TYPE
from create_open_targets_files import CHR_MAP, OUTPUT_COLUMNS, _mlog10p, _sci

# projectId -> the prefix its studyIds carry and the sub-study name prefix written here
PROJECTS = {
    "GTEx-v10": ("gtex-v10", "GTEx_v10"),
    "OTAR2057_IBDverse": ("OTAR2057_IBDverse", "IBDverse"),
}

STUDY_COLUMNS = ["studyId", "projectId", "geneId", "traitFromSource", "condition"]

PARQUET_COLUMNS = [
    "studyLocusId",
    "studyId",
    "finemappingMethod",
    "purityMinR2",
    "variantId",
    "pValueMantissa",
    "pValueExponent",
    "locus",
]

PER_STUDY_DIR = "opentargets_qtl_per_study"

EQTL_CATALOGUE_STUDIES = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "metadata", "eqtl_catalogue_studies.tsv"
)


def gtex_tissue_labels(eqtl_catalogue_studies: str) -> dict[str, str]:
    """GTEx sample group -> eQTL Catalogue tissue label, e.g. adipose_visceral -> adipose (visceral)."""
    studies = pl.read_csv(eqtl_catalogue_studies, separator="\t")
    gtex = studies.filter(pl.col("study_label") == "GTEx").select("sample_group", "tissue_label").unique()
    assert gtex["sample_group"].n_unique() == gtex.height, "GTEx sample group with two tissue labels"
    return dict(gtex.iter_rows())


def parse_study_id(study_id: str, study_prefix: str, sample_groups: list[str] | None) -> tuple[str, str]:
    """(quant_method, tissue or cell type) from a QTL studyId.

    GTEx: `gtex-v10_<quant>_<sample group>_<molecular trait>`, where the trait part may start
    with a digit, so the sample group is matched against the known ones rather than split off.
    IBDverse: `OTAR2057_IBDverse_<quant>_<cell type label>_<gene id>`.
    """
    rest = study_id.removeprefix(f"{study_prefix}_")
    quant, _, rest = rest.partition("_")
    if sample_groups is None:
        label, _, gene = rest.rpartition("_")
        if not label or not gene.startswith("ENSG"):
            raise ValueError(f"cannot parse {study_id}")
        return quant, label
    # longest first, in case one sample group is a prefix of another; studyIds are lowercased,
    # so LCL is spelled lcl there
    for group in sample_groups:
        if rest.startswith(f"{group.lower()}_"):
            return quant, group
    raise ValueError(f"no known sample group in {study_id}")


def file_safe(name: str) -> str:
    """Cell type labels carry '+', '(' and ')'; sub-study names are also file names."""
    return re.sub(r"[^A-Za-z0-9_.+-]", "_", name)


def read_studies(data_dir: str, tissue_labels: dict[str, str], gene_names: pl.DataFrame) -> pl.DataFrame:
    """One row per kept study: studyId -> sub-study, data_type, trait, trait_original, cell_type."""
    files = sorted(glob.glob(os.path.join(data_dir, "study_metadata", "*.parquet")))
    if not files:
        sys.exit(f"no parquet files under {data_dir}/study_metadata")
    studies = pl.read_parquet(files, columns=STUDY_COLUMNS).filter(pl.col("projectId").is_in(list(PROJECTS)))

    sample_groups = sorted(tissue_labels, key=len, reverse=True)
    rows = []
    for study_id, project, condition in studies.select("studyId", "projectId", "condition").iter_rows():
        study_prefix, substudy_prefix = PROJECTS[project]
        is_gtex = project == "GTEx-v10"
        quant, group = parse_study_id(study_id, study_prefix, sample_groups if is_gtex else None)
        if quant not in QUANT_DATA_TYPE:
            raise ValueError(f"unknown quantification method {quant} in {study_id}")
        tissue = tissue_labels[group] if is_gtex else group
        rows.append(
            (
                study_id,
                f"{substudy_prefix}_{file_safe(group)}_{quant}",
                QUANT_DATA_TYPE[quant],
                quant,
                f"{tissue}|{condition}".replace(" ", "_"),
            )
        )
    parsed = pl.DataFrame(
        rows, schema=["studyId", "substudy", "data_type", "quant", "cell_type"], orient="row"
    )
    # two cell type labels that differ only in the characters file_safe replaces would share a file
    assert (
        parsed.select("substudy", "cell_type").unique().height == parsed["substudy"].n_unique()
    ), "sub-study name shared by two cell types"

    return (
        studies.join(parsed, on="studyId")
        .join(gene_names, left_on="geneId", right_on="gene_key", how="left")
        .select(
            "studyId",
            "substudy",
            "data_type",
            "cell_type",
            # symbol-less genes keep their gene id, as in the eQTL Catalogue munge
            pl.coalesce("gene_name", "geneId").alias("trait"),
            pl.concat_str(["traitFromSource", "quant"], separator="|").alias("trait_original"),
        )
    )


def read_gene_names(path: str) -> pl.DataFrame:
    return (
        pl.read_csv(path, separator="\t", schema_overrides={"chromosome": pl.Utf8})
        .select(pl.col("phenotype_id").str.replace(r"\..*$", "").alias("gene_key"), "gene_name")
        .unique(subset="gene_key")
    )


def convert_pq_to_df(parquet_path: str, studies: pl.DataFrame) -> pl.DataFrame:
    """Read one credible set parquet file into one row per 95 % credible set variant of a kept study."""
    df = (
        pl.read_parquet(parquet_path, columns=PARQUET_COLUMNS)
        .filter(pl.col("studyId").is_in(studies["studyId"].implode()))
        .filter(
            pl.col("finemappingMethod").str.to_lowercase().str.contains("susie", literal=True)
            & pl.col("locus").is_not_null()
        )
        .rename(
            {
                "variantId": "lead_variant_id",
                "pValueMantissa": "cs_p_mantissa",
                "pValueExponent": "cs_p_exponent",
            }
        )
        .with_row_index("_cs_row")
        .explode("locus")
        .unnest("locus")
        .filter(pl.col("is95CredibleSet"))
        .rename({"standardError": "se"})
        .join(studies, on="studyId")
    )
    df = df.with_columns(pl.col("variantId").str.split("_").alias("_cpra"))
    return df.with_columns(
        pl.col("_cpra").list.get(0).replace(CHR_MAP).cast(pl.Int32, strict=False).alias("chr"),
        pl.col("_cpra").list.get(1).cast(pl.Int64, strict=False).alias("pos"),
        pl.col("_cpra").list.get(2).alias("ref"),
        pl.col("_cpra").list.get(3).alias("alt"),
        # every QTL credible set variant carries its own p-value in this release; the fallback
        # is kept so that a release where one does not still gets a lead p-value
        pl.when(pl.col("pValueMantissa").is_not_null() & pl.col("pValueExponent").is_not_null())
        .then(_mlog10p(pl.col("pValueMantissa"), pl.col("pValueExponent")))
        .when(pl.col("variantId") == pl.col("lead_variant_id"))
        .then(_mlog10p(pl.col("cs_p_mantissa"), pl.col("cs_p_exponent")))
        .alias("mlog10p"),
        _sci("beta"),
        _sci("se"),
        pl.col("posteriorProbability").round(4).alias("pip"),
        pl.col("studyLocusId").alias("cs_id"),
        pl.len().over("_cs_row").cast(pl.Int32).alias("cs_size"),
        pl.col("purityMinR2").round(4).alias("cs_min_r2"),
    )


def append_per_substudy(df: pl.DataFrame, dataset: str, out_dir: str) -> None:
    na = pl.lit(None, dtype=pl.String)
    out = df.with_columns(
        pl.lit(dataset).alias("dataset"),
        na.alias("aaf"),
        na.alias("most_severe"),
        na.alias("gene_most_severe"),
    )
    for part in out.partition_by("substudy"):
        path = os.path.join(out_dir, f"{part['substudy'][0]}.SUSIE.munged.tsv")
        new = not os.path.exists(path)
        with open(path, "ab") as f:
            part.select(OUTPUT_COLUMNS).write_csv(f, separator="\t", null_value="NA", include_header=new)


def sort_substudy(path: str, schema: dict) -> None:
    pl.read_csv(path, separator="\t", null_values=["NA"], schema=schema).sort(
        "chr", "pos", "ref", "alt", "trait", "trait_original"
    ).write_csv(path, separator="\t", null_value="NA")


def main(dataset: str, data_dir: str) -> None:
    out_dir = os.path.join(data_dir, PER_STUDY_DIR)
    if glob.glob(os.path.join(out_dir, "*.SUSIE.munged.tsv")):
        # the files are appended to, so a rerun over a partial run would duplicate rows
        sys.exit(f"{out_dir} already has per-study files; remove them first")
    os.makedirs(out_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(data_dir, "credible_set", "*.parquet")))
    if not files:
        sys.exit(f"no parquet files under {data_dir}/credible_set")

    studies = read_studies(
        data_dir,
        gtex_tissue_labels(EQTL_CATALOGUE_STUDIES),
        read_gene_names(os.path.join(data_dir, "gene_counts_Ensembl_105_phenotype_metadata.tsv.gz")),
    )
    print(
        f"{studies.height} studies in {studies['substudy'].n_unique()} sub-studies, "
        f"{studies.filter(pl.col('trait').str.starts_with('ENSG')).height} without a gene name"
    )

    n_rows = 0
    for i, path in enumerate(files, 1):
        df = convert_pq_to_df(path, studies)
        append_per_substudy(df, dataset, out_dir)
        n_rows += df.height
        print(f"[{i}/{len(files)}] {os.path.basename(path)}: {df.height} rows", flush=True)
        del df

    per_study = sorted(glob.glob(os.path.join(out_dir, "*.SUSIE.munged.tsv")))
    # one schema for every file, as create_open_targets_per_study_files.py does: inferred per
    # file, an all-NA column would be typed differently from one file to the next
    schema = dict(pl.read_csv(per_study[0], separator="\t", null_values=["NA"], n_rows=1000).schema)
    schema.update({c: pl.String for c in ("aaf", "most_severe", "gene_most_severe", "beta", "se")})
    for i, path in enumerate(per_study, 1):
        sort_substudy(path, schema)
        print(f"[{i}/{len(per_study)}] sorted {os.path.basename(path)}", flush=True)

    print(f"wrote {n_rows} credible set variants to {len(per_study)} files under {out_dir}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: python create_open_targets_qtl_files.py <dataset_name> <data_dir>")
    main(sys.argv[1], sys.argv[2])
