"""Tests for the Open Targets QTL munge: the studyId parsing it relies on, and one end-to-end run
over a release cut down to a few rows."""

import gzip

import polars as pl
import pytest

import create_open_targets_qtl_files as otq
from create_open_targets_qtl_files import parse_study_id

GTEX_GROUPS = sorted(["adipose_subcutaneous", "adipose_visceral", "brain_cortex", "LCL"], key=len, reverse=True)


@pytest.mark.parametrize(
    "study_id, expected",
    [
        ("gtex-v10_ge_adipose_subcutaneous_ensg00000003147", ("ge", "adipose_subcutaneous")),
        ("gtex-v10_tx_adipose_visceral_enst00000304987", ("tx", "adipose_visceral")),
        (
            "gtex-v10_exon_brain_cortex_ensg00000115750_17_2_9844342_9844477",
            ("exon", "brain_cortex"),
        ),
        # studyIds are lowercased, the sample group keeps the eQTL Catalogue spelling
        (
            "gtex-v10_txrev_lcl_ensg00000089280_grp_2_contained_enst00000487974",
            ("txrev", "LCL"),
        ),
        # leafcutter traits start with the chromosome, so the sample group cannot be split off
        # at the first digit
        ("gtex-v10_leafcutter_lcl_11_71480879_71481923_clu_60703__", ("leafcutter", "LCL")),
        (
            "gtex-v10_majiq_adipose_subcutaneous_ensg00000010610_t_6814054-6814463_6800471-6814142",
            ("majiq", "adipose_subcutaneous"),
        ),
    ],
)
def test_gtex_study_ids(study_id, expected):
    assert parse_study_id(study_id, "gtex-v10", GTEX_GROUPS) == expected


def test_gtex_study_id_with_unknown_sample_group_fails():
    with pytest.raises(ValueError):
        parse_study_id("gtex-v10_ge_kidney_medulla_ensg00000003147", "gtex-v10", GTEX_GROUPS)


@pytest.mark.parametrize(
    "study_id, expected",
    [
        ("OTAR2057_IBDverse_ge_CD4+_CRM_ENSG00000107679", ("ge", "CD4+_CRM")),
        # a label ending in an underscore keeps it
        ("OTAR2057_IBDverse_ge_Memory_B_CD27Lo__ENSG00000107679", ("ge", "Memory_B_CD27Lo_")),
        ("OTAR2057_IBDverse_ge_Colon_Goblet_ATOH1_(MUC2_Lo)_ENSG00000107679", ("ge", "Colon_Goblet_ATOH1_(MUC2_Lo)")),
    ],
)
def test_ibdverse_study_ids(study_id, expected):
    assert parse_study_id(study_id, "OTAR2057_IBDverse", None) == expected


def _locus(variant, pip, is95=True, p=(1.0, -10), beta=0.5, se=0.1):
    return {
        "is95CredibleSet": is95,
        "is99CredibleSet": True,
        "logBF": 1.0,
        "posteriorProbability": pip,
        "variantId": variant,
        "pValueMantissa": p[0],
        "pValueExponent": p[1],
        "beta": beta,
        "standardError": se,
        "r2Overall": 1.0,
    }


@pytest.fixture
def release(tmp_path, monkeypatch):
    (tmp_path / "study_metadata").mkdir()
    (tmp_path / "credible_set").mkdir()
    pl.DataFrame(
        {
            "studyId": [
                "gtex-v10_ge_adipose_subcutaneous_ensg00000000001",
                "gtex-v10_leafcutter_adipose_subcutaneous_1_100_200_clu_1__",
                "OTAR2057_IBDverse_ge_CD4+_CRM_ENSG00000000002",
                "QTD000001_ENSG00000000001",
                "GCST1",
            ],
            "projectId": ["GTEx-v10", "GTEx-v10", "OTAR2057_IBDverse", "Alasoo_2018", "GCST"],
            "geneId": ["ENSG00000000001", "ENSG00000000001", "ENSG00000000002", "ENSG00000000001", None],
            "traitFromSource": ["ENSG00000000001", "1:100:200:clu_1_+", "ENSG00000000002", "ENSG00000000001", "Height"],
            "condition": ["naive", "naive", "naive", "naive", None],
        }
    ).write_parquet(tmp_path / "study_metadata" / "part-0.parquet")

    def cs(study, locus_id, lead, loci):
        return {
            "studyLocusId": locus_id,
            "studyId": study,
            "finemappingMethod": "SuSie",
            "purityMinR2": 0.9,
            "variantId": lead,
            "pValueMantissa": 1.0,
            "pValueExponent": -10,
            "locus": loci,
        }

    rows = [
        cs(
            "gtex-v10_ge_adipose_subcutaneous_ensg00000000001",
            "L1",
            "1_300_A_G",
            [_locus("1_300_A_G", 0.7), _locus("1_100_C_T", 0.25), _locus("1_500_G_A", 0.04, is95=False)],
        ),
        cs("gtex-v10_leafcutter_adipose_subcutaneous_1_100_200_clu_1__", "L2", "X_50_A_C", [_locus("X_50_A_C", 0.99)]),
        cs("OTAR2057_IBDverse_ge_CD4+_CRM_ENSG00000000002", "L3", "2_10_T_C", [_locus("2_10_T_C", 0.96, p=(5.0, -8))]),
        cs("QTD000001_ENSG00000000001", "L4", "1_300_A_G", [_locus("1_300_A_G", 0.9)]),
        cs("GCST1", "L5", "1_300_A_G", [_locus("1_300_A_G", 0.9)]),
    ]
    pl.DataFrame(rows).write_parquet(tmp_path / "credible_set" / "part-0.parquet")

    with gzip.open(tmp_path / "gene_counts_Ensembl_105_phenotype_metadata.tsv.gz", "wt") as f:
        f.write("phenotype_id\tgene_name\tchromosome\nENSG00000000001.5\tGENE1\t1\n")

    studies_tsv = tmp_path / "eqtl_catalogue_studies.tsv"
    pl.DataFrame(
        {
            "study_label": ["GTEx", "GTEx", "Alasoo_2018"],
            "sample_group": ["adipose_subcutaneous", "adipose_subcutaneous", "macrophage_naive"],
            "tissue_label": ["adipose", "adipose", "macrophage"],
        }
    ).write_csv(studies_tsv, separator="\t")
    monkeypatch.setattr(otq, "EQTL_CATALOGUE_STUDIES", str(studies_tsv))
    return tmp_path


def _read(path):
    return pl.read_csv(path, separator="\t", null_values=["NA"], infer_schema_length=0)


def test_end_to_end(release):
    otq.main("Open_Targets_QTL_26.09", str(release))
    out = release / otq.PER_STUDY_DIR

    files = sorted(p.name for p in out.glob("*.SUSIE.munged.tsv"))
    # only the allow-listed projects, one file per project x tissue x quantification method
    assert files == [
        "GTEx_v10_adipose_subcutaneous_ge.SUSIE.munged.tsv",
        "GTEx_v10_adipose_subcutaneous_leafcutter.SUSIE.munged.tsv",
        "IBDverse_CD4+_CRM_ge.SUSIE.munged.tsv",
    ]

    ge = _read(out / files[0])
    assert ge.columns == otq.OUTPUT_COLUMNS
    # sorted by position, and the variant outside the 95 % set is dropped
    assert ge["pos"].to_list() == ["100", "300"]
    row = ge.row(1, named=True)
    assert row["dataset"] == "Open_Targets_QTL_26.09"
    assert row["data_type"] == "eQTL"
    assert row["trait"] == "GENE1"
    assert row["trait_original"] == "ENSG00000000001|ge"
    assert row["cell_type"] == "adipose|naive"
    assert row["cs_id"] == "L1"
    assert row["cs_size"] == "2"
    assert row["mlog10p"] == "10.0"
    assert row["beta"] == "5.000e-01"
    assert row["se"] == "1.000e-01"
    assert row["aaf"] is None and row["most_severe"] is None

    sqtl = _read(out / files[1]).row(0, named=True)
    assert sqtl["data_type"] == "sQTL"
    assert sqtl["trait_original"] == "1:100:200:clu_1_+|leafcutter"
    assert sqtl["chr"] == "23"

    ibd = _read(out / files[2]).row(0, named=True)
    # a gene without a name keeps its gene id
    assert ibd["trait"] == "ENSG00000000002"
    assert ibd["cell_type"] == "CD4+_CRM|naive"

    # no stats: nothing reads them for this dataset
    assert sorted(p.name for p in out.iterdir()) == files


def test_rerun_over_existing_output_refuses(release):
    otq.main("Open_Targets_QTL_26.09", str(release))
    with pytest.raises(SystemExit):
        otq.main("Open_Targets_QTL_26.09", str(release))
