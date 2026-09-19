"""Unit tests for the parts of `munge_brava.py` whose wrong answer still produces a file that
loads and indexes: which test rows are kept, how the mask becomes an annotation, and what a
quantitative phenotype's sample count is called."""

import polars as pl

from munge_brava import build_output, keep_burden_ivw

GENCODE = pl.DataFrame(
    {
        "gene_id_base": ["ENSG00000000001"],
        "gene": ["GENE1"],
        "gene_chr": [1],
        "gene_start_pos": [100],
        "gene_end_pos": [200],
    }
)


def _source() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "Region": ["ENSG00000000001"] * 4,
            "Group": ["pLoF;damaging_missense_or_protein_altering", "pLoF", "pLoF", "pLoF"],
            "max_MAF": ["1e-04", "0.001", "0.001", "1e-04"],
            "Pvalue": [1e-10, 0.5, 0.5, 0.5],
            "type": ["Inverse variance weighted", "Inverse variance weighted",
                     "Stouffer", "Inverse variance weighted"],
            "BETA_Burden": [0.1, 0.2, 0.3, 0.4],
            "SE_Burden": [0.01, 0.02, 0.03, 0.04],
            "class": ["Burden", "Burden", "Burden", "SKAT-O"],
        }
    )


def test_only_the_inverse_variance_burden_rows_are_kept():
    kept = keep_burden_ivw(_source())
    assert kept.height == 2
    assert kept["Group"].to_list() == ["pLoF;damaging_missense_or_protein_altering", "pLoF"]


def test_annotation_pastes_the_mask_and_max_maf_as_the_source_wrote_them():
    out, _, _ = build_output(keep_burden_ivw(_source()), GENCODE,
                             {"phenocode": "AFib", "phenostring": "Atrial fibrillation",
                              "num_cases": 10, "num_controls": 90})
    assert out["annotation"].to_list() == [
        "pLoF|damaging_missense_or_protein_altering|MAF<1e-04",
        "pLoF|MAF<0.001",
    ]
    assert out["n_cases"].to_list() == [10, 10]
    assert out["n_controls"].to_list() == [90, 90]


def test_a_quantitative_phenotype_puts_its_sample_count_in_n_cases():
    out, _, _ = build_output(keep_burden_ivw(_source()), GENCODE,
                             {"phenocode": "LDLC", "phenostring": "LDLC", "num_samples": 1000})
    # the genebass convention for a trait with no cases: n_cases is N and n_controls is missing
    assert out["n_cases"].to_list() == [1000, 1000]
    assert out["n_controls"].to_list() == [None, None]


def test_pvalue_zero_recovers_mlog10p_from_beta_se():
    # LDLC APOB pLoF|MAF<1e-04, the flagship cell the recovery path exists for
    source = pl.DataFrame(
        {
            "Region": ["ENSG00000000001"] * 2,
            "Group": ["pLoF", "pLoF"],
            "max_MAF": ["1e-04", "1e-04"],
            "Pvalue": [0.0, 1e-10],
            "type": ["Inverse variance weighted"] * 2,
            "BETA_Burden": [-0.0515217619899695, 0.1],
            "SE_Burden": [0.0010842697450283, 0.01],
            "class": ["Burden", "Burden"],
        }
    )
    out, _, _ = build_output(keep_burden_ivw(source), GENCODE,
                             {"phenocode": "LDLC", "phenostring": "LDLC", "num_samples": 1000})
    mlog10p = out["mlog10p_burden"].to_list()
    assert mlog10p[0] == 492.0742
    assert mlog10p[1] == 10.0  # p > 0 still goes through -log10(p), not the recovery path
