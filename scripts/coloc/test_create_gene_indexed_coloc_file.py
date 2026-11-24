"""How a coloc pair's QTL side is resolved to a gene, and which pairs the gene-indexed copy holds."""

import polars as pl

from create_gene_indexed_coloc_file import GENE_COLS, gene_indexed

GENES = pl.DataFrame(
    {
        "ensg": ["ENSG01", "ENSG02", "ENSG03", "ENSG04"],
        "gene_name": ["IL23R", "PARK7", "NADK", "H6PD"],
        "gene_chr": [1, 1, 1, 1],
        "gene_start": [100, 200, 300, 400],
        "gene_end": [150, 250, 350, 450],
    },
    schema_overrides={"gene_chr": pl.Int8, "gene_start": pl.Int32, "gene_end": pl.Int32},
)
TRAIT_GENE_MAP = pl.DataFrame({"trait_original": ["seq.1.2"], "mapped_symbol": ["PARK7"]})


def _pairs(*sides2, type1="GWAS") -> pl.DataFrame:
    """Pairs with a GWAS on side 1 and (data_type, trait, trait_original) on side 2."""
    return pl.DataFrame(
        {
            "#dataset1": ["FinnGen_R14"] * len(sides2),
            "dataset2": ["D"] * len(sides2),
            "data_type1": [type1] * len(sides2),
            "data_type2": [s[0] for s in sides2],
            "trait1": ["T1D"] * len(sides2),
            "trait1_original": ["T1D"] * len(sides2),
            "trait2": [s[1] for s in sides2],
            "trait2_original": [s[2] for s in sides2],
            "chr": ["12"] * len(sides2),
        }
    )


def _genes_of(coloc: pl.DataFrame, trait_gene_map=TRAIT_GENE_MAP) -> list:
    return gene_indexed(coloc, GENES, trait_gene_map)["qtl_gene"].to_list()


def test_each_resolution_rule():
    coloc = _pairs(
        ("eQTL", "anything", "ENSG01.7_1_100_150|exon"),  # Ensembl id in trait_original
        ("eQTL", "ENSG02", "ENST0001|tx"),  # Ensembl id as the trait
        ("pQTL", "seq.1.2", "seq.1.2"),  # reagent id through the map
        ("pQTL", "PARK7", "PARK7"),  # symbol
        ("pQTL", "NADK_cross_batch_normalised_4144", "NADK_cross_batch_normalised_4144"),
        ("pQTL", "H6PD.7161.25.3..1", "H6PD.7161.25.3..1"),  # SYMBOL.<aptamer>
        ("sQTL", "IL23R", "1:100:120:clu_1_+|leafcutter"),
    )
    assert _genes_of(coloc) == ["IL23R", "PARK7", "PARK7", "PARK7", "NADK", "H6PD", "IL23R"]


def test_ensembl_id_wins_over_a_symbol_that_also_matches():
    assert _genes_of(_pairs(("eQTL", "PARK7", "ENSG01|ge"))) == ["IL23R"]


def test_unresolved_trait_keeps_null_gene_columns():
    out = gene_indexed(_pairs(("pQTL", "NTproBNP", "NTproBNP")), GENES, None)
    assert out.select(GENE_COLS).drop("qtl_side").row(0) == (None, None, None, None)


def test_reagent_id_needs_the_map():
    assert _genes_of(_pairs(("pQTL", "seq.1.2", "seq.1.2")), trait_gene_map=None) == [None]


def test_pairs_without_a_gene_qtl_side_are_left_out():
    coloc = _pairs(("GWAS", "T2D", "T2D"), ("metaboQTL", "ALA", "ALA"), ("caQTL", "chr1-1-2", "chr1-1-2"))
    assert gene_indexed(coloc, GENES, None).height == 0


def test_a_pair_of_two_gene_qtls_is_emitted_under_each_side():
    coloc = _pairs(("pQTL", "PARK7", "PARK7"), type1="eQTL").with_columns(
        pl.lit("IL23R").alias("trait1"), pl.lit("ENSG01").alias("trait1_original")
    )
    out = gene_indexed(coloc, GENES, None).sort("qtl_side")
    assert out["qtl_side"].to_list() == [1, 2]
    assert out["qtl_gene"].to_list() == ["IL23R", "PARK7"]
    assert out["gene_start"].to_list() == [100, 200]


def test_original_columns_pass_through_unchanged_and_first():
    coloc = _pairs(("pQTL", "PARK7", "PARK7"))
    out = gene_indexed(coloc, GENES, None)
    assert out.columns == coloc.columns + GENE_COLS
    assert out.select(coloc.columns).equals(coloc)
