"""
Build the phenotype metadata JSON for the EstBB-UKBB NMR metabolic trait fine-mapping
(Tambets et al. 2026), mapping each of the 249 Nightingale biomarker codes the results carry
to the trait name its authors registered in the GWAS Catalog.

    python3 scripts/nmr_meta_phenotypes.py \
      --input UKBB_EUR_fine_mapping_with_meta_EUR_lead_variants.parquet \
      --output configs/nmr_meta_pheno.json

Output is the pheweb-shaped list the results platform's `quantitative_pheweb` metadata
harmonizer reads: `phenocode`, `phenostring`, `num_samples`.

THE NAMES ARE NOT WRITTEN HERE, THEY ARE MATCHED. Typing 249 biomarker descriptions by hand
would put a second, unverifiable vocabulary next to the authors' own, and the errors would be
invisible -- "Cholesteryl esters in large HDL" and "Cholesteryl esters in very large HDL" differ
by one word. Instead the GWAS Catalog's 2,241 study records for this publication (249 traits x 9
ancestry groups) supply the strings, and the mapping onto the codes is derived:

  - the lipoprotein codes are compositional, `[size_]class_measure[_pct]`, and the description
    follows from the parts. That grammar is `_compose` below;
  - the rest -- amino acids, fatty acids, ratios, the whole-serum totals -- are not
    compositional and are listed in SPECIAL;
  - every generated description is then matched, case- and punctuation-insensitively, against
    the catalogue, and the script FAILS unless the result is a bijection: 249 codes onto 249
    distinct catalogue traits, none left over on either side.

So a wrong SPECIAL entry or a gap in the grammar cannot pass silently -- it takes a catalogue
name some other code needed, or leaves one unclaimed, and either breaks the bijection. What
reaches the output is always the catalogue's string, never the generated one.

The three class-level cholesterol codes are in SPECIAL rather than the grammar because the
catalogue names them "HDL cholesterol levels" while the same measure one size down is
"Cholesterol in large HDL"; `IDL_C` does follow the grammar. That inconsistency is the
catalogue's, and pinning it here is what keeps the bijection honest.
"""

import argparse
import json
import re
import sys
import urllib.request
from collections import Counter

GWAS_CATALOG_PUBMED_ID = "42162431"
GWAS_CATALOG_URL = (
    "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByPublicationIdPubmedId"
    "?pubmedId={pubmed_id}&size=500&page={page}"
)

SIZES = {
    "XXL": "chylomicrons and extremely large VLDL",
    "XL": "very large",
    "L": "large",
    "M": "medium",
    "S": "small",
    "XS": "very small",
}

CLASSES = ("VLDL", "LDL", "IDL", "HDL")

MEASURES = {
    "C": "Cholesterol",
    "CE": "Cholesteryl esters",
    "FC": "Free cholesterol",
    "PL": "Phospholipids",
    "TG": "Triglycerides",
    "L": "Total lipids",
}

SPECIAL = {
    "Acetate": "Acetate levels",
    "Acetoacetate": "Acetoacetate levels",
    "Acetone": "Acetone levels",
    "Ala": "Alanine levels",
    "Albumin": "Albumin levels",
    "ApoA1": "Apolipoprotein A1 levels",
    "ApoB": "Apolipoprotein B levels",
    "ApoB_by_ApoA1": "Ratio of apolipoprotein B to apolipoprotein A1 levels",
    "Cholines": "Total cholines levels",
    "Citrate": "Citrate levels",
    "Clinical_LDL_C": "Clinical LDL cholesterol levels",
    "Creatinine": "Creatinine levels",
    "DHA": "Docosahexaenoic acid levels",
    "DHA_pct": "Ratio of docosahexaenoic acid to total fatty acid levels",
    "Gln": "Glutamine levels",
    "Glucose": "Glucose levels",
    "Gly": "Glycine levels",
    "GlycA": "Glycoprotein acetyls levels",
    "HDL_C": "HDL cholesterol levels",
    "His": "Histidine levels",
    "Ile": "Isoleucine levels",
    "LA": "Linoleic acid levels",
    "LA_pct": "Ratio of linoleic acid to total fatty acids",
    "LDL_C": "LDL cholesterol levels",
    "Lactate": "Lactate levels",
    "Leu": "Leucine levels",
    "MUFA": "Monounsaturated fatty acid levels",
    "MUFA_pct": "Ratio of monounsaturated fatty acids to total fatty acids",
    "Omega_3": "Omega-3 fatty acids levels",
    "Omega_3_pct": "Ratio of omega-3 fatty acids to total fatty acids",
    "Omega_6": "Omega-6 fatty acids levels",
    "Omega_6_by_Omega_3": "Ratio of omega-6 fatty acids to omega-3 fatty acids",
    "Omega_6_pct": "Ratio of omega-6 fatty acids to total fatty acids",
    "PUFA": "Polyunsaturated fatty acid levels",
    "PUFA_by_MUFA": "Ratio of polyunsaturated fatty acids to monounsaturated fatty acids",
    "PUFA_pct": "Ratio of polyunsaturated fatty acids to total fatty acids",
    "Phe": "Phenylalanine levels",
    "Phosphatidylc": "Phosphatidylcholine levels",
    "Phosphoglyc": "Phosphoglycerides levels",
    "Pyruvate": "Pyruvate levels",
    "Remnant_C": "Remnant cholesterol (non-HDL, non-LDL -cholesterol)",
    "SFA": "Saturated fatty acids levels",
    "SFA_pct": "Ratio of saturated fatty acids to total fatty acids",
    "Sphingomyelins": "Sphingomyelins levels",
    "TG_by_PG": "Ratio of triglycerides to phosphoglycerides",
    "Total_BCAA": (
        "Total concentration of branched-chain amino acids (leucine + isoleucine + valine)"
    ),
    "Total_C": "Total cholesterol levels",
    "Total_CE": "Total esterified cholesterol levels",
    "Total_FA": "Total fatty acid levels",
    "Total_FC": "Total free cholesterol levels",
    "Total_L": "Total lipid levels in lipoprotein particles",
    "Total_P": "Total concentration of lipoprotein particles",
    "Total_PL": "Total phospholipids in lipoprotein particles",
    "Total_TG": "Total triglycerides levels",
    "Tyr": "Tyrosine levels",
    "Unsaturation": "Degree of unsaturation",
    "VLDL_C": "VLDL cholesterol levels",
    "Val": "Valine levels",
    "bOHbutyrate": "3-Hydroxybutyrate levels",
    "non_HDL_C": "Total cholesterol minus HDL-C levels",
}

LIPOPROTEIN_PATTERN = re.compile(
    r"^(?:(?P<size>XXL|XL|XS|L|M|S)_)?(?P<klass>VLDL|LDL|IDL|HDL)_"
    r"(?P<measure>CE|FC|PL|TG|C|L|P|size)(?P<pct>_pct)?$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Zenodo fine-mapping parquet")
    parser.add_argument("--output", required=True, help="Path to write the phenotype JSON to")
    parser.add_argument(
        "--studies",
        nargs="*",
        help="Cached GWAS Catalog API responses for the publication; fetched live if omitted, "
        "which takes several minutes per page",
    )
    parser.add_argument(
        "--n-samples",
        type=int,
        default=413897,
        help="UKBB_EUR sample size the fine-mapping was run on; every trait shares it",
    )
    return parser.parse_args()


def normalize(name: str) -> str:
    """Fold the differences that are not information: case, punctuation, repeated spaces.

    The catalogue is internally inconsistent about capitalisation of the very strings this
    matches on -- "Total Lipids in Medium LDL" next to "Total lipids in small LDL" -- so an
    exact comparison would fail on the catalogue's own typography.
    """
    return re.sub(r"[^a-z0-9]+", " ", name.lower()).strip()


def _compose(code: str) -> str | None:
    """Describe a compositional lipoprotein code, or return None if it is not one."""
    match = LIPOPROTEIN_PATTERN.match(code)
    if match is None:
        return None

    size, klass, measure, pct = match.group("size", "klass", "measure", "pct")
    # the XXL label already names the lipoprotein class, so it replaces the pair
    subclass = SIZES["XXL"] if size == "XXL" else f"{SIZES[size]} {klass}" if size else klass

    if measure == "size":
        return f"Average diameter for {klass} particles"
    if measure == "P":
        return f"Concentration of {subclass} particles"
    if pct:
        return f"{MEASURES[measure]} to total lipids ratio in {subclass}"
    return f"{MEASURES[measure]} in {subclass}"


def describe(code: str) -> str | None:
    return SPECIAL.get(code) or _compose(code)


def read_trait_codes(path: str) -> list[str]:
    import polars as pl

    return sorted(pl.read_parquet(path, columns=["molecular_trait_id"])["molecular_trait_id"].unique())


def fetch_catalog_traits(paths: list[str] | None) -> set[str]:
    """Collect the distinct trait names the GWAS Catalog holds for this publication."""
    pages: list[dict] = []
    if paths:
        pages = [json.load(open(path)) for path in paths]
    else:
        page = 0
        while True:
            url = GWAS_CATALOG_URL.format(pubmed_id=GWAS_CATALOG_PUBMED_ID, page=page)
            print(f"  fetching page {page}...")
            with urllib.request.urlopen(url, timeout=900) as response:
                body = json.load(response)
            pages.append(body)
            page += 1
            if page >= body.get("page", {}).get("totalPages", 0):
                break

    traits = {
        (study.get("diseaseTrait") or {}).get("trait")
        for page_body in pages
        for study in page_body["_embedded"]["studies"]
    }
    traits.discard(None)
    return traits


def match(codes: list[str], catalog: set[str]) -> dict[str, str]:
    """Pair every code with a distinct catalogue trait, or fail saying which ones did not.

    The bijection is the whole check: nothing here verifies a description against the biology,
    only that the derivation accounts for exactly the traits the authors registered.
    """
    by_normalized: dict[str, str] = {normalize(trait): trait for trait in catalog}
    if len(by_normalized) != len(catalog):
        duplicates = [t for t, n in Counter(normalize(t) for t in catalog).items() if n > 1]
        sys.exit(f"catalogue traits collide once normalized: {duplicates}")

    mapping: dict[str, str] = {}
    unmatched: list[str] = []
    for code in codes:
        described = describe(code)
        trait = by_normalized.get(normalize(described)) if described else None
        if trait is None:
            unmatched.append(f"{code} -> {described!r}")
        else:
            mapping[code] = trait

    claimed = Counter(mapping.values())
    collisions = [trait for trait, n in claimed.items() if n > 1]
    leftover = sorted(set(catalog) - set(mapping.values()))

    if unmatched or collisions or leftover:
        for line in unmatched:
            print(f"  no catalogue trait for {line}", file=sys.stderr)
        for trait in collisions:
            print(f"  several codes claim {trait!r}", file=sys.stderr)
        for trait in leftover:
            print(f"  no code claims {trait!r}", file=sys.stderr)
        sys.exit(
            f"not a bijection: {len(unmatched)} unmatched codes, {len(collisions)} collisions, "
            f"{len(leftover)} unclaimed catalogue traits"
        )

    return mapping


def main() -> None:
    args = parse_args()

    print(f"Reading trait codes from {args.input}...")
    codes = read_trait_codes(args.input)
    print(f"  {len(codes)} traits")

    print("Reading GWAS Catalog studies...")
    catalog = fetch_catalog_traits(args.studies)
    print(f"  {len(catalog)} distinct trait names")

    mapping = match(codes, catalog)
    print(f"  matched all {len(mapping)} codes onto distinct catalogue traits")

    with open(args.output, "w") as out:
        json.dump(
            [
                {
                    "phenocode": code,
                    "phenostring": mapping[code],
                    "num_samples": args.n_samples,
                }
                for code in codes
            ],
            out,
            indent=4,
        )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
