#!/usr/bin/env python3
"""Build the per-gene evidence bundle an AI gene summary is written from.

Every gene gets the same sections in the same order, cut down by fixed rules, so a summary
can be checked against its bundle and regenerated only when the bundle's hash changes.
An empty section is kept and says why it is empty: "no rows" and "source unavailable" are
different statements and the summary must not confuse them.

The data comes from the routes the gene page's tabs read, so the summary and the tabs
cannot disagree.
"""

import argparse
import csv
import hashlib
import io
import json
import math
import statistics
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BUNDLE_VERSION = 3

BURDEN_MLOG10P_MIN = 4  # the gene page's burden cut
# exome-wide significance for a single-variant exome association, p < 1e-6
EXOME_MLOG10P_MIN = 6
MLOG10P_FLOOR = 324
GWS_MLOG10P = 7.30103
CODING_PIP_MIN = 0.05
COLOC_PP_H4_MIN = 0.8
# a QTL signal this close to its gene counts as cis
CIS_WINDOW = 1_000_000
CELL_MIN_CELLS = 500
CELL_MIN_DONORS = 3

# exome variant rows of these datasets carry no control count, and their traits are all binary
BINARY_EXOME_DATASETS = {"IBD_exome"}

CAPS = {
    "burden": 20,
    "burden_per_dataset": 10,
    "exome_variants": 12,
    "exome_traits_per_variant": 3,
    "coding_variants": 12,
    "traits_per_coding_variant": 5,
    "gwas_sets": 20,
    "qtl_sets_per_type": 6,
    "coloc": 15,
    "coloc_trans": 8,
    "coloc_trans_loci": 4,
    "gtex_tissues": 6,
    "hpa": 15,
    "cell_types": 10,
    "diseases": 12,
}

# mirror of bff/coding.ts
CODING = {
    "missense", "frameshift", "inframe_insertion", "inframe_deletion", "transcript_ablation",
    "stop_gained", "stop_lost", "start_lost", "splice_acceptor", "splice_donor",
    "incomplete_terminal_codon", "protein_altering",
}

GENCC_RANK = ["Definitive", "Strong", "Moderate", "Supportive", "Limited",
              "Disputed Evidence", "Refuted Evidence", "No Known Disease Relationship"]


class Unavailable(Exception):
    pass


def fetch(base: str, path: str, tries: int = 3, timeout: int = 180) -> str:
    err = None
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(f"{base}{path}", timeout=timeout) as r:
                return r.read().decode()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return ""
            err = e
        except (urllib.error.URLError, TimeoutError) as e:
            err = e
        time.sleep(2 * (attempt + 1))
    raise Unavailable(f"{path}: {err}")


def tsv(text: str) -> list[dict]:
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    return [{k: (None if v in ("NA", "", None) else v) for k, v in r.items()} for r in rows]


def num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def sig(v, digits: int = 3):
    if v is None:
        return None
    if v == 0:
        return 0
    out = round(v, digits - 1 - math.floor(math.log10(abs(v))))
    return int(out) if float(out).is_integer() and abs(out) >= 100 else out


def pval(mlog10p):
    """A p-value string from -log10(p); a string because the largest ones underflow a float."""
    if mlog10p is None:
        return None
    exp = math.ceil(mlog10p)
    mant = 10 ** (exp - mlog10p)
    if round(mant, 1) >= 10:
        mant, exp = mant / 10, exp - 1
    return f"{mant:.1f}e-{exp}"


def stats(mlog10p, beta, binary: bool) -> dict:
    out = {"p": pval(mlog10p), "mlog10p": sig(mlog10p, 4), "beta": sig(beta)}
    if binary and beta is not None and abs(beta) < 20:
        out["odds_ratio"] = sig(math.exp(beta))
    return out


def distance_kb(gene: dict, variant: str):
    """Distance from a chr:pos:ref:alt variant to the gene in kb, 0 inside it, None on
    another chromosome."""
    chrom, pos = variant.split(":")[:2]
    if chrom != str(gene["chr"]):
        return None
    pos = int(pos)
    gap = max(gene["start"] - pos, pos - gene["end"], 0)
    return sig(gap / 1000) if gap else 0


def location(gene: dict, variant: str) -> dict:
    d = distance_kb(gene, variant)
    if d is None:
        return {"location": "another chromosome (chr" + variant.split(":")[0] + ")"}
    return {"location": "inside the gene" if d == 0 else "outside the gene", "distance_to_gene_kb": d}


def q(s: str) -> str:
    return urllib.parse.quote(s, safe="")


# ---------------------------------------------------------------- sections


def gene_section(base: str, symbol: str) -> dict:
    hits = json.loads(fetch(base, f"/search?q={q(symbol)}"))
    hit = next((h for h in hits if h.get("type") == "gene" and h.get("symbol") == symbol), None)
    if hit is None:
        raise Unavailable(f"{symbol}: not found by /search")
    return {
        "symbol": symbol,
        "name": hit["name"],
        "aliases": hit.get("aliases") or [],
        "hgnc_id": hit["hgnc_id"],
        "ensembl_id": hit["ensembl_id"],
        "chr": hit["chrom"],
        "start": hit["gene_start"],
        "end": hit["gene_end"],
    }


def burden_section(base: str, symbol: str) -> dict:
    rows = tsv(fetch(base, f"/gene_based/{q(symbol)}"))
    best: dict[tuple, dict] = {}
    masks = Counter()
    for r in rows:
        m = num(r["mlog10p_burden"])
        if m is None:
            continue
        key = (r["dataset"], r["trait"])
        if m >= BURDEN_MLOG10P_MIN:
            masks[key] += 1
        if key not in best or m > num(best[key]["mlog10p_burden"]):
            best[key] = r

    def row(r):
        key = (r["dataset"], r["trait"])
        binary = r["n_controls"] is not None
        return {
            "dataset": r["dataset"],
            "trait": r["trait"],
            "mask": r["annotation"],
            **stats(num(r["mlog10p_burden"]), num(r["beta"]), binary),
            "n_cases": num(r["n_cases"]) and int(num(r["n_cases"])),
            "n_controls": num(r["n_controls"]) and int(num(r["n_controls"])),
            "significant_masks": masks[key],
        }

    ranked = sorted(best.values(), key=lambda r: -num(r["mlog10p_burden"]))
    hits_all = [r for r in ranked if num(r["mlog10p_burden"]) >= BURDEN_MLOG10P_MIN]
    per_ds = Counter()
    shown = []
    for r in hits_all:
        if per_ds[r["dataset"]] < CAPS["burden_per_dataset"] and len(shown) < CAPS["burden"]:
            per_ds[r["dataset"]] += 1
            shown.append(row(r))
    hit_datasets = {r["dataset"] for r in hits_all}
    # datasets whose rows are served unfiltered: the best row of each is a real null result
    nulls = {}
    for r in ranked:
        if r["dataset"] not in hit_datasets and r["dataset"] not in nulls:
            nulls[r["dataset"]] = row(r)
    return {
        "note": f"one row per dataset and trait, its strongest mask; hits have p <= 1e-{BURDEN_MLOG10P_MIN}",
        "n_hits": len(hits_all),
        "n_hits_by_dataset": dict(Counter(r["dataset"] for r in hits_all)),
        "hits": shown,
        "tested_without_hit": list(nulls.values()),
    }


def exome_section(api: str, symbol: str) -> dict:
    by_variant = defaultdict(list)
    for r in tsv(fetch(api, f"/exome_results_by_gene/{q(symbol)}")):
        # a p-value that underflowed to 0 arrives as an infinite -log10(p); 324 is the floor
        # the IBD exome source stores, so both read the same way
        m = min(float(r["mlog10p"]), MLOG10P_FLOOR)
        if m >= EXOME_MLOG10P_MIN:
            by_variant[f"{r['chr']}:{r['pos']}:{r['ref']}:{r['alt']}"].append({**r, "mlog10p": m})
    variants = []
    for variant, rs in by_variant.items():
        rs.sort(key=lambda r: -r["mlog10p"])
        top = rs[0]
        variants.append({
            "variant": variant,
            "annotation": top["annotation"],
            "allele_frequency": sig(num(top["af_overall"])),
            "n_significant_traits": len(rs),
            "datasets": sorted({r["dataset"] for r in rs}),
            "traits": [{
                "dataset": r["dataset"],
                "trait": r["trait"],
                **stats(r["mlog10p"], num(r["beta"]),
                        r["n_controls"] is not None or r["dataset"] in BINARY_EXOME_DATASETS),
            } for r in rs[:CAPS["exome_traits_per_variant"]]],
            "_rank": top["mlog10p"],
        })
    variants.sort(key=lambda v: -v.pop("_rank"))
    return {
        "note": (
            f"exome-wide significant (p < 1e-{EXOME_MLOG10P_MIN}) single-variant exome "
            "associations, grouped by variant. "
            "A p of 1.0e-324 is the smallest the source stores, read it as p < 1e-323"
        ),
        "n_significant_variants": len(variants),
        "n_by_dataset": dict(Counter(d for v in variants for d in v["datasets"])),
        "variants": variants[:CAPS["exome_variants"]],
    }


def credible_set_section(base: str, gene: dict) -> dict:
    symbol = gene["symbol"]
    d = json.loads(fetch(base, f"/gene_association/{q(symbol)}"))
    v, t, s, m = d["variants"], d["traits"], d["sets"], d["members"]

    def vid(i):
        return f"{d['chr']}:{v['pos'][i]}:{v['ref'][i]}:{v['alt'][i]}"

    def cons(i):
        c = v["consequence"][i]
        return d["consequences"][c] if c >= 0 else None

    def vgene(i):
        g = v["gene"][i]
        return d["genes"][g] if g >= 0 else None

    sets = []
    coding = defaultdict(list)
    offset = 0
    for i in range(len(s["trait"])):
        n = s["memberCount"][i]
        members = list(zip(m["variant"][offset:offset + n], m["pip"][offset:offset + n]))
        offset += n
        ti = s["trait"][i]
        binary = t["type"][ti] == 1
        lead = s["lead"][i]
        lead_pip = next((p for vi, p in members if vi == lead), None)
        entry = {
            "trait": t["name"][ti].replace("_", " "),
            "dataset": d["datasets"][s["dataset"][i]],
            "kind": "pseudo credible set" if s["pseudo"][i] else "fine-mapped",
            "set_size": s["size"][i],
            "lead_variant": vid(lead),
            **location(gene, vid(lead)),
            "lead_consequence": cons(lead),
            "lead_gene": vgene(lead),
            "lead_pip": None if lead_pip is None else lead_pip / 10000,
            **stats(s["leadMlog10p"][i], s["leadBeta"][i], binary),
            "significance": (
                "genome-wide significant" if (s["leadMlog10p"][i] or 0) >= GWS_MLOG10P
                else "below genome-wide significance"
            ),
        }
        sets.append(entry)
        for vi, p in members:
            c = cons(vi)
            if (c and c.replace("_variant", "") in CODING and vgene(vi) == symbol
                    and p / 10000 > CODING_PIP_MIN):
                coding[vi].append({
                    "trait": entry["trait"], "dataset": entry["dataset"], "kind": entry["kind"],
                    "pip": p / 10000, "set_size": entry["set_size"],
                    "set_lead_p": entry["p"],
                })

    coding_rows = []
    for vi, traits in coding.items():
        traits.sort(key=lambda x: -x["pip"])
        coding_rows.append({
            "variant": vid(vi),
            "consequence": cons(vi),
            "n_traits": len(traits),
            "max_pip": traits[0]["pip"],
            "traits": traits[:CAPS["traits_per_coding_variant"]],
        })
    coding_rows.sort(key=lambda r: (-r["max_pip"], -r["n_traits"]))

    best_per_trait = {}
    for e in sets:
        k = (e["trait"], e["dataset"])
        if k not in best_per_trait or (e["mlog10p"] or 0) > (best_per_trait[k]["mlog10p"] or 0):
            best_per_trait[k] = e
    ranked = sorted(best_per_trait.values(), key=lambda e: -(e["mlog10p"] or 0))
    return {
        "note": (
            "GWAS credible sets with at least one member variant inside the gene; the lead "
            "variant itself can lie outside it. A pseudo credible set is approximated from "
            "summary statistics and LD, not fine-mapped. Genome-wide significance is p < 5e-8"
        ),
        "n_sets": len(sets),
        "n_traits": len(best_per_trait),
        "n_sets_by_dataset": dict(Counter(e["dataset"] for e in sets)),
        "coding_variants_note": (
            f"coding variants of this gene with PIP > {CODING_PIP_MIN} in a credible set"
        ),
        "n_coding_variants": len(coding_rows),
        "coding_variants": coding_rows[:CAPS["coding_variants"]],
        "strongest_sets": ranked[:CAPS["gwas_sets"]],
    }


def qtl_section(base: str, gene: dict) -> dict:
    symbol = gene["symbol"]
    rows = tsv(fetch(base, f"/credible_sets_by_qtl_gene/{q(symbol)}"))
    sets = {}
    for r in rows:
        k = (r["resource"], r["dataset"], r["data_type"], r["trait_original"], r["cell_type"], r["cs_id"])
        pip = num(r["pip"]) or 0
        if k not in sets or pip > sets[k]["_pip"]:
            sets[k] = {**r, "_pip": pip}
    by_type = defaultdict(list)
    for r in sets.values():
        by_type[r["data_type"]].append(r)
    out = {}
    for dtype, rs in sorted(by_type.items()):
        rs.sort(key=lambda r: -(num(r["mlog10p"]) or 0))
        out[dtype] = {
            "n_credible_sets": len(rs),
            "n_tissues_or_cell_types": len({r["cell_type"] for r in rs}),
            "n_datasets": len({r["dataset"] for r in rs}),
            "strongest": [{
                "dataset": r["dataset"],
                "tissue_or_cell_type": r["cell_type"],
                "molecular_trait": r["trait_original"] if dtype == "caQTL" else r["trait"],
                "lead_variant": f"{r['chr']}:{r['pos']}:{r['ref']}:{r['alt']}",
                **location(gene, f"{r['chr']}:{r['pos']}:{r['ref']}:{r['alt']}"),
                "lead_pip": sig(r["_pip"]),
                "set_size": num(r["cs_size"]) and int(num(r["cs_size"])),
                **stats(num(r["mlog10p"]), num(r["beta"]), False),
            } for r in rs[:CAPS["qtl_sets_per_type"]]],
        }
    return {
        "note": "molecular QTL credible sets whose target is this gene; lead = highest-PIP member",
        "n_credible_sets": len(sets),
        "by_type": out,
    }


def coloc_sections(api: str, names: dict, gene: dict) -> tuple[dict, dict]:
    """Cis and trans QTL-GWAS colocalizations of the gene, as two sections."""
    cis, trans = defaultdict(list), defaultdict(list)
    bundle_gene = gene
    for r in tsv(fetch(api, f"/colocalization_by_gene/{q(gene['symbol'])}")):
        s = r["qtl_side"]
        g = "1" if s == "2" else "2"
        if r[f"data_type{g}"] != "GWAS" or float(r["PP.H4.abf"]) <= COLOC_PP_H4_MIN:
            continue
        qtl_pos = int(r[f"hit{s}"].split(":")[1])
        is_cis = (
            r["chr"] == r["gene_chr"]
            and int(r["gene_start"]) - CIS_WINDOW <= qtl_pos <= int(r["gene_end"]) + CIS_WINDOW
        )
        row = {
            "trait": r[f"trait{g}"], "dataset": r[f"dataset{g}"],
            "qtl_type": r[f"data_type{s}"], "qtl_dataset": r[f"dataset{s}"],
            "cell_type": r[f"cell_type{s}"], "pp_h4": float(r["PP.H4.abf"]),
            "gwas_hit": r[f"hit{g}"], "qtl_hit": r[f"hit{s}"],
            "gwas_beta": num(r[f"hit{g}_beta"]), "qtl_beta": num(r[f"hit{s}_beta"]),
            "gwas_mlog10p": num(r[f"hit{g}_mlog10p"]) or 0,
        }
        (cis if is_cis else trans)[(row["dataset"], row["trait"])].append(row)

    def per_trait(groups: dict, cap: int) -> tuple[list, int]:
        out = []
        for (dataset, trait), rs in groups.items():
            rs.sort(key=lambda r: -r["pp_h4"])
            b = rs[0]
            same = b["gwas_hit"] == b["qtl_hit"] and b["gwas_beta"] and b["qtl_beta"]
            out.append({
                "gwas_trait": names.get(trait, trait).replace("_", " "),
                "gwas_dataset": dataset,
                "n_qtl_contexts": len(rs),
                "qtl_types": sorted({r["qtl_type"] for r in rs}),
                "best": {
                    "qtl_type": b["qtl_type"],
                    "qtl_dataset": b["qtl_dataset"],
                    "tissue_or_cell_type": b["cell_type"],
                    "pp_h4": sig(b["pp_h4"]),
                    "gwas_lead_variant": b["gwas_hit"],
                    "qtl_lead_variant": b["qtl_hit"],
                    "qtl_lead": location(bundle_gene, b["qtl_hit"]),
                    "gwas_lead": location(bundle_gene, b["gwas_hit"]),
                    "gwas_p": pval(b["gwas_mlog10p"]),
                    "direction": "lead variants differ" if not same else (
                        "higher QTL level goes with higher trait/risk"
                        if b["gwas_beta"] * b["qtl_beta"] > 0
                        else "higher QTL level goes with lower trait/risk"
                    ),
                },
                "_rank": (b["gwas_mlog10p"], len(rs)),
            })
        n_gws = sum(1 for o in out if o["_rank"][0] >= GWS_MLOG10P)
        out.sort(key=lambda o: o.pop("_rank"), reverse=True)
        return out[:cap], n_gws

    cis_rows, cis_gws = per_trait(cis, CAPS["coloc"])
    trans_rows, trans_gws = per_trait(trans, CAPS["coloc_trans"])
    loci = Counter()
    for (_d, _t), rs in trans.items():
        for hit in {r["qtl_hit"] for r in rs}:
            loci[hit] += 1
    shared = (
        f"GWAS signals colocalizing (PP.H4 > {COLOC_PP_H4_MIN}) with an eQTL, pQTL or sQTL of "
        "this gene; one row per GWAS trait, its best QTL context. "
    )
    return (
        {
            "note": shared + f"CIS: the QTL signal lies within {CIS_WINDOW // 1_000_000} Mb of the gene",
            "n_gwas_traits": len(cis),
            "n_gwas_traits_genome_wide_significant": cis_gws,
            "traits": cis_rows,
        },
        {
            "note": shared + (
                "TRANS: the QTL signal lies elsewhere in the genome, so the shared variant acts "
                "on this gene's product from a distance. `loci` counts the traits behind each "
                "trans QTL lead variant: one locus behind many traits is usually a pleiotropic "
                "hub rather than evidence about this gene"
            ),
            "n_gwas_traits": len(trans),
            "n_gwas_traits_genome_wide_significant": trans_gws,
            "loci": [{"qtl_lead_variant": v, **location(gene, v), "n_gwas_traits": n}
                     for v, n in loci.most_common(CAPS["coloc_trans_loci"])],
            "traits": trans_rows,
        },
    )


def gtex_section(base: str, ensembl_id: str) -> dict:
    d = json.loads(fetch(base, f"/gtex_expression?ensembl_id={q(ensembl_id)}"))
    med = sorted(
        ((t["tissue"].replace("_", " "), statistics.median(t["values"])) for t in d["tissues"] if t["values"]),
        key=lambda x: -x[1],
    )
    return {
        "note": f"bulk RNA-seq, median {d['unit']} per tissue, GTEx {d['release']}",
        "n_tissues": len(med),
        "n_tissues_median_above_1": sum(1 for _, m in med if m > 1),
        "median_across_tissues": sig(statistics.median(m for _, m in med)) if med else None,
        "top_tissues": [{"tissue": t, "median_tpm": sig(m)} for t, m in med[:CAPS["gtex_tissues"]]],
    }


def hpa_section(base: str, symbol: str) -> dict:
    rows = [r for r in tsv(fetch(base, f"/expression_by_gene/{q(symbol)}")) if r["resource"] == "hpa"]
    order = {"High": 0, "Medium": 1, "Low": 2}
    shown = sorted((r for r in rows if r["level"] in ("High", "Medium")), key=lambda r: order[r["level"]])
    return {
        "note": "Human Protein Atlas protein level by tissue and cell type (immunohistochemistry)",
        "n_by_level": dict(Counter(r["level"] for r in rows)),
        "high_or_medium": [{
            "tissue": r["tissue_cell"].split("|")[0].replace("_", " "),
            "cell_type": r["tissue_cell"].split("|")[-1].replace("_", " "),
            "level": r["level"],
        } for r in shown[:CAPS["hpa"]]],
    }


def cell_type_section(base: str, symbol: str) -> dict:
    rows = tsv(fetch(base, f"/cell_type_expression/{q(symbol)}"))
    ok = [r for r in rows
          if (num(r["n_cells"]) or 0) >= CELL_MIN_CELLS and (num(r["n_donors"]) or 0) >= CELL_MIN_DONORS]
    ok.sort(key=lambda r: -(num(r["mean_log1p_cp10k"]) or 0))
    return {
        "note": (
            "CELLxGENE single-cell expression by cell type and tissue; groups with at least "
            f"{CELL_MIN_CELLS} cells and {CELL_MIN_DONORS} donors. tau is a 0-1 specificity "
            "index, 1 = expressed in one group only"
        ),
        "n_groups": len(ok),
        "tau": sig(num(rows[0]["tau"])) if rows else None,
        "top_groups": [{
            "cell_type": r["cell_type_rollup"],
            "tissue": r["tissue_general"],
            "percent_cells_expressing": sig(100 * (num(r["frac_expressing"]) or 0)),
            "mean_log1p_cp10k": sig(num(r["mean_log1p_cp10k"])),
            "n_donors": int(num(r["n_donors"])),
        } for r in ok[:CAPS["cell_types"]]],
    }


def mendelian_section(base: str, symbol: str) -> dict:
    rows = tsv(fetch(base, f"/gene_disease/{q(symbol)}"))
    by_disease = defaultdict(list)
    for r in rows:
        by_disease[r["disease_title"]].append(r)

    def rank(c):
        return GENCC_RANK.index(c) if c in GENCC_RANK else len(GENCC_RANK)

    out = []
    for disease, rs in by_disease.items():
        gencc = [r for r in rs if r["resource"] == "gencc"]
        best = min((r["classification"] for r in gencc), key=rank, default=None)
        out.append({
            "disease": disease,
            "best_gencc_classification": best,
            "inheritance": sorted({r["mode_of_inheritance"] for r in rs
                                   if r["mode_of_inheritance"] not in (None, "Unknown")}),
            "n_gencc_submitters": len({r["submitter"] for r in gencc}),
            "in_monarch": any(r["resource"] == "monarch" for r in rs),
        })
    out.sort(key=lambda o: (rank(o["best_gencc_classification"]), -o["n_gencc_submitters"]))
    return {
        "note": "curated Mendelian gene-disease relationships (GenCC submissions, Monarch)",
        "n_diseases": len(out),
        "diseases": out[:CAPS["diseases"]],
    }


def constraint_section(base: str, ensembl_id: str) -> dict:
    d = json.loads(fetch(base, f"/gene_constraint?ensembl_id={q(ensembl_id)}", tries=4, timeout=60))
    if "error" in d:
        raise Unavailable(d.get("message", "constraint error"))
    return {
        "note": "gnomAD constraint (LOEUF below 0.6 is commonly read as LoF-intolerant) and ClinVar counts",
        "pLI": sig(d.get("pLI")),
        "loeuf": sig(d.get("loeuf")),
        "flags": d.get("flags") or [],
        "clinvar_variants": d.get("clinvar"),
    }


# ---------------------------------------------------------------- assembly

ID_PREFIX = {
    ("burden", "hits"): "B", ("burden", "tested_without_hit"): "N",
    ("exome_variants", "variants"): "X",
    ("gwas_credible_sets", "coding_variants"): "C", ("gwas_credible_sets", "strongest_sets"): "G",
    ("colocalization_cis", "traits"): "L", ("colocalization_trans", "traits"): "R",
    ("expression_gtex", "top_tissues"): "T", ("expression_hpa", "high_or_medium"): "H",
    ("expression_cell_types", "top_groups"): "S",
    ("mendelian", "diseases"): "M",
}


def add_ids(bundle: dict) -> None:
    for (section, key), prefix in ID_PREFIX.items():
        rows = bundle[section].get(key) if bundle[section].get("status") == "ok" else None
        for i, row in enumerate(rows or [], 1):
            bundle[section][key][i - 1] = {"id": f"{prefix}{i}", **row}
    qtl = bundle["qtl"]
    if qtl.get("status") == "ok":
        n = 0
        for dtype in qtl["by_type"].values():
            for i, row in enumerate(dtype["strongest"]):
                n += 1
                dtype["strongest"][i] = {"id": f"Q{n}", **row}
    for section, short in (("constraint", "K1"), ("expression_gtex", "T0"), ("expression_cell_types", "S0")):
        if bundle[section].get("status") == "ok":
            bundle[section] = {"id": short, **bundle[section]}


def build(base: str, api: str, symbol: str, names: dict) -> dict:
    gene = gene_section(base, symbol)
    try:
        coloc = coloc_sections(api, names, gene)
    except (Unavailable, KeyError, ValueError) as e:
        coloc = e

    def coloc_part(i: int):
        if isinstance(coloc, Exception):
            raise Unavailable(str(coloc))
        return coloc[i]

    sections = {
        "burden": lambda: burden_section(base, symbol),
        "exome_variants": lambda: exome_section(api, symbol),
        "gwas_credible_sets": lambda: credible_set_section(base, gene),
        "qtl": lambda: qtl_section(base, gene),
        "colocalization_cis": lambda: coloc_part(0),
        "colocalization_trans": lambda: coloc_part(1),
        "expression_gtex": lambda: gtex_section(base, gene["ensembl_id"]),
        "expression_hpa": lambda: hpa_section(base, symbol),
        "expression_cell_types": lambda: cell_type_section(base, symbol),
        "mendelian": lambda: mendelian_section(base, symbol),
        "constraint": lambda: constraint_section(base, gene["ensembl_id"]),
    }
    bundle = {"bundle_version": BUNDLE_VERSION, "gene": gene}
    for name, fn in sections.items():
        try:
            bundle[name] = {"status": "ok", **fn()}
        except (Unavailable, KeyError, ValueError) as e:
            bundle[name] = {"status": "unavailable", "reason": str(e)[:200]}
    add_ids(bundle)
    canonical = json.dumps(bundle, sort_keys=True, separators=(",", ":"))
    bundle["bundle_sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
    return bundle


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--genes", required=True, help="file with one gene symbol per line (first column)")
    ap.add_argument("--out", required=True, help="directory for <SYMBOL>.json bundles")
    ap.add_argument("--api", default="http://127.0.0.1:5001/api/v1", help="the browser BFF's /api/v1")
    ap.add_argument(
        "--results-api",
        help="results-api /api/v1 for the exome and colocalization by-gene routes (default: --api)",
    )
    ap.add_argument("--workers", type=int, default=3)
    args = ap.parse_args()

    genes = [ln.split()[0] for ln in Path(args.genes).read_text().splitlines() if ln.strip()]
    api = args.results_api or args.api
    names = {r["trait"]: r["name"] for r in tsv(fetch(args.api, "/trait_name_mapping")) if r["name"]}
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    def one(symbol: str):
        try:
            b = build(args.api, api, symbol, names)
        except Exception as e:  # one gene's failure is reported, not allowed to stop the run
            return symbol, None, f"{type(e).__name__}: {e}"
        (out / f"{symbol}.json").write_text(json.dumps(b, indent=1))
        bad = [k for k, v in b.items() if isinstance(v, dict) and v.get("status") == "unavailable"]
        return symbol, len(json.dumps(b, separators=(",", ":"))), ",".join(bad)

    failed = 0
    with ThreadPoolExecutor(args.workers) as pool:
        for symbol, size, note in pool.map(one, genes):
            failed += size is None
            print(f"{symbol}\t{size}\t{note}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
