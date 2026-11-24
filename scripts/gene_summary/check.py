#!/usr/bin/env python3
"""Check gene summaries against their evidence bundles, and file the ones written by hand.

A summary is rejected when it cites a row id the bundle does not have, states a number or a
variant that is not in the bundle, or breaks the length rules. `summarize.py` runs the same
check on what the API returns; this entry point is for summaries written outside it, one
`{headline, overview, findings, gaps}` JSON per gene, and with `--out` it files each as the
record `summarize.py` would have written.
"""

import argparse
import json
import re
import sys
from pathlib import Path

TABS = ["gene burden", "variant association", "qtl", "expression", "mendelian", "annotation"]

# a hyphen after a letter is part of a name (COVID-19, IL-12), not a minus sign
NUMBER = re.compile(
    r"(?<![\w.:])(?:(?<![\w-])-)?\d+(?:,\d{3})*(?:\.\d+)?(?:e-?\d+)?(?![\w:]|\.\d)", re.I
)
# a chromosome read off a variant id: "chromosome 12", "chromosomes 1, 9 and 12"
CHROMOSOMES = re.compile(r"\bchromosomes?\s+\d{1,2}(?:\s*(?:,|and|or)\s*\d{1,2})*", re.I)
VARIANT = re.compile(r"\b(?:\d{1,2}|X|Y):\d+:[ACGT]+:[ACGT]+\b")
# cut-offs the prompt and the bundle notes state in words, and small counts a writer derives
# by reading ("two masks", "three diseases")
FREE_NUMBERS = {"5e-8", "1e-6", "1e-4", "0.05", "0.6", "0.8", "1e-323"}


def bundle_numbers(node, out: set, ids: dict, variants: set) -> None:
    if isinstance(node, dict):
        if isinstance(node.get("id"), str):
            ids[node["id"]] = node
        for k, v in node.items():
            if k != "id":
                bundle_numbers(v, out, ids, variants)
    elif isinstance(node, list):
        for v in node:
            bundle_numbers(v, out, ids, variants)
    elif isinstance(node, bool) or node is None:
        return
    elif isinstance(node, (int, float)):
        out.add(canon(repr(node)))
    elif isinstance(node, str):
        variants.update(VARIANT.findall(node))
        for m in NUMBER.findall(VARIANT.sub(" ", node)):
            out.add(canon(m))


def canon(token: str) -> str:
    t = token.replace(",", "").lower()
    try:
        f = float(t)
    except ValueError:
        return t
    if "e" in t:
        mant, exp = t.split("e")
        return f"{float(mant):g}e{int(exp)}"
    return f"{f:g}"


def validate(summary: dict, bundle: dict) -> list[str]:
    numbers, ids, variants = set(), {}, set()
    bundle_numbers(bundle, numbers, ids, variants)
    numbers |= {canon(n) for n in FREE_NUMBERS}
    problems = []

    def check_text(where: str, text: str) -> None:
        for v in VARIANT.findall(text):
            if v not in variants:
                problems.append(f"{where}: variant {v} is not in the bundle")
        for tok in NUMBER.findall(CHROMOSOMES.sub(" ", VARIANT.sub(" ", text))):
            c = canon(tok)
            try:
                small_count = float(c).is_integer() and 0 <= float(c) <= 10
            except ValueError:
                small_count = False
            if c not in numbers and not small_count:
                problems.append(f"{where}: number {tok} is not in the bundle")

    check_text("headline", summary["headline"])
    check_text("overview", summary["overview"])
    for i, gap in enumerate(summary["gaps"]):
        check_text(f"gaps[{i}]", gap)
    if not 3 <= len(summary["findings"]) <= 7:
        problems.append(f"findings: {len(summary['findings'])} items, expected 3 to 7")
    for i, f in enumerate(summary["findings"]):
        if not f["evidence"]:
            problems.append(f"findings[{i}]: cites no evidence")
        for e in f["evidence"]:
            if e not in ids:
                problems.append(f"findings[{i}]: cites unknown id {e}")
        check_text(f"findings[{i}]", f["statement"])
    words = len(summary["overview"].split())
    if not 40 <= words <= 160:
        problems.append(f"overview: {words} words")
    return problems

def shape_problems(summary) -> list[str]:
    """What the API's JSON schema would have refused; empty when the shape is right."""
    if not isinstance(summary, dict) or set(summary) != {"headline", "overview", "findings", "gaps"}:
        return ["top level must have exactly headline, overview, findings, gaps"]
    problems = []
    if not isinstance(summary["headline"], str) or not isinstance(summary["overview"], str):
        problems.append("headline and overview must be strings")
    if not isinstance(summary["gaps"], list) or not all(isinstance(g, str) for g in summary["gaps"]):
        problems.append("gaps must be a list of strings")
    if not isinstance(summary["findings"], list):
        return problems + ["findings must be a list"]
    for i, f in enumerate(summary["findings"]):
        if not isinstance(f, dict) or set(f) != {"statement", "evidence", "tab"}:
            problems.append(f"findings[{i}]: must have exactly statement, evidence, tab")
        elif (
            not isinstance(f["statement"], str)
            or not isinstance(f["evidence"], list)
            or not all(isinstance(e, str) for e in f["evidence"])
        ):
            problems.append(f"findings[{i}]: statement is a string, evidence a list of strings")
        elif f["tab"] not in TABS:
            problems.append(f"findings[{i}]: tab {f['tab']!r} is not one of {TABS}")
    return problems



def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bundles", required=True, help="directory of <SYMBOL>.json bundles")
    ap.add_argument("--summaries", required=True, help="directory of <SYMBOL>.json bare summaries")
    ap.add_argument("--genes", nargs="*", help="only these symbols")
    ap.add_argument("--out", help="directory for <SYMBOL>.json records")
    ap.add_argument("--model", help="what wrote the summaries, for the record (with --out)")
    ap.add_argument("--prompt-version", type=int, help="for the record (with --out)")
    args = ap.parse_args()

    paths = sorted(Path(args.summaries).glob("*.json"))
    if args.genes:
        missing = set(args.genes) - {p.stem for p in paths}
        for symbol in sorted(missing):
            print(f"{symbol}\tMISSING")
        paths = [p for p in paths if p.stem in set(args.genes)]
    else:
        missing = set()
    if args.out:
        Path(args.out).mkdir(parents=True, exist_ok=True)

    invalid = len(missing)
    for path in paths:
        bundle = json.loads((Path(args.bundles) / path.name).read_text())
        sha = bundle.pop("bundle_sha256")
        try:
            summary = json.loads(path.read_text())
            problems = shape_problems(summary) or validate(summary, bundle)
        except json.JSONDecodeError as e:
            summary, problems = None, [f"not JSON: {e}"]
        invalid += bool(problems)
        print(f"{path.stem}\t{'INVALID' if problems else 'ok'}")
        for p in problems:
            print(f"  - {p}")
        if args.out:
            record = {
                "gene": bundle["gene"]["symbol"],
                "bundle_sha256": sha,
                "prompt_version": args.prompt_version,
                "model": args.model,
                "valid": not problems,
                "summary": summary,
                "problems": problems,
            }
            (Path(args.out) / path.name).write_text(json.dumps(record, indent=1))
    return 1 if invalid else 0


if __name__ == "__main__":
    sys.exit(main())
