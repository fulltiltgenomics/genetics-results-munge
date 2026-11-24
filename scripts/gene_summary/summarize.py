#!/usr/bin/env python3
"""Write an AI gene summary from each evidence bundle, then check it against the bundle.

One model call per gene, no tools: the bundle is the only source. The summary is JSON whose
findings each cite bundle row ids, and the validator rejects a summary that cites an id the
bundle does not have, or states a number or a variant that is not in it.
"""

import argparse
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import anthropic

from check import TABS, validate

PROMPT_VERSION = 3
MODEL = "claude-opus-5-5"
EFFORT = "medium"

SYSTEM = """\
You write the "AI summary" shown at the top of a gene's page in a human genetics results \
browser. Readers are geneticists and biologists deciding whether the gene is worth their time \
and which evidence tab to open next. They can see the full tables in the other tabs; your job \
is to tell them what the evidence adds up to.

You are given one JSON evidence bundle for one gene. It was assembled by fixed rules from the \
browser's own data: gene burden tests, single-variant exome associations, GWAS credible sets, \
molecular QTLs, QTL-GWAS colocalization, expression (GTEx, Human Protein Atlas, CELLxGENE), \
curated Mendelian gene-disease relationships and gnomAD constraint.

The bundle is your only source. Do not add anything you know about the gene from the \
literature or elsewhere: no function, pathway, mechanism, drug or disease role unless the \
bundle shows it. The page will carry a separate literature section later, and readers must be \
able to trust that every statement here traces to a row they can open. You may interpret what \
the rows show together (for example that rare coding burden, a fine-mapped missense variant \
and common-variant signals point at the same disease), and you should say when evidence types \
disagree or when one is absent.

How to read the bundle:
- Rows carry an `id`. Each section's `n_*` fields are the full counts; the listed rows are the \
strongest ones, capped. When you say how many there are, use the counts, not the list length.
- A section with `"status": "unavailable"` could not be fetched. Say nothing about it except, \
if it matters, that it was unavailable. Never present it as an absence of evidence.
- A section that is ok with zero rows is a real absence, and worth stating when a reader would \
expect a signal.
- `burden.tested_without_hit` lists datasets where the gene was tested and nothing passed the \
cut: a real null result.
- `kind: "pseudo credible set"` means the set was approximated from summary statistics and LD \
rather than fine-mapped. A high PIP there is weaker evidence than in a fine-mapped set; keep \
the distinction when it affects the conclusion.
- `odds_ratio` is given for binary traits. For quantitative traits `beta` is in the dataset's \
own units and is not comparable across datasets; give its direction rather than its size.
- Significance terms are tied to the test. "Exome-wide significant" is for single-variant \
exome associations (p < 1e-6); every variant listed in that section meets it. "Genome-wide \
significant" is for GWAS credible sets (p < 5e-8), and those rows carry a `significance` \
label. Gene burden results have no such label: give their p-value and do not call them \
exome-wide or genome-wide significant.
- `distance_to_gene_kb` is how far a lead variant lies from the gene (0 inside it). When \
distance matters to the argument, quote that number rather than estimating one from positions.
- Colocalization `direction` is only stated when both studies share a lead variant.
- Cis colocalization (the QTL signal is at the gene) is evidence that a disease or trait \
signal acts through this gene. Trans colocalization (the QTL signal is elsewhere in the \
genome) says a distant locus changes this gene's product and also associates with the trait: \
weaker, and when one trans locus stands behind many traits it is a pleiotropic hub, to be \
mentioned in a sentence rather than treated as evidence about this gene.
- Trait names from large biobanks are often near-duplicates of one phenotype. Group them \
rather than listing them.

Numbers: copy every number exactly as it appears in the bundle (p-values as the `p` string, \
odds ratios, PIPs, counts, TPM). Do not convert, re-round, combine or derive numbers, and do not \
state a cut-off or threshold the bundle does not give: say "suggestive" or "did not pass the \
cut" instead of inventing one. Write variants exactly as given (chr:pos:ref:alt).

Weigh the evidence rather than recite it. A genome-wide significant coding variant with high \
PIP in a fine-mapped set, or an exome-wide significant burden result, is strong evidence that \
the gene itself is involved. A non-coding credible set that merely overlaps the gene is not, \
on its own, evidence for this gene over its neighbours. Many weak phenome-wide hits are \
usually noise; say so briefly instead of listing them.

Output fields:
- `headline`: one sentence, at most 35 words, the single most important thing the data says \
about this gene. If the data says little, say that.
- `overview`: one paragraph, 60 to 130 words, plain prose, no lists, no ids.
- `findings`: 3 to 7 items, most important first. Each is one or two sentences making one \
point, with `evidence` listing the ids of the rows it rests on (at least one) and `tab` naming \
the page tab where the reader can see it.
- `gaps`: 0 to 3 short sentences on what is absent, null or unavailable and matters for \
interpreting the rest.

Write for a reader who has not seen the bundle and does not know it exists: full sentences, \
no field names, no arrows or shorthand, and never the word "bundle" (say "the data" or name \
the evidence)."""

SCHEMA = {
    "type": "object",
    "properties": {
        "headline": {"type": "string"},
        "overview": {"type": "string"},
        "findings": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "statement": {"type": "string"},
                    "evidence": {"type": "array", "items": {"type": "string"}},
                    "tab": {"type": "string", "enum": TABS},
                },
                "required": ["statement", "evidence", "tab"],
                "additionalProperties": False,
            },
        },
        "gaps": {"type": "array", "items": {"type": "string"}},
    },
    "required": ["headline", "overview", "findings", "gaps"],
    "additionalProperties": False,
}

# ---------------------------------------------------------------- generation


def summarize(client: anthropic.Anthropic, bundle: dict, feedback: list[str] | None = None) -> dict:
    user = "Evidence bundle:\n\n" + json.dumps(bundle, separators=(",", ":"))
    if feedback:
        user += (
            "\n\nA previous summary of this bundle was rejected by an automatic check for "
            "these reasons. Write a new one that avoids them:\n- " + "\n- ".join(feedback)
        )
    t0 = time.time()
    response = client.messages.create(
        model=MODEL,
        max_tokens=16000,
        thinking={"type": "adaptive"},
        output_config={"effort": EFFORT, "format": {"type": "json_schema", "schema": SCHEMA}},
        system=[{"type": "text", "text": SYSTEM, "cache_control": {"type": "ephemeral"}}],
        messages=[{"role": "user", "content": user}],
    )
    u = response.usage
    meta = {
        "model": response.model,
        "stop_reason": response.stop_reason,
        "seconds": round(time.time() - t0, 1),
        "input_tokens": u.input_tokens,
        "cache_read_input_tokens": u.cache_read_input_tokens or 0,
        "cache_creation_input_tokens": u.cache_creation_input_tokens or 0,
        "output_tokens": u.output_tokens,
    }
    if response.stop_reason != "end_turn":
        return {"summary": None, "meta": meta}
    text = next(b.text for b in response.content if b.type == "text")
    return {"summary": json.loads(text), "meta": meta}


def run_gene(client, bundle_path: Path, out_dir: Path, retries: int) -> dict:
    full = json.loads(bundle_path.read_text())
    sha = full.pop("bundle_sha256")
    symbol = full["gene"]["symbol"]
    attempts, problems, result = [], None, None
    for _ in range(1 + retries):
        try:
            result = summarize(client, full, problems)
        except anthropic.APIStatusError as e:
            attempts.append({"error": f"{e.status_code} {e.message}"[:300]})
            break
        except anthropic.APIConnectionError as e:
            attempts.append({"error": f"connection: {e}"[:300]})
            break
        problems = (
            validate(result["summary"], full) if result["summary"]
            else [f"stop_reason {result['meta']['stop_reason']}"]
        )
        attempts.append({**result["meta"], "problems": problems, "summary": result["summary"]})
        if not problems or result["summary"] is None:
            break
    record = {
        "gene": symbol,
        "bundle_sha256": sha,
        "prompt_version": PROMPT_VERSION,
        "model": MODEL,
        "effort": EFFORT,
        "valid": bool(result and result["summary"] and not problems),
        "summary": result and result["summary"],
        "attempts": attempts,
    }
    (out_dir / f"{symbol}.json").write_text(json.dumps(record, indent=1))
    return record


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bundles", required=True, help="directory of <SYMBOL>.json bundles")
    ap.add_argument("--out", required=True, help="directory for <SYMBOL>.json summaries")
    ap.add_argument("--genes", nargs="*", help="only these symbols")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--retries", type=int, default=1, help="regenerations after a failed check")
    ap.add_argument("--force", action="store_true", help="rewrite summaries that already exist")
    args = ap.parse_args()
    if not os.environ.get("ANTHROPIC_API_KEY"):
        print("ANTHROPIC_API_KEY is not set", file=sys.stderr)
        return 2

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    paths = sorted(Path(args.bundles).glob("*.json"))
    if args.genes:
        paths = [p for p in paths if p.stem in set(args.genes)]
    if not args.force:
        paths = [p for p in paths if not (out / p.name).exists()]
    client = anthropic.Anthropic()

    invalid = 0
    with ThreadPoolExecutor(args.workers) as pool:
        for rec in pool.map(lambda p: run_gene(client, p, out, args.retries), paths):
            last = rec["attempts"][-1]
            invalid += not rec["valid"]
            print(
                f"{rec['gene']}\tvalid={rec['valid']}\tattempts={len(rec['attempts'])}\t"
                f"in={last.get('input_tokens')}\tout={last.get('output_tokens')}\t"
                f"{'; '.join(last.get('problems') or [last.get('error', '')])[:200]}",
                flush=True,
            )
    return 1 if invalid else 0


if __name__ == "__main__":
    sys.exit(main())
