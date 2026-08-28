# genetics-results-munge

Turns external genetics results (GWAS sumstats, exome/burden results, QTLs, open
chromatin, colocalisation) into bgzipped, tabix-indexed TSVs that the results API
serves and BigQuery loads.


# Read the code first, this file second

The nearest existing script is the specification. Before writing or changing a munge,
read one that handles the same *shape* of input and follow it — its flags, its function
split, its output columns:

| input shape | read |
|---|---|
| GWAS sumstats aligned against gnomAD | `scripts/munge_aih.py` |
| exome variant / gene burden results | `scripts/munge_ibd_exome.py` |
| peak / open-chromatin BED-like tables | `scripts/munge_catlas.py` |
| FinnGen-format fine-mapping at scale | `wdl/munge_finngen_finemapping_results.wdl` |

Where the facts live: `README.md` says what exists and how to run it; `docs/` holds the
long-form derivations; each script's **header docstring** holds its source files, format
assumptions and flags; each script's final `select()` **is** its column definition.
Don't restate any of that here — this file carries only what reading the code can't
tell you, which is mostly the things that fail *silently*.


# Before munging a new input

Read the actual bytes of the file. Source papers and READMEs are routinely wrong about
their own delivery. Establish these, then write what you found into the script's header
docstring (see `munge_aih.py` for the shape):

1. genome build — output is always GRCh38; lift over if it isn't (`munge_calderon.py`)
2. which allele the effect size refers to, and whether it is the reference allele
3. whether p underflows to 0, and whether beta/se are there to recover it from
4. chromosome coding, duplicate rows, indel representation, missing-value token
5. for anything carrying gene names: the gencode version, determined by match rate
   against `/mnt/disks/data/gencode.v*.annotation.genes.tsv` — never assumed

If any of these stays ambiguous, ask. A wrong guess here produces a file that loads,
indexes and queries fine while being wrong.


# Invariants

Break one of these and nothing errors — the data is just quietly wrong or unfindable.

**Statistics**
- `mlog10p` only, never a raw p-value; rounded to 4 decimals
- when p underflows to 0, recompute from beta/se with `scipy.special.log_ndtr`
- `beta` and `se` formatted `:.3e`

**Coordinates**
- GRCh38, `#chr` as an integer, **X is 23** (exome munges also map Y→24, MT→26)
- rows whose chr is null after the cast are malformed — drop them

**Alleles**
- `ref`/`alt` follow gnomAD, and `beta`/`af` are flipped when the source's effect allele
  turned out to be the reference. Each script's docstring states how it resolved this
  and what it dropped; indel representation is where this goes wrong, so state the
  decision explicitly rather than reusing another script's.

**Names** — `INFO` uppercase; `neff` lowercase; `n_cases`/`n_controls`,
`af_cases`/`af_controls`, `rsid`. The API maps these by exact string.

**Gene burden tabix is a point index** (`-s5 -b6 -e6`, both on `gene_start_pos`), so the
API finds a gene only if the file's coordinates come from exactly the `gencode_version`
configured for that dataset in the API's `gene_based_results.py`. Any other version
returns nothing, with no error. A new version also needs adding to `genes.py`'s
`gencode_versions` plus a regenerated mapping from
`create_gene_name_mapping_across_gencode_versions.py`.


# Output contract

- bgzipped TSV + tabix index, plus a `mlog10p > 4` filtered companion with its own index
- write it with `write_sumstat_output()` / `write_exome_output()` from
  `scripts/sumstat_utils.py`; only write your own bgzip pipe when the data can't fit a
  DataFrame (`munge_ibd.py`, `scripts/split_burden_per_trait.py`)
- gene burden datasets also ship unfiltered per-trait files, so that a gene's *null*
  result in a trait is still retrievable
- `--output` accepts `gs://`; staging to GCS is opt-in (`--stage`), never the default
- verify before staging: row counts at each filter step, a tabix query that returns a
  known variant, and the AF-vs-gnomAD plot (`--gnomad-af-plot`) where alleles were
  aligned. The peak munges take `--sample` for an end-to-end run on synthetic input.


# What changing an output touches

- **API** — column mappings live in `genetics-results-api/app/config/profiles/{finngen,daly}/summary_stats.py`.
  Both need updating for a new dataset or a renamed column.
- **BigQuery** — loaders and schemas live in `genetics-results-db`; check it in the same
  session rather than assuming it follows. `genetics-results-suite` is the spec of
  record for the suite as a whole.
- **Docs in this repo** — a doc is stale the moment it *enumerates* something the code
  no longer matches, so re-derive lists from the code instead of trusting them:

| changed | update |
|---|---|
| `scripts/munge_*`, `scripts/create_*` | `README.md` — dataset list, run instructions, flags |
| `scripts/sumstat_utils.py` | `CLAUDE.md` — the output contract above |
| `scripts/coloc/*` | `scripts/coloc/R14_UPDATE.md` |
| `wdl/munge_finngen_finemapping_results*`, `wdl/qtl_file.wdl` | `README.md` |
| `wdl/create_pseudo_credible_sets*`, `wdl/autoreporting_*.json` | `docs/pseudo-credible-sets.md` |

`scripts/check-doc-drift.sh` warns (never blocks) on commits that skip this; it runs
from `.beads/hooks/pre-commit`, outside the beads-managed fence that beads rewrites.


# Working agreements

- Don't guess. Ask rather than infer an assumption about the data.
- Don't change code unrelated to the request; flag it instead.
- Re-read a file before editing it — it may have changed since you last read it.
- Interactive mode: stage changes and propose a commit message, don't commit.
  Autonomous mode (ralph wiggum): commit each completed task.
- When fixing a bug, explain what happened and *why* before changing anything.
