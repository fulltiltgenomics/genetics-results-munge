====


QUALITY CODING RULES


# Code changes

1. If you find errors or suggestions in code which are not DIRECTLY related to user's current request, never change it without asking first.
2. Before suggesting changes to files, always assume user might have changed the file since your last read and consider reading the file again.


# Security

1. Never commit sensitive files (.env, credentials, API keys)
2. Use environment variables for API keys and credentials
3. Keep API keys and credentials out of logs and output


# Project Specifications

1. `README.md` is the overview of project purpose, structure and logic. Reading it should often be your first step in understanding a task.
2. Longer-form documentation lives in the `docs/` folder. Currently that is `docs/pseudo-credible-sets.md`, which defines how pseudo credible sets are derived from external summary statistics.
3. Create other files under `docs/` if necessary.
4. Maintain `README.md` and everything under `docs/` to be up to date with the project.
5. Per-script details (source files, format assumptions, command line flags) belong in the script's own header comment, not in the README — the README points at them.


# Documentation ownership

Changing a path on the left makes the doc on the right wrong until it is updated in
the same commit. `scripts/check-doc-drift.sh` warns (never blocks) on commits that
violate this; it runs from the `pre-commit` hook.

| changed path | doc to update | what to check |
|---|---|---|
| `scripts/munge_*.py`, `scripts/munge_*.sh` | `README.md` | the per-script dataset list under "other datasets", `--stage` and input/output flags |
| `scripts/create_*.py`, `scripts/create_*.sh` | `README.md` | output products, the per-script run instructions and their arguments |
| `scripts/sumstat_utils.py` | `CLAUDE.md` | shared helper list, required sumstat columns, tabix and GCS output rules |
| `scripts/coloc/*.py`, `scripts/coloc/*.sh` | `scripts/coloc/R14_UPDATE.md` | input layout, metadata files, the invocation runbook |
| `wdl/munge_finngen_finemapping_results*`, `wdl/qtl_file.wdl` | `README.md` | WDL pipeline inputs, the Cromwell examples, output columns |
| `wdl/create_pseudo_credible_sets*`, `wdl/autoreporting_*.json` | `docs/pseudo-credible-sets.md` | thresholds, script defaults vs production flags, output columns, per-dataset table |

A doc is stale the moment it *enumerates* something the code no longer matches.
Counts and lists rot silently — dataset lists, column tables, helper-function lists,
per-dataset invocation tables — so re-derive them from the code rather than trusting them.

The hook lives in `.beads/hooks/pre-commit`, which **is** version controlled, so a fresh
clone gets it as soon as `core.hooksPath` points at `.beads/hooks` (`bd init` sets that;
`bd` re-points it on checkout). The doc-drift call sits after the
`# --- END BEADS INTEGRATION ... ---` marker, never inside the beads-managed fence, which
beads rewrites. It only needs to run `scripts/check-doc-drift.sh || true`.


# Cross-repo documentation

`genetics-results-suite` is the spec of record for the suite as a whole; this repo
documents only its own munging. Changing a munge output that feeds a BigQuery table
means `genetics-results-db` (loaders, schemas) and that repo's docs may need updating
too — check them in the same session rather than assuming they follow.


# Software Development Behavior Guidelines

1. Don't guess and do things which you are not certain about. Ask the user instead.
2. Don't add or modify code unrelated to the specific request and context at the moment.
3. In interactive mode: only use git when asked, stage changes and propose a commit message for user review. In autonomous/orchestrator mode (e.g. ralph wiggum): commit after each completed task with a descriptive message.
4. **Always** prior to finishing a task and considering it completed, revise all the changes and update Project Specification files.
5. When trying to fix any bug or error **ALWAYS** think carefully and analyze in detail what happened and WHY? Explain and confirm with user.


# Code Conventions

1. Project structure:
   - `scripts/` - Munging scripts (one per dataset)
   - `scripts/sumstat_utils.py` - Shared utilities for sumstat munging
   - `scripts/coloc/`, `scripts/genebass/` - Munging of colocalization and Genebass results
   - `wdl/` - WDL pipelines and their per-dataset input and Cromwell option JSONs
   - `metadata/` - Dataset and study metadata read by the scripts
   - `docs/` - Project documentation
2. Code should be self-descriptive
   - Only add comments for tricky or complex parts of the code (explaining WHY something is done)
   - NO redundant and trivial comments that simply restate what the code does
3. Private fields and methods should be prefixed with underscore
4. Git commit messages should be concise and descriptive


# Sumstat Munging Rules

## Output format
- bgzipped TSV (`.munged.tsv.gz`) with tabix index (`.tsv.gz.tbi`)
- tabix indexed with `-s1 -b2 -e2` (chr in col 1, pos in col 2)
- always produce a filtered companion file `.mlog10p_gt4.tsv.gz` (mlog10p > 4) with its own tabix index
- use `write_sumstat_output()` from `sumstat_utils.py` for the write+filter+upload pattern. `munge_ibd.py` is the exception: it streams chromosome by chromosome into two long-lived bgzip pipes instead of holding a whole DataFrame, so it indexes and uploads itself

## Required columns (all sumstat files must have these)
- `#chr` — chromosome as integer, **X chromosome is always 23**
- `pos` — genomic position, GRCh38
- `ref` — reference allele
- `alt` — alternative (effect) allele
- `beta` — effect size (log-OR for binary traits), formatted as `:.3e`
- `se` — standard error, formatted as `:.3e`
- `mlog10p` — **always use -log10(p)**, never raw p-values. Round to 4 decimals. When p underflows to 0, compute from beta/se using `log_ndtr`
- `af` — allele frequency (from gnomAD or study)

## Column naming conventions
- `INFO` — imputation quality, always **uppercase**
- `neff` — effective sample size, always **lowercase**
- `n_cases` / `n_controls` — case/control counts (not `nca`/`nco`, not `Nca`/`Nco`)
- `af_cases` / `af_controls` — allele frequencies in cases/controls
- `rsid` — dbSNP rs identifier

## Chromosome handling
- read chr as string from input, map `X` → `23`, cast to Int32
- drop rows where chr is null after casting (malformed lines)
- IBD script iterates per-chromosome files and uses `CHR_MAP = {"X": 23}`

## GCS output
- all scripts accept `gs://` paths in `--output`
- write locally to temp dir, tabix-index, then upload both `.tsv.gz` and `.tsv.gz.tbi` via `gcloud storage cp`
- use `upload_to_gcs()` from `sumstat_utils.py`

## Shared code (`scripts/sumstat_utils.py`)
- gnomAD constants: `GNOMAD_DEFAULT`, `GNOMAD_AF_COLS`, `GNOMAD_KEEP_COLS`
- `upload_to_gcs()` — upload file + .tbi to GCS
- `write_bgzip()` — bgzip + tabix a DataFrame
- `write_sumstat_output()` — full + filtered write with GCS support
- `write_exome_output()` — the same for exome results, with the tabix columns and the mlog10p column name given by the caller (gene burden files are indexed on the gene locus, `-s5 -b6 -e6`)
- `read_gnomad_filtered()` — read filtered gnomAD TSV (plain or gzipped)
- `build_rsid_set()` — extract rsid set from DataFrame
- `stream_gnomad_by_rsid()` — stream gnomAD keeping rsid matches

## gnomAD allele alignment
- PGC/BIP: join on rsid, classify A1/A2 vs ref/alt orientation, flip beta and AFs when A1 is ref, drop mismatches, and take chr/pos from gnomAD (build 38)
- GP2: join on chr:pos:ref:alt (both build 38), flip beta and AF for SNPs where effect_allele is ref; never flip indels, whose representation may differ
- IBD: gnomAD is streamed by (chr, pos) and joined on chr:pos:ref:alt in both orientations, for the AF-AF plot only — the output is never flipped
- COVID: ALT is already the effect allele, so the chr:pos:ref:alt join only adds rsid and the `most_severe` / `gene_most_severe` annotation
- `--gnomad-filtered` on every one of these skips the streaming pass and reads a previously saved filtered gnomAD TSV instead; that file is what the drivers reuse across phenotypes of the same study


# Exome / Gene Burden Output Formats

Column definitions live with the scripts that write them: gene burden results in the final
`select()` of `scripts/genebass/convert_genebass_gene_results.py`, exome variant results in the
header comment of `scripts/genebass/convert_genebass_variant_results.sh`.

The non-Genebass exome munges (`munge_schema*.py`, `munge_bipex.py`, `munge_ibd_exome.py`,
`munge_ibd_supp_*.py`) reproduce those two layouts in their own `build_output()` and write them
with `write_exome_output()`, filling columns the source does not provide with NA. Gene burden
files are tabix indexed on the gene locus (`-s5 -b6 -e6`), variant files on the variant position
(`-s2 -b3 -e3`).

`munge_als.py` and `munge_asmqtl.py` are variant-level but not in that layout — ALS follows
`scripts/genebass/cleanup_genebass_variant_results.py` with `most_severe`/`gene_most_severe`
added, and ASM-QTL has its own methylation-target columns. Check their `select()` before
assuming a shared schema.


# API Integration

Sumstat column mappings are defined in:
- `genetics-results-api/app/config/profiles/finngen/summary_stats.py`
- `genetics-results-api/app/config/profiles/daly/summary_stats.py`

The API maps file columns to internal names (e.g. `INFO` → `info`, `neff` → `n_eff`). Update both profile files when adding new sumstat datasets or changing column names.


====

**Don't forget any of the 'QUALITY CODING RULES' above!!!**


<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:ca08a54f -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

## Session Completion

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   bd dolt push
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
<!-- END BEADS INTEGRATION -->
