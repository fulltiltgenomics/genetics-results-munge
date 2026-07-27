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

1. Project documentation is maintained in files in `docs/` folder.
2. `docs/project-spec.md` is an overview of project purpose, structure and logic.
3. Create other files under `docs/` if necessary.
4. Maintain `docs/project-spec.md` and any other generated files to be up to date with the project.
5. Reread `docs/project-spec.md` often and whenever you need to refresh your context with what the project is about and implementation logic.
6. This should often be your first step in understanding a task.


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
- use `write_sumstat_output()` from `sumstat_utils.py` for the write+filter+upload pattern

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
- `read_gnomad_filtered()` — read filtered gnomAD TSV (plain or gzipped)
- `build_rsid_set()` — extract rsid set from DataFrame
- `stream_gnomad_by_rsid()` — stream gnomAD keeping rsid matches

## gnomAD allele alignment
- PGC/BIP: join on rsid, classify A1/A2 vs ref/alt orientation, flip beta and AFs when A1 is ref
- GP2: join on chr:pos:ref:alt (both build 38), flip SNPs where effect_allele is ref
- IBD: join on chr:pos for AF-AF plot only (no allele flipping in output)


# Exome / Gene Burden Output Formats

Column definitions live with the scripts that write them: gene burden results in the final
`select()` of `scripts/genebass/convert_genebass_gene_results.py`, exome variant results in the
header comment of `scripts/genebass/convert_genebass_variant_results.sh`.


# API Integration

Sumstat column mappings are defined in:
- `genetics-results-api/app/config/profiles/finngen/summary_stats.py`
- `genetics-results-api/app/config/profiles/daly/summary_stats.py`

The API maps file columns to internal names (e.g. `INFO` → `info`, `neff` → `n_eff`). Update both profile files when adding new sumstat datasets or changing column names.


====

**Don't forget any of the 'QUALITY CODING RULES' above!!!**
