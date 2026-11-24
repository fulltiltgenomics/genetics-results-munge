#!/usr/bin/env sh
# Warns when a commit changes something the docs describe but leaves the doc
# untouched. Mappings mirror the "Documentation ownership" table in CLAUDE.md.
#
# The rules are this repo's own; `check` and the rest of the mechanism live in
# scripts/doc-drift-lib.sh, which is the same file in every repo of the suite.

. "$(dirname "$0")/doc-drift-lib.sh"

DOC_README='^README\.md$'

check '^scripts/munge_[^/]*\.(py|sh)$' "$DOC_README" \
    'scripts/munge_*.{py,sh} -> README.md (the per-script dataset list under "other datasets", --stage and input/output flags)'

check '^scripts/create_[^/]*\.(py|sh)$' "$DOC_README" \
    'scripts/create_*.{py,sh} -> README.md (output products, the per-script run instructions and their arguments)'

check '^scripts/(sumstat|peak|exome)_utils\.py$' '^CLAUDE\.md$' \
    'scripts/{sumstat,peak,exome}_utils.py -> CLAUDE.md (shared helper list, required columns, index mode, tabix and GCS output rules)'

check '^scripts/peak_utils\.py$' "$DOC_README" \
    'scripts/peak_utils.py -> README.md (the open-chromatin dataset entry names it as the peak family shared writer)'

check '^scripts/coloc/[^/]*\.(py|sh)$' '^scripts/coloc/R14_UPDATE\.md$' \
    'scripts/coloc/ -> scripts/coloc/R14_UPDATE.md (input layout, metadata files, the invocation runbook)'

check '^scripts/munge_gene_disease\.(py|sh)$' '^docs/gene-disease-associations\.md$' \
    'scripts/munge_gene_disease.{py,sh} -> docs/gene-disease-associations.md (the source files taken, the version-in-the-filename cutover, the drop and dedup rules)'

check '^scripts/(munge_rcnv|rcnv_liftover_windows)\.(py|sh)$' '^docs/rcnv-[^/]*\.md$' \
    'scripts/munge_rcnv.{py,sh}, scripts/rcnv_liftover_windows.py -> docs/rcnv-dosage-sensitivity.md, docs/rcnv-sliding-windows.md (per-product column mappings, drop rules, the liftOver measurement)'

check '^wdl/(munge_finngen_finemapping_results|qtl_file)' "$DOC_README" \
    'wdl/munge_finngen_finemapping_results*, wdl/qtl_file.wdl -> README.md (WDL pipeline inputs, the Cromwell examples, output columns)'

check '^scripts/(build_gene_sets\.(py|sh)|gene_set_bundle\.py)$' '^docs/gene-sets\.md$' \
    'scripts/build_gene_sets.{py,sh}, scripts/gene_set_bundle.py -> docs/gene-sets.md (collections taken and refused, the licence table, id mapping rules, bundle layout, staging)'

check '^scripts/protein_structure_contract\.py$' '^docs/protein-structure\.md$' \
    'scripts/protein_structure_contract.py -> docs/protein-structure.md (release layout, model file id, window ownership and pair rules, PAE encoding, variant absence, view-spec)'

check '^scripts/(munge_decode_pqtl|decode_pqtl_phenotypes)\.py$' '^docs/pseudo-credible-sets\.md$' \
    'scripts/munge_decode_pqtl.py, scripts/decode_pqtl_phenotypes.py -> docs/pseudo-credible-sets.md (step 0 of the deCODE run: per-aptamer file columns, input-array selection, trait names)'

check '^wdl/(create_pseudo_credible_sets|autoreporting)' '^docs/pseudo-credible-sets\.md$' \
    'wdl/create_pseudo_credible_sets*, wdl/autoreporting_*.json -> docs/pseudo-credible-sets.md (thresholds, script defaults vs production flags, output columns, per-dataset table)'

check '^scripts/r14_meta_phenotypes\.py$' '^docs/r14-meta-serving\.md$' \
    'scripts/r14_meta_phenotypes.py -> docs/r14-meta-serving.md (which keys are omitted for quantitative traits, the contributing_cohorts naming, header-derived cohort presence and the TSV/header disagreement result)'

doc_drift_finish
