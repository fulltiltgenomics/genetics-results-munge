#!/usr/bin/env sh
# Warns when a commit changes something the docs describe but leaves the doc
# untouched. Mappings mirror the "Documentation ownership" table in CLAUDE.md.
#
# This never blocks. A warning that is occasionally ignored beats a gate that
# gets bypassed with --no-verify, because a bypassed gate is both absent and
# assumed present.

set -u

root=$(git rev-parse --show-toplevel 2>/dev/null) || exit 0
cd "$root" || exit 0

staged=$(git diff --cached --name-only --diff-filter=ACMRD)
[ -n "$staged" ] || exit 0

hit() {
    printf '%s\n' "$staged" | grep -qE "$1"
}

found=0
check() {
    if hit "$1" && ! hit "$2"; then
        if [ "$found" -eq 0 ]; then
            printf '\ndoc-drift warning — this commit changes code the docs describe:\n\n' >&2
            found=1
        fi
        printf '  %s\n' "$3" >&2
    fi
}

DOC_README='^README\.md$'

check '^scripts/munge_[^/]*\.(py|sh)$' "$DOC_README" \
    'scripts/munge_*.{py,sh} -> README.md (the per-script dataset list under "other datasets", --stage and input/output flags)'

check '^scripts/create_[^/]*\.(py|sh)$' "$DOC_README" \
    'scripts/create_*.{py,sh} -> README.md (output products, the per-script run instructions and their arguments)'

check '^scripts/sumstat_utils\.py$' '^CLAUDE\.md$' \
    'scripts/sumstat_utils.py -> CLAUDE.md (shared helper list, required sumstat columns, tabix and GCS output rules)'

check '^scripts/coloc/[^/]*\.(py|sh)$' '^scripts/coloc/R14_UPDATE\.md$' \
    'scripts/coloc/ -> scripts/coloc/R14_UPDATE.md (input layout, metadata files, the invocation runbook)'

check '^wdl/(munge_finngen_finemapping_results|qtl_file)' "$DOC_README" \
    'wdl/munge_finngen_finemapping_results*, wdl/qtl_file.wdl -> README.md (WDL pipeline inputs, the Cromwell examples, output columns)'

check '^wdl/(create_pseudo_credible_sets|autoreporting)' '^docs/pseudo-credible-sets\.md$' \
    'wdl/create_pseudo_credible_sets*, wdl/autoreporting_*.json -> docs/pseudo-credible-sets.md (thresholds, script defaults vs production flags, output columns, per-dataset table)'

if [ "$found" -eq 1 ]; then
    printf '\n  Update the doc in this commit, or note why it does not apply.\n' >&2
    printf '  Not blocking. Mappings live in CLAUDE.md > Documentation ownership.\n\n' >&2
fi

exit 0
