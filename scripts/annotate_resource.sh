#!/bin/bash
# Stamp gnomAD consequence onto one served credible-set resource. The work is done by
# annotate_resource.py, whose header docstring is the reference; this wrapper only makes
# sure nothing lands on the root disk and the tools the run needs are there.

set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
python=${PYTHON:-python3}

staging=
dry_run=
prev=
for arg in "$@"; do
    case $prev in --staging) staging=$arg ;; esac
    case $arg in
        --staging=*) staging=${arg#--staging=} ;;
        --dry-run) dry_run=1 ;;
    esac
    prev=$arg
done

if [ -z "$dry_run" ]; then
    for tool in bgzip tabix; do
        command -v "$tool" >/dev/null || { echo "$0: $tool not on PATH" >&2; exit 2; }
    done
    if [ -n "$staging" ]; then
        # bgzip, sort and gcloud all honour TMPDIR, and the default is the root disk
        mkdir -p "$staging/tmp"
        export TMPDIR="$staging/tmp"
    fi
fi

exec "$python" "$here/annotate_resource.py" "$@"
