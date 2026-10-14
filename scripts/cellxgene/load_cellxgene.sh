#!/bin/bash
# driver.py's load step (LOAD_CMD): load_cellxgene.sh <release_uri>, run after publish and
# before releases/CURRENT is written. Exit 0 lets CURRENT flip; anything else makes the driver
# exit 6 with CURRENT untouched.
#
# BigQuery loaders and schemas live in genetics-results-db, so this only hands over to its
# scripts/load_cellxgene.sh, found through GENETICS_RESULTS_DB: set by the cellxgene image,
# which copies that loader in; elsewhere a checkout next to this repo.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB="${GENETICS_RESULTS_DB:-$(cd "${HERE}/../../.." && pwd)/genetics-results-db}"
LOADER="${DB}/scripts/load_cellxgene.sh"
if [ ! -x "${LOADER}" ]; then
  echo "load_cellxgene.sh: ${LOADER} not found; set GENETICS_RESULTS_DB to a genetics-results-db checkout" >&2
  exit 2
fi
exec "${LOADER}" "$@"
