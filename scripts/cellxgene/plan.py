#!/usr/bin/env python3
"""CELLxGENE planner: which datasets still need aggregating, and on what machine.

    plan.py --out plan.json [--run-id ID] [--approve-large]
            [--max-datasets 200] [--max-gb 300] [--listing listing.json]

Reads the curation API listing, keeps what `aggregate.py` would process, and subtracts the
datasets whose `_SUCCESS.json` already exists under
`gs://<bucket>/cellxgene/agg/v<FORMAT_VERSION>/<dataset_version_id>/`. The completed set
comes from one prefix listing, not a probe per dataset. `driver.py` consumes the JSON.

In scope: organism Homo sapiens (NCBITaxon:9606), `primary_cell_count` > 0 (a dataset with no primary cells
never enters the plan, so it is not downloaded), an H5AD asset, and at least one assay in
`aggregate.ASSAYS` — the worker's own allow-list, imported so the two cannot drift. The
datasets that merely contain an allowed assay alongside a disallowed one are in scope; the
worker drops the other assays' cells.

Machine class (one class, one exception): `e2-standard-2` (8 GB). `e2-highmem-2` (16 GB) when
`est_nnz` >= HIGHMEM_NNZ or `cell_count` >= HIGHMEM_CELLS. `est_nnz` = cell_count x
mean_genes_per_cell from the listing. The worker caps its accumulators at 4 GB however large
the dataset; beside them it holds per-cell obs arrays (roughly 8 bytes per cell each, over
ALL cells), so `cell_count` is the threshold with a memory basis: the heaviest dataset (4.1 M
cells, 15.6 B nnz) peaked at 5.4 GB with a 3 GB accumulator, and 3 M cells is about three
quarters of it. nnz sets runtime, not memory (chunks are bounded by CHUNK_NNZ); its threshold
is a conservative heuristic that moves the long runs onto the bigger machine for headroom.
Items are ordered largest `est_nnz` first.

Plan JSON
    run_id, aggregation_format_version, generated_at, listing_url, agg_prefix,
    items[]: dataset_id, dataset_version_id, collection_id, title, url, filesize (bytes),
        cell_count, primary_cell_count, est_nnz, machine_class,
        assays[]: {ontology_term_id, label} of the allowed assays present
    superseded[]: dataset_id, dataset_version_id (done, now stale), current_dataset_version_id
    totals: listed, in_scope, done, todo, todo_bytes, todo_gb, by_machine_class
        {class: {datasets, gb}}, superseded, skipped_no_primary
    breaker: tripped, max_datasets, max_gb, approved

Exit codes: 0 plan written and cleared to run, 2 error, 4 plan exceeds the circuit breaker
and --approve-large was not given (the plan file is still written; the driver must not run).
"""
import argparse
import datetime as dt
import json
import os
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggregate  # noqa: E402

API_URL = "https://api.cellxgene.cziscience.com/curation/v1/datasets"
BUCKET = "daly-expression-data"
AGG_ROOT = "cellxgene/agg"
STANDARD, HIGHMEM = "e2-standard-2", "e2-highmem-2"
HIGHMEM_NNZ = 5e9
HIGHMEM_CELLS = 3_000_000
BREAKER_EXIT = 4
MARKER = "_SUCCESS.json"
HUMAN_TAXON = "NCBITaxon:9606"


def fetch_listing(url=API_URL):
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.load(r)


def agg_prefix(fmt=aggregate.FORMAT_VERSION):
    return f"{AGG_ROOT}/v{fmt}/"


def _client():
    from google.cloud import storage
    try:
        return storage.Client()
    except OSError:
        # object reads and prefix listings bill nothing to a project, but the client
        # refuses to construct without one when none can be inferred from the environment
        return storage.Client(project="_")


def list_done(bucket=BUCKET, fmt=aggregate.FORMAT_VERSION):
    """dataset_version_ids with a completion marker, from a single prefix listing."""
    prefix = agg_prefix(fmt)
    blobs = _client().list_blobs(bucket, prefix=prefix)
    return {n.split("/")[-2] for n in (b.name for b in blobs)
            if n.endswith("/" + MARKER) and n.count("/") == prefix.count("/") + 1}


def read_marker_dataset_id(version_id, bucket=BUCKET, fmt=aggregate.FORMAT_VERSION, client=None):
    blob = (client or _client()).bucket(bucket).blob(f"{agg_prefix(fmt)}{version_id}/{MARKER}")
    return json.loads(blob.download_as_bytes()).get("dataset_id")


def allowed_assays(ds):
    return [{"ontology_term_id": a["ontology_term_id"], "label": a.get("label")}
            for a in ds.get("assay") or [] if a["ontology_term_id"] in aggregate.ASSAYS]


def h5ad_asset(ds):
    return next((a for a in ds.get("assets") or [] if a.get("filetype") == "H5AD"), None)


def classify(ds):
    """(item-or-None, reason). Reason is 'no_primary' only for the case counted in totals."""
    if not any(o.get("ontology_term_id") == HUMAN_TAXON for o in ds.get("organism") or []):
        return None, "organism"
    if not (ds.get("primary_cell_count") or 0) > 0:
        return None, "no_primary"
    assays, asset = allowed_assays(ds), h5ad_asset(ds)
    if not assays or not asset:
        return None, "assay" if not assays else "no_h5ad"
    cells = ds.get("cell_count") or 0
    nnz = round(cells * (ds.get("mean_genes_per_cell") or 0))
    return {"dataset_id": ds["dataset_id"], "dataset_version_id": ds["dataset_version_id"],
            "collection_id": ds["collection_id"], "title": ds.get("title"),
            "url": asset["url"], "filesize": asset["filesize"], "cell_count": cells,
            "primary_cell_count": ds["primary_cell_count"], "est_nnz": nnz,
            "machine_class": machine_class(nnz, cells), "assays": assays}, None


def machine_class(est_nnz, cell_count, nnz_limit=HIGHMEM_NNZ, cells_limit=HIGHMEM_CELLS):
    return HIGHMEM if est_nnz >= nnz_limit or cell_count >= cells_limit else STANDARD


def build_plan(listing, done, superseded_dataset_ids=None, run_id=None, max_datasets=200,
               max_gb=300.0, approved=False, fmt=aggregate.FORMAT_VERSION):
    """superseded_dataset_ids: done version id -> dataset_id, for done ids not in the listing."""
    in_scope, skipped = [], {}
    for ds in listing:
        item, why = classify(ds)
        if item:
            in_scope.append(item)
        else:
            skipped[why] = skipped.get(why, 0) + 1
    items = sorted((i for i in in_scope if i["dataset_version_id"] not in done),
                   key=lambda i: (-i["est_nnz"], i["dataset_version_id"]))
    current = {ds["dataset_id"]: ds["dataset_version_id"] for ds in listing}
    superseded = sorted(
        ({"dataset_id": did, "dataset_version_id": vid, "current_dataset_version_id": current[did]}
         for vid, did in (superseded_dataset_ids or {}).items()
         if did in current and current[did] != vid),
        key=lambda s: s["dataset_version_id"])
    by_class = {}
    for i in items:
        c = by_class.setdefault(i["machine_class"], {"datasets": 0, "gb": 0.0})
        c["datasets"] += 1
        c["gb"] = round(c["gb"] + i["filesize"] / 1e9, 3)
    todo_bytes = sum(i["filesize"] for i in items)
    tripped = len(items) > max_datasets or todo_bytes / 1e9 > max_gb
    return {
        "run_id": run_id or dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
        "aggregation_format_version": fmt,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "listing_url": API_URL, "agg_prefix": f"gs://{BUCKET}/{agg_prefix(fmt)}",
        "items": items, "superseded": superseded,
        "totals": {"listed": len(listing), "in_scope": len(in_scope),
                   "done": len(in_scope) - len(items), "todo": len(items),
                   "todo_bytes": todo_bytes, "todo_gb": round(todo_bytes / 1e9, 3),
                   "by_machine_class": by_class, "superseded": len(superseded),
                   "skipped_no_primary": skipped.get("no_primary", 0)},
        "breaker": {"tripped": tripped, "max_datasets": max_datasets, "max_gb": max_gb,
                    "approved": approved},
    }


def summary(plan):
    t, b = plan["totals"], plan["breaker"]
    lines = [f"run {plan['run_id']}: {t['todo']} to do of {t['in_scope']} in scope "
             f"({t['done']} done), {t['todo_gb']:.1f} GB, {t['superseded']} superseded"]
    lines += [f"  {c}: {v['datasets']} datasets, {v['gb']:.1f} GB"
              for c, v in sorted(t["by_machine_class"].items())]
    if b["tripped"] and not b["approved"]:
        lines.append(f"PLAN TOO LARGE (limit {b['max_datasets']} datasets / {b['max_gb']:g} GB): "
                     "not cleared to run; review the plan and rerun with --approve-large")
    return "\n".join(lines)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--out", required=True)
    p.add_argument("--run-id")
    p.add_argument("--bucket", default=BUCKET)
    p.add_argument("--listing", help="read this saved API listing instead of calling the API")
    p.add_argument("--max-datasets", type=int, default=200)
    p.add_argument("--max-gb", type=float, default=300.0)
    p.add_argument("--approve-large", action="store_true")
    a = p.parse_args(argv)
    try:
        if a.listing:
            with open(a.listing) as fh:
                listing = json.load(fh)
        else:
            listing = fetch_listing()
        done = list_done(a.bucket)
        listed = {ds["dataset_version_id"] for ds in listing}
        client = _client() if done - listed else None
        stale = {v: read_marker_dataset_id(v, a.bucket, client=client) for v in sorted(done - listed)}
        plan = build_plan(listing, done, stale, a.run_id, a.max_datasets, a.max_gb, a.approve_large)
    except Exception as e:
        print(f"plan.py: {type(e).__name__}: {e}", file=sys.stderr)
        return 2
    with open(a.out, "w") as fh:
        json.dump(plan, fh, indent=1)
    print(summary(plan))
    return BREAKER_EXIT if plan["breaker"]["tripped"] and not a.approve_large else 0


if __name__ == "__main__":
    sys.exit(main())
