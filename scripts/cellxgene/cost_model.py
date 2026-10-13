#!/usr/bin/env python3
"""Initial-run and weekly-refresh cost model for the CELLxGENE epic (throwaway).

  cost_model.py PORTAL_DATASETS_JSON OUT_ROOT LOGS_DIR

Throughput comes from the proving_ground.py metrics under OUT_ROOT and the
download logs; prices are the Cloud Billing catalog SKUs listed in PRICES
(read 2026-10-13, USD, us-central1). Boot time is assumed, not measured.
"""
import datetime as dt
import glob
import json
import os
import sys

import numpy as np

# SKU id -> price; read from cloudbilling.googleapis.com/v1/services/*/skus on 2026-10-13
PRICES = {
    "e2_core": 0.02181159, "e2_ram": 0.00292353,            # CF4E-A0C7-E3BF, F449-33EC-A5EF
    "e2_core_spot": 0.01309, "e2_ram_spot": 0.001754,       # F179-E1EA-D97A, 9B1F-1E62-4061
    "n2_core": 0.031611, "n2_ram": 0.004237,                # BB77-5FDA-69D9, 5B01-D157-A097
    "n2_core_spot": 0.01896, "n2_ram_spot": 0.002542,       # 2856-91B1-406A, 6276-ABCE-E764
    "pd_balanced_gib_mo": 0.10,                             # 6AE1-525F-8B80
    "coldline_gib_mo": 0.004,                               # D716-3E59-667B Coldline Storage Iowa
    "standard_gib_mo": 0.02,                                # E5F0-6A5D-7BAD Standard Storage US Regional (first 5 GiB free)
    "bq_active_logical_gib_mo": 0.023,                      # 0018-A5A0-9D6D (us-central1, first 10 GiB free)
}
MACHINES = {  # name: (family, vCPU, GiB)
    "e2-standard-2": ("e2", 2, 8), "e2-highmem-4": ("e2", 4, 32), "n2-highmem-8": ("n2", 8, 64)}
BOOT_S = 150  # Batch VM provision + image pull: assumption until the proving ground (i8om.1) measures it
HOURS_PER_MONTH = 730
GiB = 2 ** 30


def size_class(filesize):
    if filesize <= 2e9:
        return "e2-standard-2"
    if filesize <= 20e9:
        return "e2-highmem-4"
    return "n2-highmem-8"


def hourly(machine, spot):
    fam, cpu, ram = MACHINES[machine]
    sfx = "_spot" if spot else ""
    return cpu * PRICES[f"{fam}_core{sfx}"] + ram * PRICES[f"{fam}_ram{sfx}"]


def main(listing, out_root, logs):
    mets = [json.load(open(p)) for p in glob.glob(f"{out_root}/*/metrics.json")]
    dls = [json.load(open(p)) for p in glob.glob(f"{logs}/*.download.json") if os.path.getsize(p)]
    big = [m for m in mets if m["nnz"] and m["nnz"] > 1e8]
    # per-nnz cost of every nnz-proportional phase, from the larger runs (small ones are overhead-dominated)
    def spn(m):
        return sum(p["s"] for k, p in m["phases"].items() if k in ("main_pass", "pass1_libsize")) / m["nnz"]
    # the multi-GB files ran ~30% slower per nnz than the small ones (larger scatter targets), so each class uses its own rate
    s_per_nnz = np.median([spn(m) for m in big if m["file_bytes"] <= 2e9])
    s_per_nnz_large = max(spn(m) for m in big if m["file_bytes"] > 2e9)
    s_per_cell = np.median([m["phases"]["obs_var"]["s"] / m["n_cells"] for m in big])
    write_s = max(m["phases"]["write"]["s"] for m in mets)
    big_dl = [d for d in dls if d["bytes"] > 1e9]
    bw = np.median([d["bytes"] / d["download_s"] for d in big_dl])
    bw_low = min(d["bytes"] / d["download_s"] for d in big_dl)
    out_per_cell = np.median([sum(os.path.getsize(f) for f in glob.glob(f"{out_root}/{m['dataset_id']}/*")) / m["n_cells"] for m in big])
    rows_per_group = np.median([m["group_gene_rows"] / m["n_groups"] for m in mets])
    print(f"model inputs: {1/s_per_nnz/1e6:.1f} M nnz/s (<=2 GB files), {1/s_per_nnz_large/1e6:.1f} M nnz/s (>2 GB, slowest measured), obs {s_per_cell*1e6:.1f} us/cell, download median {bw/1e6:.0f} MB/s "
          f"(min {bw_low/1e6:.0f}), boot {BOOT_S}s assumed, rows/group {rows_per_group:.0f}")

    d = json.load(open(listing))
    human = [x for x in d if any(o["ontology_term_id"] == "NCBITaxon:9606" for o in x["organism"])]
    work = [x for x in human if (x.get("primary_cell_count") or 0) > 0
            and not any("snmC" in a["label"] or "methyl" in a["label"].lower() for a in x["assay"])]
    now = dt.datetime(2026, 10, 13, tzinfo=dt.timezone.utc)

    def est(x, bandwidth):
        fs = next(a["filesize"] for a in x["assets"] if a["filetype"] == "H5AD")
        nnz = (x["cell_count"] or 0) * (x.get("mean_genes_per_cell") or 0)
        t = BOOT_S + fs / bandwidth + x["cell_count"] * s_per_cell + nnz * (s_per_nnz if fs <= 2e9 else s_per_nnz_large) + write_s
        return fs, nnz, t

    def run(sel, bandwidth, label):
        tot = {}
        crit = []
        for x in sel:
            fs, nnz, t = est(x, bandwidth)
            m = size_class(fs)
            disk_gib = max(30, fs * 1.1 / GiB + 20)
            r = tot.setdefault(m, {"n": 0, "h": 0.0, "gb": 0.0, "spot": 0.0, "od": 0.0, "disk": 0.0})
            r["n"] += 1; r["h"] += t / 3600; r["gb"] += fs / 1e9
            r["spot"] += t / 3600 * hourly(m, True)
            r["od"] += t / 3600 * hourly(m, False)
            r["disk"] += disk_gib * PRICES["pd_balanced_gib_mo"] / HOURS_PER_MONTH * t / 3600
            crit.append((t, fs, nnz, x["dataset_id"], x["title"][:50]))
        crit.sort(reverse=True)
        print(f"\n[{label}] {len(sel)} datasets, {sum(r['gb'] for r in tot.values())/1e3:.2f} TB")
        print("| class | n | GB | worker-h | spot $ | on-demand $ | pd-balanced $ |")
        print("|---|---|---|---|---|---|---|")
        for m, r in sorted(tot.items()):
            print(f"| {m} | {r['n']} | {r['gb']:.0f} | {r['h']:.1f} | {r['spot']:.2f} | {r['od']:.2f} | {r['disk']:.2f} |")
        s = {k: sum(r[k] for r in tot.values()) for k in ("n", "h", "spot", "od", "disk", "gb")}
        print(f"| total | {s['n']} | {s['gb']:.0f} | {s['h']:.1f} | {s['spot']:.2f} | {s['od']:.2f} | {s['disk']:.2f} |")
        for t, fs, nnz, i, title in crit[:3]:
            print(f"  critical path: {t/3600:.2f} h  {fs/1e9:.1f} GB  {nnz/1e9:.1f}B nnz  {i}  {title}")
        return s

    for bwl, b in (("median download", bw), ("slowest measured download", bw_low)):
        s = run(work, b, f"initial run, {bwl}")
    tb = sum(next(a["filesize"] for a in x["assets"] if a["filetype"] == "H5AD") for x in work)
    cells = sum(x["cell_count"] for x in work)
    agg_gib = cells * out_per_cell / GiB
    def card(x):
        return len(x["cell_type"]) * max(1, len(x["tissue"])) * max(1, len(x["disease"])) * max(1, len(x["assay"]))
    by_id = {x["dataset_id"]: x for x in work}
    calib = np.median([m["n_groups"] / card(by_id[m["dataset_id"]]) for m in mets if m["dataset_id"] in by_id])
    n_groups = calib * sum(card(x) for x in work)
    print(f"groups per dataset: listing-cardinality product x {calib:.2f} (median measured ratio)")
    bq_rows = n_groups * rows_per_group
    rows_per_cell = np.median([m["group_gene_rows"] / m["n_cells"] for m in big])
    bq_rows_low = rows_per_cell * cells
    print(f"BigQuery rows, lower estimate: {rows_per_cell:.1f} rows/cell (median measured) x {cells/1e6:.0f} M cells = {bq_rows_low/1e9:.2f}B rows "
          f"~{bq_rows_low*200/GiB:.0f} GiB -> ${max(0, bq_rows_low*200/GiB-10)*PRICES['bq_active_logical_gib_mo']:.2f}/month")
    bq_gib = bq_rows * 200 / GiB  # ~200 logical bytes/row: ids and labels as strings plus five numerics
    print(f"\nstorage/month: Coldline raw {tb/GiB:.0f} GiB -> ${tb/GiB*PRICES['coldline_gib_mo']:.2f}; "
          f"aggregates {agg_gib:.0f} GiB Standard -> ${max(0, agg_gib-5)*PRICES['standard_gib_mo']:.2f}; "
          f"BigQuery ~{bq_rows/1e9:.2f}B rows ~{bq_gib:.0f} GiB logical "
          f"-> ${max(0, bq_gib-10)*PRICES['bq_active_logical_gib_mo']:.2f}")

    def last(x):
        return max(dt.datetime.fromisoformat(v) for v in (x.get("published_at"), x.get("revised_at")) if v)
    wk13 = [x for x in work if last(x) >= now - dt.timedelta(weeks=13)]
    w = run(wk13, bw, "13-week churn (divide by 13 for a week)")
    print(f"weekly refresh: {len(wk13)/13:.1f} datasets/week, ${w['spot']/13:.2f} spot compute + ${w['disk']/13:.2f} disk per week "
          f"(on-demand ${w['od']/13:.2f})")


if __name__ == "__main__":
    main(*sys.argv[1:])
