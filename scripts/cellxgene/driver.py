#!/usr/bin/env python3
"""CELLxGENE driver: plan, fan the plan out to Cloud Batch, then publish, load and flip CURRENT.

    driver.py run  [--root cellxgene] [--image REGISTRY/cellxgene:TAG] [--approve-large]
                   [--only VID,...] [--limit N] [--allow-partial] [--no-publish]
                   [--publish-cmd TEMPLATE] [--load-cmd TEMPLATE] [--allow-missing-load]
                   [--allow-shrink]
    driver.py task                      # the Batch worker entrypoint; reads its item from env

Bucket layout, all under `gs://<bucket>/<root>/` (`--root cellxgene-staging` for rehearsals):
- `agg/v<FORMAT>/<dataset_version_id>/` worker outputs; `_SUCCESS.json` (done) and
  `_REJECTED.json` (done, not retried: aggregate.py exit 3) there are the ONLY resume state.
  `runs/` is deleted by object age, so nothing that decides work lives in it
- `raw/<dataset_version_id>.h5ad` the source file, uploaded as COLDLINE with
  ifGenerationMatch=0: overwriting a Coldline object younger than 90 days is billed as early
  deletion, and spot retries would otherwise re-upload it
- `runs/<run_id>/` plan and listing copies, per-job item files, per-dataset outcome records
  (`outcomes/<vid>.json`, written by the task), `report.json`
- `releases/<release>/` publish output, `<release>` = `<UTC date>-<HHMMSS>` of the run, unique
  per run; the driver refuses to publish into a release prefix that already has objects
- `releases/CURRENT` written last, one small JSON object, read by results-api (.12):
    release            the release name; `uri` is gs://<bucket>/<root>/releases/<release>/
    run_id             the driver run that wrote it (its report is runs/<run_id>/report.json)
    written_at         UTC ISO time of the flip
    datasets           how many dataset versions the release holds
    bigquery_loaded    false only when the load command was absent and --allow-missing-load
                       let the flip through: the tabix files are current, BigQuery is not
    publish_key        sha256 over `publish_inputs`; equal key = an identical release, so a
                       rerun with nothing new does not republish
    publish_inputs     datasets_sha256 (sorted version ids), aggregation_format_version,
                       rollup_sha256 (per pinned roll-up table), publish_sha256 (publish.py),
                       image (worker image digest), publish_cmd (template, carries thresholds)
- `driver.lock` the concurrency lock

Steps
1. Lock: `driver.lock` created with ifGenerationMatch=0, holding owner, host, run_id, time.
   Stale-lock policy: a lock whose object is older than `--lock-stale-hours` (72) is taken
   over (deleted with ifGenerationMatch=<its generation>, so of two drivers racing to break
   it only one wins); a younger one stops this driver with exit 75. 72 h exceeds the driver
   Batch job's own maxRunDuration (48 h, driver-batch-job.json), so a lock that old can only
   belong to a driver that is dead. A driver run by hand has no such cap: do not run one for
   longer than that.
2. Plan: `plan.py` on the listing (exit 4 = circuit breaker: stop unless --approve-large).
   `--only`/`--limit` restrict the LISTING before planning, so the breaker sees the restricted
   plan; a restricted run never publishes, since a release built from a partial listing
   would drop every other dataset.
3. Submit: one Batch job per (machine class, boot-disk tier); one task per dataset; task
   index -> item through `runs/<run_id>/jobs/<job>.items.json`. Batch accepts one lifecycle
   policy per task: FAIL_TASK on exit 3; everything else, spot preemption (50001) included,
   is retried up to maxRetryCount. Tasks of a job share VMs:
   `parallelism` caps the VMs and each runs tasks back to back, which matters because VM
   start (76 s measured) dwarfs the work on a small file. The image is pinned by digest.
4. Poll to completion. Past --timeout-hours, still-running jobs are cancelled. Cancelling is
   asynchronous, so the driver keeps the lock and polls them until each is terminal, for at
   most --cancel-wait-minutes (15); only then is the lock released. Jobs that are still not
   terminal at that bound are reported as `cancel_unconfirmed` and the lock is released
   anyway: their tasks MAY still be running and uploading, so a successor can race them
   (harmless to raw/ and markers, which are no-clobber and written last, but not ruled out).
   Per dataset: ok (marker
   present), rejected (worker exit 3: the task writes `_REJECTED.json` and deletes the raw
   copy it uploaded, which nothing will read), failed (anything else after retries). Failed
   datasets stop the publish unless --allow-partial; rejected ones do not.
5. Superseded: customTime=now on `raw/<old>.h5ad` only once the NEW version has its marker,
   and only if no customTime is set yet (it can only move forward; the lifecycle rule
   deletes 90 days after it). Datasets withdrawn from the listing are reported, not acted on.
   Decision for a new version that never succeeds (fails, or is rejected): the old version's
   outputs and raw file are kept, but the dataset drops out of the release, because publish
   reads only the listing's CURRENT versions; it is reported every run as
   `waiting_for_new_version` until the new version has a marker.
6. Publish runs in this process, not as another Batch job: the scheduled driver already IS a
   Batch job on a VM sized for it, so a second job would add a VM start, an image pull and a
   second poll for nothing. Then the load command; then `releases/CURRENT`. Shrink guard:
   if the release would hold < 90% of CURRENT's `datasets`, nothing is published, loaded or
   flipped unless --allow-shrink (checked before publish, so BigQuery is not loaded with a
   release that will not be served).
SIGTERM (Batch maxRunDuration, a cancelled job) exits through the same path, so the lock is
released.

Exit codes: 0 done, 2 error, 4 breaker, 5 failed datasets (nothing published),
6 publish or load failed, 7 shrink guard, 75 locked by another driver.
"""
import argparse
import datetime as dt
import json
import os
import hashlib
import shlex
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import aggregate  # noqa: E402
import plan  # noqa: E402

PROJECT, REGION = "daly-finngenie", "us-central1"
BUCKET, ROOT = plan.BUCKET, "cellxgene"
NETWORK = f"projects/{PROJECT}/global/networks/expression-batch-network"
SUBNET = f"projects/{PROJECT}/regions/{REGION}/subnetworks/expression-batch-subnet"
WORKER_SA = f"expression-batch-worker@{PROJECT}.iam.gserviceaccount.com"
IMAGE = f"{REGION}-docker.pkg.dev/{PROJECT}/genetics-results/cellxgene:latest"
TEMPLATE = HERE / "worker-batch-job.json"
REJECTED_EXIT = 3
EXIT_BREAKER, EXIT_FAILED, EXIT_PUBLISH, EXIT_SHRINK, EXIT_LOCKED = 4, 5, 6, 7, 75
SHRINK_FLOOR = 0.9
# memoryMib must fit what Batch leaves allocatable on the VM, or the job never schedules
MEMORY_MIB = {plan.STANDARD: 6144, plan.HIGHMEM: 14336}
LOCK_STALE_HOURS = 72
TERMINAL = {"SUCCEEDED", "FAILED", "CANCELLED", "DELETION_IN_PROGRESS"}
PUBLISH_CMD = (f"{sys.executable} {HERE}/publish.py --agg {{agg}} --out {{releases}} "
               "--release-date {release} --datasets {datasets} --listing {listing}")
LOAD_CMD = f"{HERE}/load_cellxgene.sh {{release_uri}}"


def now():
    return dt.datetime.now(dt.timezone.utc)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- sizing

def disk_gb(filesize):
    """2x the file (download + staged outputs) + 20 GB for the image and OS, in tiers so
    similar files share a job: 30, 50, then multiples of 50."""
    need = 2 * filesize / 1e9 + 20
    return 30 if need <= 30 else 50 if need <= 50 else -(-int(need) // 50) * 50


def task_timeout_s(filesize):
    # 6 min/GB is ~2x the slowest measured transfer (38 MB/s, ~26 s/GB each way) plus
    # aggregation at ~10 M nnz/s; the 1 h floor covers the VM start and image pull
    return int(min(3600 + 360 * filesize / 1e9, 24 * 3600))


def group_items(items):
    """(machine_class, disk_gb) -> items, largest file first within each group."""
    groups = {}
    for i in sorted(items, key=lambda i: -i["filesize"]):
        groups.setdefault((i["machine_class"], disk_gb(i["filesize"])), []).append(i)
    return groups


def job_id(run_id, machine_class, disk):
    return f"cxg-{run_id.lower().replace('_', '-')}-{machine_class.split('-')[1]}-{disk}g"[:63].rstrip("-")


def build_job(items_uri, items, machine_class, disk, image_digest, run_id, bucket, root,
              parallelism, template=TEMPLATE):
    job = json.loads(Path(template).read_text())
    tg = job["taskGroups"][0]
    tg["taskCount"] = len(items)
    tg["parallelism"] = min(len(items), parallelism)
    spec = tg["taskSpec"]
    spec["runnables"][0]["container"]["imageUri"] = image_digest
    spec["computeResource"]["memoryMib"] = MEMORY_MIB[machine_class]
    spec["maxRunDuration"] = f"{max(task_timeout_s(i['filesize']) for i in items)}s"
    spec["environment"] = {"variables": {
        "RUN_ID": run_id, "ITEMS_URI": items_uri, "BUCKET": bucket, "ROOT": root,
        "WORKER_IMAGE_DIGEST": image_digest}}
    pol = job["allocationPolicy"]["instances"][0]["policy"]
    pol["machineType"] = machine_class
    pol["bootDisk"]["sizeGb"] = str(disk)
    job["labels"]["run-id"] = run_id.lower()
    return job


# ---------------------------------------------------------------- GCS

def storage_client():
    from google.cloud import storage
    return storage.Client(project=PROJECT)


def precondition_failed():
    from google.api_core.exceptions import PreconditionFailed
    return PreconditionFailed


def not_found():
    from google.api_core.exceptions import NotFound
    return NotFound


class Layout:
    def __init__(self, bucket, root, fmt=aggregate.FORMAT_VERSION):
        self.bucket, self.root, self.fmt = bucket, root.strip("/"), fmt

    def agg(self, vid=""):
        return f"{self.root}/agg/v{self.fmt}/{vid}"

    def marker(self, vid):
        return f"{self.agg(vid)}/{plan.MARKER}"

    def rejected(self, vid):
        return f"{self.agg(vid)}/{plan.REJECTED}"

    # raw/ is keyed by dataset_version_id only, not by FORMAT_VERSION like the aggregates: after a
    # format bump, a stricter check that rejects a version deletes the raw of a version that
    # succeeded under the old format (recoverable while CZI still serves that version)
    def raw(self, vid):
        return f"{self.root}/raw/{vid}.h5ad"

    def run(self, run_id, rest=""):
        return f"{self.root}/runs/{run_id}/{rest}"

    def lock(self):
        return f"{self.root}/driver.lock"

    def current(self):
        return f"{self.root}/releases/CURRENT"

    def uri(self, path):
        return f"gs://{self.bucket}/{path}"


def put_json(bkt, path, obj, **kw):
    bkt.blob(path).upload_from_string(json.dumps(obj, indent=1), content_type="application/json", **kw)


def markers(client, lay):
    """(done, rejected) version ids."""
    m = plan.list_markers(lay.bucket, lay.fmt, lay.root, client)
    return m[plan.MARKER], m[plan.REJECTED]


class Locked(Exception):
    pass


def acquire_lock(bkt, lay, run_id, stale_hours=LOCK_STALE_HOURS, clock=now):
    """Return the lock's generation. Raises Locked when a live driver holds it."""
    body = {"owner": os.environ.get("BATCH_JOB_UID") or os.environ.get("USER") or "unknown",
            "host": socket.gethostname(), "pid": os.getpid(), "run_id": run_id,
            "acquired_at": clock().isoformat(timespec="seconds")}
    for _ in range(2):
        blob = bkt.blob(lay.lock())
        try:
            blob.upload_from_string(json.dumps(body), if_generation_match=0)
            return blob.generation
        except precondition_failed():
            pass
        held = bkt.get_blob(lay.lock())
        if held is None:
            continue
        age_h = (clock() - held.time_created).total_seconds() / 3600
        if age_h < stale_hours:
            raise Locked(f"{lay.uri(lay.lock())} held for {age_h:.1f} h: {held.download_as_text()}")
        log(f"breaking stale lock ({age_h:.1f} h >= {stale_hours} h): {held.download_as_text()}")
        try:
            held.delete(if_generation_match=held.generation)
        except (precondition_failed(), not_found()):
            pass
    raise Locked(f"{lay.uri(lay.lock())}: lost the race to break a stale lock")


def release_lock(bkt, lay, generation):
    try:
        bkt.blob(lay.lock()).delete(if_generation_match=generation)
    except (precondition_failed(), not_found()):
        log("lock was no longer ours at release; left alone")


# ---------------------------------------------------------------- worker task

def upload_raw(bkt, lay, item, path):
    """No-clobber Coldline upload; 'present' when the object already exists."""
    blob = bkt.blob(lay.raw(item["dataset_version_id"]))
    if blob.exists():
        return "present"
    blob.storage_class = "COLDLINE"
    blob.metadata = {"dataset_id": item["dataset_id"], "source_url": item["url"]}
    try:
        blob.upload_from_filename(path, if_generation_match=0, timeout=3600)
    except precondition_failed():
        return "present"
    return "uploaded"


def vm_name():
    """The GCE instance, so a report shows which tasks shared a VM (the hostname is the container's)."""
    import requests
    try:
        return requests.get("http://metadata.google.internal/computeMetadata/v1/instance/name",
                            headers={"Metadata-Flavor": "Google"}, timeout=2).text
    except Exception:
        return None


def run_task(env=os.environ, client=None, aggregate_run=None):
    """One Batch task: item BATCH_TASK_INDEX of ITEMS_URI. Returns the process exit code."""
    client = client or storage_client()
    lay = Layout(env["BUCKET"], env["ROOT"])
    bkt = client.bucket(lay.bucket)
    items_path = env["ITEMS_URI"].split(f"gs://{lay.bucket}/", 1)[1]
    item = json.loads(bkt.blob(items_path).download_as_bytes())[int(env["BATCH_TASK_INDEX"])]
    vid = item["dataset_version_id"]
    rec = {"dataset_version_id": vid, "dataset_id": item["dataset_id"], "vm": vm_name(),
           "attempt": int(env.get("BATCH_TASK_RETRY_ATTEMPT", 0)), "started_at": now().isoformat()}

    def outcome(status, **kw):
        put_json(bkt, lay.run(env["RUN_ID"], f"outcomes/{vid}.json"),
                 rec | {"status": status, "finished_at": now().isoformat()} | kw)

    if bkt.blob(lay.marker(vid)).exists():
        outcome("skipped_done")
        return 0
    if bkt.blob(lay.rejected(vid)).exists():
        outcome("rejected", reason="_REJECTED.json already present")
        return REJECTED_EXIT
    raw = {}
    original = aggregate.download

    def download_and_keep(url, dest, *a, **kw):
        # aggregate's own downloader, so the marker keeps its ETag-pinned source record;
        # the raw copy goes up before the aggregation spends the time
        src = original(url, dest, *a, **kw)
        t = time.time()
        raw["raw"] = upload_raw(bkt, lay, item, dest)
        raw["raw_upload_s"] = round(time.time() - t, 1)
        return src

    aggregate.download = download_and_keep
    try:
        (aggregate_run or aggregate.run)(argparse.Namespace(
            input=item["url"], out=lay.uri(lay.agg(vid)), dataset_id=item["dataset_id"],
            dataset_version_id=vid, workdir=None, chunk_nnz=aggregate.CHUNK_NNZ,
            max_accum_bytes=aggregate.MAX_ACCUM_BYTES))
    except aggregate.InputRejected as e:
        put_json(bkt, lay.rejected(vid), {
            "dataset_id": item["dataset_id"], "dataset_version_id": vid, "reason": str(e),
            "worker_image": env.get("WORKER_IMAGE_DIGEST"), "run_id": env["RUN_ID"],
            "rejected_at": now().isoformat(timespec="seconds")})
        held = bkt.get_blob(lay.raw(vid))
        if held is not None:
            try:
                held.delete(if_generation_match=held.generation)
                raw["raw"] = "deleted"
            except (precondition_failed(), not_found()):
                pass
        outcome("rejected", reason=str(e), **raw)
        return REJECTED_EXIT
    except Exception as e:
        outcome("error", reason=f"{type(e).__name__}: {e}", **raw)
        raise
    finally:
        aggregate.download = original
    outcome("ok", **raw)
    return 0


# ---------------------------------------------------------------- Batch

class BatchApi:
    def __init__(self, project=PROJECT, region=REGION):
        from google.cloud import batch_v1
        self.batch_v1, self.client = batch_v1, batch_v1.BatchServiceClient()
        self.parent = f"projects/{project}/locations/{region}"

    def submit(self, job_id, spec):
        job = self.batch_v1.Job.from_json(json.dumps(spec))
        return self.client.create_job(parent=self.parent, job=job, job_id=job_id).name

    def cancel(self, name):
        self.client.cancel_job(name=name)

    def status(self, name):
        j = self.client.get_job(name=name)
        counts = {k: dict(v.counts) for k, v in j.status.task_groups.items()}
        events = [(e.event_time.isoformat() if e.event_time else None, e.description)
                  for e in j.status.status_events]
        return {"state": j.status.state.name, "counts": counts, "events": events,
                "run_duration_s": j.status.run_duration.total_seconds() if j.status.run_duration else None}


def resolve_digest(image, session=None):
    """repo:tag -> repo@sha256:..., read from the registry so every task runs one image."""
    if "@sha256:" in image:
        return image
    import google.auth
    import google.auth.transport.requests
    import requests
    repo, _, tag = image.rpartition(":")
    host, _, path = repo.partition("/")
    creds, _ = google.auth.default(scopes=["https://www.googleapis.com/auth/cloud-platform"])
    creds.refresh(google.auth.transport.requests.Request())
    accept = ", ".join(["application/vnd.oci.image.index.v1+json",
                        "application/vnd.docker.distribution.manifest.list.v2+json",
                        "application/vnd.oci.image.manifest.v1+json",
                        "application/vnd.docker.distribution.manifest.v2+json"])
    r = (session or requests).head(f"https://{host}/v2/{path}/manifests/{tag}", timeout=60,
                                   headers={"Authorization": f"Bearer {creds.token}", "Accept": accept})
    r.raise_for_status()
    return f"{repo}@{r.headers['Docker-Content-Digest']}"


def wait(batch, names, poll_s, timeout_h, sleep=time.sleep, clock=time.time):
    deadline = clock() + timeout_h * 3600
    while True:
        st = {n: batch.status(n) for n in names}
        if all(s["state"] in TERMINAL for s in st.values()) or clock() > deadline:
            return st
        log(" ".join(f"{n.rsplit('/', 1)[-1]}={s['state']}" for n, s in st.items()))
        sleep(poll_s)


# ---------------------------------------------------------------- driver

def check_restriction(only, limit):
    """An empty --only or a zero --limit is a mistake, never 'unrestricted'."""
    if only is not None and not [x for x in only.split(",") if x.strip()]:
        raise ValueError("--only names no dataset")
    if limit is not None and limit < 1:
        raise ValueError("--limit must be >= 1")


def restrict_listing(listing, only=None, limit=None):
    if only is not None:
        keep = set(only)
        listing = [d for d in listing if d["dataset_version_id"] in keep or d["dataset_id"] in keep]
    if limit is not None:
        ok = [d for d in listing if plan.classify(d)[0]]
        listing = sorted(ok, key=lambda d: plan.h5ad_asset(d)["filesize"])[:limit]
    return listing


def outcomes(bkt, lay, run_id, items, done, rejected):
    out = {}
    for i in items:
        vid = i["dataset_version_id"]
        rec = bkt.get_blob(lay.run(run_id, f"outcomes/{vid}.json"))
        rec = json.loads(rec.download_as_bytes()) if rec else {}
        status = "ok" if vid in done else "rejected" if vid in rejected else "failed"
        out[vid] = {"status": status, **{k: rec.get(k) for k in ("vm", "attempt", "reason", "raw",
                                                                   "started_at", "finished_at") if k in rec}}
    return out


def mark_superseded(bkt, lay, superseded, done, clock=now):
    """customTime=now on the old raw file, only after its replacement has a marker."""
    # only versions with a _SUCCESS marker are ever retired here: the raw file of a version that
    # failed (rather than was rejected) and was later superseded stays until someone deletes it
    res = []
    for s in superseded:
        old, new = s["dataset_version_id"], s["current_dataset_version_id"]
        if new not in done:
            res.append({**s, "action": "waiting_for_new_version"})
            continue
        blob = bkt.get_blob(lay.raw(old))
        if blob is None:
            action = "no_raw_object"
        elif blob.custom_time is not None:
            action = "already_marked"
        else:
            blob.custom_time = clock()
            blob.patch()
            action = "custom_time_set"
        res.append({**s, "action": action})
    return res


def run_cmd(template, **fields):
    """None when the command or a script it names is absent; else the exit code."""
    argv = shlex.split(template.format(**fields))
    scripts = [t for t in argv if t.endswith((".py", ".sh"))]
    if (not shutil.which(argv[0]) and not Path(argv[0]).exists()) or not all(Path(t).exists() for t in scripts):
        return None
    log("+", " ".join(argv)[:500])
    return subprocess.run(argv).returncode


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def publish_inputs(in_scope, image_digest, publish_cmd, rollup_dir=HERE / "rollup"):
    """Everything that changes a release's content besides the per-dataset outputs themselves."""
    return {"datasets_sha256": hashlib.sha256(",".join(sorted(in_scope)).encode()).hexdigest(),
            "aggregation_format_version": aggregate.FORMAT_VERSION,
            "rollup_sha256": {p.name: sha256_file(p) for p in sorted(Path(rollup_dir).glob("*.tsv"))},
            "publish_sha256": sha256_file(HERE / "publish.py") if (HERE / "publish.py").exists() else None,
            "image": image_digest, "publish_cmd": publish_cmd}


# raising here releases the lock but does not cancel the Batch worker jobs already submitted;
# they run on to their own end, and the next driver run finds their markers
def _sigterm(signum, frame):
    raise SystemExit(128 + signum)


def drive(a, client=None, batch=None, sleep=time.sleep):
    try:
        check_restriction(a.only, a.limit)
    except ValueError as e:
        log(f"driver.py: {e}")
        return 2
    try:
        prev = signal.signal(signal.SIGTERM, _sigterm)
    except ValueError:      # not the main thread: nothing to install
        prev = None
    client = client or storage_client()
    lay = Layout(a.bucket, a.root)
    bkt = client.bucket(lay.bucket)
    run_id = a.run_id or now().strftime("%Y%m%dt%H%M%Sz")
    try:
        gen = acquire_lock(bkt, lay, run_id, a.lock_stale_hours)
    except Locked as e:
        log(f"driver.py: locked: {e}")
        if prev is not None:
            signal.signal(signal.SIGTERM, prev)
        return EXIT_LOCKED
    try:
        return _drive(a, client, bkt, lay, run_id, batch, sleep)
    finally:
        release_lock(bkt, lay, gen)
        if prev is not None:
            signal.signal(signal.SIGTERM, prev)


def _drive(a, client, bkt, lay, run_id, batch, sleep):
    work = Path(tempfile.mkdtemp(prefix="cxg-driver-"))
    full = plan.fetch_listing() if not a.listing else json.loads(Path(a.listing).read_text())
    restricted = a.only is not None or a.limit is not None
    listing = restrict_listing(full, a.only.split(",") if a.only is not None else None, a.limit)
    (work / "listing.json").write_text(json.dumps(listing))
    bkt.blob(lay.run(run_id, "listing.json")).upload_from_filename(str(work / "listing.json"))

    rc = plan.main(["--out", str(work / "plan.json"), "--run-id", run_id, "--bucket", lay.bucket,
                    "--root", lay.root,
                    "--listing", str(work / "listing.json"), "--max-datasets", str(a.max_datasets),
                    "--max-gb", str(a.max_gb)] + (["--approve-large"] if a.approve_large else []))
    if (work / "plan.json").exists():
        bkt.blob(lay.run(run_id, "plan.json")).upload_from_filename(str(work / "plan.json"))
    if rc == plan.BREAKER_EXIT:
        log("driver.py: plan exceeds the circuit breaker; nothing submitted. Review "
            f"{lay.uri(lay.run(run_id, 'plan.json'))} and rerun with --approve-large")
        return EXIT_BREAKER
    if rc != 0:
        return 2
    p = json.loads((work / "plan.json").read_text())
    items = p["items"]

    jobs, cancelled, unconfirmed, digest = {}, [], [], None
    if items:
        batch = batch or BatchApi()
        digest = resolve_digest(a.image)
        log(f"image {digest}")
        for (mc, disk), group in sorted(group_items(items).items()):
            jid = job_id(run_id, mc, disk)
            items_path = lay.run(run_id, f"jobs/{jid}.items.json")
            put_json(bkt, items_path, group)
            spec = build_job(lay.uri(items_path), group, mc, disk, digest, run_id, lay.bucket,
                             lay.root, a.parallelism)
            name = batch.submit(jid, spec)
            jobs[name] = {"machine_class": mc, "disk_gb": disk, "tasks": len(group)}
            log(f"submitted {name}: {len(group)} tasks, {mc}, {disk} GB disk")
        for name, st in wait(batch, list(jobs), a.poll_seconds, a.timeout_hours, sleep).items():
            jobs[name] |= st
            if st["state"] not in TERMINAL:
                log(f"timeout: cancelling {name} ({st['state']})")
                batch.cancel(name)
                jobs[name]["cancelled"] = True
                cancelled.append(name)
        if cancelled:
            # the lock stays held until the cancels are confirmed, or the bound runs out
            for name, st in wait(batch, cancelled, a.poll_seconds, a.cancel_wait_minutes / 60, sleep).items():
                jobs[name]["state_after_cancel"] = st["state"]
                if st["state"] not in TERMINAL:
                    unconfirmed.append(name)
            if unconfirmed:
                log(f"cancel not confirmed within {a.cancel_wait_minutes} min, lock released anyway: "
                    + " ".join(unconfirmed))

    done, rejected_markers = markers(client, lay)
    per = outcomes(bkt, lay, run_id, items, done, rejected_markers)
    failed = sorted(v for v, o in per.items() if o["status"] == "failed")
    rejected = sorted(v for v, o in per.items() if o["status"] == "rejected")
    listed = {d["dataset_version_id"] for d in listing}
    superseded = mark_superseded(bkt, lay, p["superseded"], done)
    withdrawn = None if restricted else sorted(
        (done - listed) - {s["dataset_version_id"] for s in p["superseded"]})
    report = {"run_id": run_id, "plan_totals": p["totals"], "jobs": jobs, "datasets": per,
              "failed": failed, "rejected": rejected, "previously_rejected": p["rejected"],
              "cancelled_jobs": cancelled, "cancel_unconfirmed": unconfirmed, "superseded": superseded,
              "withdrawn_not_acted_on": withdrawn, "publish": None}

    def finish(code, **kw):
        report["publish"] = kw
        put_json(bkt, lay.run(run_id, "report.json"), report)
        log(json.dumps({k: report[k] for k in ("run_id", "failed", "rejected", "publish")}))
        log(f"datasets: {sum(o['status'] == 'ok' for o in per.values())} ok, {len(rejected)} rejected, "
            f"{len(failed)} failed; report {lay.uri(lay.run(run_id, 'report.json'))}")
        return code

    if failed and not a.allow_partial:
        return finish(EXIT_FAILED, skipped="failed datasets; rerun resumes from markers, "
                                           "or pass --allow-partial")
    if a.no_publish or restricted:
        return finish(0, skipped="--no-publish" if a.no_publish else "restricted listing (--only/--limit)")

    # the plan's listing, not every marker under agg/: a withdrawn dataset must drop out
    in_scope = sorted(d["dataset_version_id"] for d in listing if plan.classify(d)[0]
                      and d["dataset_version_id"] in done)
    inputs = publish_inputs(in_scope, digest or resolve_digest(a.image), a.publish_cmd)
    key = hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()
    cur = bkt.get_blob(lay.current())
    cur = json.loads(cur.download_as_bytes()) if cur is not None else {}
    if cur.get("publish_key") == key:
        return finish(0, skipped="CURRENT was built from exactly these inputs")
    if cur.get("datasets") and len(in_scope) < SHRINK_FLOOR * cur["datasets"] and not a.allow_shrink:
        return finish(EXIT_SHRINK, error=f"release would hold {len(in_scope)} datasets, CURRENT holds "
                                         f"{cur['datasets']} (< {SHRINK_FLOOR:.0%}); nothing published "
                                         "(--allow-shrink)")
    t = now()
    release = a.release or f"{t.date().isoformat()}-{t.strftime('%H%M%S')}"
    releases = lay.uri(f"{lay.root}/releases")
    if any(True for _ in client.list_blobs(lay.bucket, prefix=f"{lay.root}/releases/{release}/")):
        return finish(EXIT_PUBLISH, error=f"{releases}/{release}/ already exists; a release is never rewritten")
    fields = {"agg": lay.uri(lay.agg()).rstrip("/"), "releases": releases, "release": release,
              "release_uri": f"{releases}/{release}", "datasets": ",".join(in_scope),
              "listing": str(work / "listing.json"), "run_id": run_id}
    rc = run_cmd(a.publish_cmd, **fields)
    if rc is None:
        return finish(EXIT_PUBLISH, error=f"publish command not found: {a.publish_cmd}")
    if rc != 0:
        return finish(EXIT_PUBLISH, error=f"publish exited {rc}")
    rc = run_cmd(a.load_cmd, **fields)
    if rc is None:
        log(f"load command not found ({a.load_cmd}); BigQuery not loaded")
        if not a.allow_missing_load:
            return finish(EXIT_PUBLISH, release=release, error="load command absent; CURRENT not "
                                                               "written (--allow-missing-load)")
    elif rc != 0:
        return finish(EXIT_PUBLISH, release=release, error=f"load exited {rc}; CURRENT not written")
    put_json(bkt, lay.current(), {"release": release, "uri": f"{releases}/{release}/",
                                  "run_id": run_id, "written_at": now().isoformat(timespec="seconds"),
                                  "bigquery_loaded": rc is not None, "datasets": len(in_scope),
                                  "publish_key": key, "publish_inputs": inputs})
    return finish(0, release=release, current=lay.uri(lay.current()))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("task", help="Batch worker entrypoint (env: RUN_ID ITEMS_URI BUCKET ROOT BATCH_TASK_INDEX)")
    r = sub.add_parser("run")
    r.add_argument("--bucket", default=BUCKET)
    r.add_argument("--root", default=ROOT, help="bucket prefix; cellxgene-staging for rehearsals")
    r.add_argument("--image", default=IMAGE, help="worker image, tag or @sha256 digest")
    r.add_argument("--run-id")
    r.add_argument("--listing", help="saved curation API listing instead of calling the API")
    r.add_argument("--only", help="comma-separated dataset_version_ids or dataset_ids")
    r.add_argument("--limit", type=int, help="the N smallest in-scope datasets of the listing")
    r.add_argument("--approve-large", action="store_true")
    r.add_argument("--max-datasets", type=int, default=200)
    r.add_argument("--max-gb", type=float, default=300.0)
    r.add_argument("--parallelism", type=int, default=100, help="max VMs per Batch job")
    r.add_argument("--poll-seconds", type=int, default=60)
    r.add_argument("--timeout-hours", type=float, default=40)
    r.add_argument("--cancel-wait-minutes", type=float, default=15,
                   help="after a timeout, how long to hold the lock waiting for cancels to land")
    r.add_argument("--lock-stale-hours", type=float, default=LOCK_STALE_HOURS)
    r.add_argument("--allow-partial", action="store_true")
    r.add_argument("--no-publish", action="store_true")
    r.add_argument("--release", help="release name (default: <UTC date>-<HHMMSS>)")
    r.add_argument("--publish-cmd", default=PUBLISH_CMD)
    r.add_argument("--load-cmd", default=LOAD_CMD)
    r.add_argument("--allow-missing-load", action="store_true")
    r.add_argument("--allow-shrink", action="store_true",
                   help=f"flip CURRENT even when the release holds < {SHRINK_FLOOR:.0%} of its datasets")
    a = p.parse_args(argv)
    if a.cmd == "task":
        return run_task()
    try:
        return drive(a)
    except Exception as e:
        log(f"driver.py: {type(e).__name__}: {e}")
        raise


if __name__ == "__main__":
    sys.exit(main())
