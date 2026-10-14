"""driver.py against in-memory GCS and Batch fakes; no network.

Run (needs google-cloud-storage, so inside the pipeline image):
    docker run --rm -v "$PWD:/app" cellxgene-worker:local python -m pytest -p no:cacheprovider tests/cellxgene/test_driver.py
"""
import argparse
import datetime as dt
import json
import os
import re
import signal
import sys

import pytest
from google.api_core.exceptions import NotFound, PreconditionFailed

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts", "cellxgene"))
import aggregate  # noqa: E402
import driver  # noqa: E402
import plan  # noqa: E402

with open(os.path.join(HERE, "listing_fixture.json")) as fh:
    LISTING = json.load(fh)
IN_SCOPE = [d for d in LISTING if plan.classify(d)[0]]
T0 = dt.datetime(2026, 10, 13, 12, tzinfo=dt.timezone.utc)


# ---------------------------------------------------------------- fakes

class Blob:
    def __init__(self, store, name):
        self.store, self.name = store, name
        self.storage_class = self.metadata = self.custom_time = None
        self.generation = self.time_created = None

    def _put(self, data, if_generation_match=None):
        cur = self.store.objs.get(self.name)
        if if_generation_match is not None and (cur["generation"] if cur else 0) != if_generation_match:
            raise PreconditionFailed("precondition")
        self.store.gen += 1
        self.generation = self.store.gen
        self.store.objs[self.name] = {"data": data, "generation": self.generation,
                                      "storage_class": self.storage_class or "STANDARD",
                                      "metadata": self.metadata, "custom_time": None,
                                      "time_created": self.store.clock()}
        self.store.log.append(("put", self.name))

    def upload_from_string(self, s, content_type=None, if_generation_match=None):
        self._put(s.encode() if isinstance(s, str) else s, if_generation_match)

    def upload_from_filename(self, path, if_generation_match=None, timeout=None):
        with open(path, "rb") as fh:
            self._put(fh.read(), if_generation_match)

    def exists(self):
        return self.name in self.store.objs

    def download_as_bytes(self):
        if self.name not in self.store.objs:
            raise NotFound(self.name)
        return self.store.objs[self.name]["data"]

    def download_as_text(self):
        return self.download_as_bytes().decode()

    def delete(self, if_generation_match=None):
        cur = self.store.objs.get(self.name)
        if cur is None:
            raise NotFound(self.name)
        if if_generation_match is not None and cur["generation"] != if_generation_match:
            raise PreconditionFailed("precondition")
        del self.store.objs[self.name]

    def patch(self):
        self.store.objs[self.name]["custom_time"] = self.custom_time


class Bucket:
    def __init__(self, store):
        self.store = store

    def blob(self, name):
        return Blob(self.store, name)

    def get_blob(self, name):
        o = self.store.objs.get(name)
        if o is None:
            return None
        b = Blob(self.store, name)
        b.generation, b.time_created, b.custom_time = o["generation"], o["time_created"], o["custom_time"]
        b.storage_class = o["storage_class"]
        return b


class Storage:
    def __init__(self):
        self.objs, self.gen, self.log = {}, 0, []
        self.clock = lambda: T0

    def bucket(self, name):
        return Bucket(self)

    def list_blobs(self, bucket, prefix=""):
        return [argparse.Namespace(name=n) for n in list(self.objs) if n.startswith(prefix)]

    def json(self, name):
        return json.loads(self.objs[name]["data"])


class FakeBatch:
    """Runs every task of a submitted job synchronously through `outcome(item)`."""
    def __init__(self, store, outcome, state="SUCCEEDED"):
        self.store, self.outcome, self.jobs, self.state, self.cancelled = store, outcome, {}, state, []

    def submit(self, job_id, spec):
        env = spec["taskGroups"][0]["taskSpec"]["environment"]["variables"]
        items = json.loads(self.store.objs[env["ITEMS_URI"].split("/", 3)[3]]["data"])
        lay = driver.Layout(env["BUCKET"], env["ROOT"])
        for i in items:
            res = self.outcome(i)
            vid = i["dataset_version_id"]
            if res == "ok":
                self.store.objs[lay.marker(vid)] = {"data": b"{}", "generation": 1, "custom_time": None,
                                                    "storage_class": "STANDARD", "time_created": T0}
            if res == "rejected":
                Bucket(self.store).blob(lay.rejected(vid)).upload_from_string("{}")
            if res in ("ok", "rejected"):
                Bucket(self.store).blob(lay.run(env["RUN_ID"], f"outcomes/{vid}.json")).upload_from_string(
                    json.dumps({"status": res}))
        self.jobs[job_id] = spec
        return f"projects/p/locations/r/jobs/{job_id}"

    def status(self, name):
        state = self.state() if callable(self.state) else self.state
        return {"state": state, "counts": {}, "events": [], "run_duration_s": 1}

    def cancel(self, name):
        self.cancelled.append(name)


def args(**kw):
    base = dict(bucket="b", root="cellxgene-staging", image="r/cellxgene@sha256:abc", run_id="r1",
                listing=None, only=None, limit=None, approve_large=False, max_datasets=200,
                max_gb=300.0, parallelism=100, poll_seconds=0, timeout_hours=1,
                lock_stale_hours=72, allow_partial=False, no_publish=False, release="2026-10-13",
                publish_cmd="true", load_cmd="true", allow_missing_load=False, allow_shrink=False,
                cancel_wait_minutes=0)
    return argparse.Namespace(**(base | kw))


@pytest.fixture
def listing(tmp_path):
    p = tmp_path / "listing.json"
    p.write_text(json.dumps(LISTING))
    return str(p)


@pytest.fixture
def store(monkeypatch):
    s = Storage()
    monkeypatch.setattr(plan, "_client", lambda: s)
    return s


def run_driver(store, outcome=lambda i: "ok", state="SUCCEEDED", **kw):
    batch = FakeBatch(store, outcome, state)
    return driver.drive(args(**kw), client=store, batch=batch, sleep=lambda s: None), batch


# ---------------------------------------------------------------- sizing and job spec

@pytest.mark.parametrize("gb,disk", [(0.1, 30), (5, 30), (5.1, 50), (15, 50), (16, 100),
                                     (40, 100), (45, 150), (53, 150)])
def test_disk_is_twice_the_file_plus_20_in_tiers(gb, disk):
    assert driver.disk_gb(gb * 1e9) == disk


def test_items_group_by_machine_class_and_disk_tier():
    items = [{"machine_class": plan.STANDARD, "filesize": 1e9}, {"machine_class": plan.STANDARD, "filesize": 2e9},
             {"machine_class": plan.STANDARD, "filesize": 45e9}, {"machine_class": plan.HIGHMEM, "filesize": 1e9}]
    g = driver.group_items(items)
    assert {k: len(v) for k, v in g.items()} == {(plan.STANDARD, 30): 2, (plan.STANDARD, 150): 1,
                                                 (plan.HIGHMEM, 30): 1}
    assert [i["filesize"] for i in g[(plan.STANDARD, 30)]] == [2e9, 1e9]


def test_job_spec_spot_external_ip_worker_sa_retry_policy():
    items = [{"filesize": 2e9}, {"filesize": 1e9}, {"filesize": 3e9}]
    j = driver.build_job("gs://b/x.json", items, plan.HIGHMEM, 30, "r/c@sha256:d", "R1", "b", "root", 2)
    tg = j["taskGroups"][0]
    assert (tg["taskCount"], tg["parallelism"]) == (3, 2)
    spec = tg["taskSpec"]
    assert spec["runnables"][0]["container"]["imageUri"] == "r/c@sha256:d"
    assert spec["environment"]["variables"]["WORKER_IMAGE_DIGEST"] == "r/c@sha256:d"
    assert spec["maxRetryCount"] >= 1
    assert spec["maxRunDuration"] == f"{driver.task_timeout_s(3e9)}s"
    # Batch accepts ONE lifecycle policy; spot preemption (50001) matches none and is retried
    # by maxRetryCount, while an input rejection is failed at once
    assert spec["lifecyclePolicies"] == [{"action": "FAIL_TASK", "actionCondition": {"exitCodes": [driver.REJECTED_EXIT]}}]
    ap = j["allocationPolicy"]
    inst = ap["instances"][0]["policy"]
    assert (inst["machineType"], inst["provisioningModel"], inst["bootDisk"]["sizeGb"]) == (plan.HIGHMEM, "SPOT", "30")
    nic = ap["network"]["networkInterfaces"][0]
    assert (nic["network"], nic["subnetwork"], nic["noExternalIpAddress"]) == (driver.NETWORK, driver.SUBNET, False)
    assert ap["serviceAccount"] == {"email": driver.WORKER_SA,
                                    "scopes": ["https://www.googleapis.com/auth/cloud-platform"]}
    assert j["labels"]["run-id"] == "r1"


def test_job_id_is_a_valid_batch_id():
    import re
    jid = driver.job_id("20261013t192000z", plan.HIGHMEM, 150)
    assert re.fullmatch(r"[a-z]([a-z0-9-]{0,61}[a-z0-9])?", jid), jid


def test_driver_job_spec_is_the_driver_sa_without_a_job_id():
    j = json.loads((driver.HERE / "driver-batch-job.json").read_text())
    assert j["allocationPolicy"]["serviceAccount"]["email"].startswith("expression-batch-driver@")
    assert "name" not in j and "uid" not in j
    cmd = j["taskGroups"][0]["taskSpec"]["runnables"][0]["container"]["commands"]
    assert cmd[:2] == ["/app/scripts/cellxgene/driver.py", "run"]
    stale_s = driver.LOCK_STALE_HOURS * 3600
    assert int(j["taskGroups"][0]["taskSpec"]["maxRunDuration"].rstrip("s")) < stale_s


def test_resolve_digest_reads_the_registry_header(monkeypatch):
    import google.auth
    creds = argparse.Namespace(token="t", refresh=lambda req: None)
    monkeypatch.setattr(google.auth, "default", lambda scopes: (creds, "p"))
    seen = {}

    class S:
        def head(self, url, timeout, headers):
            seen["url"] = url
            return argparse.Namespace(raise_for_status=lambda: None,
                                      headers={"Docker-Content-Digest": "sha256:feed"})
    out = driver.resolve_digest("us-docker.pkg.dev/p/repo/cellxgene:dev-1", session=S())
    assert out == "us-docker.pkg.dev/p/repo/cellxgene@sha256:feed"
    assert seen["url"] == "https://us-docker.pkg.dev/v2/p/repo/cellxgene/manifests/dev-1"
    assert driver.resolve_digest("x/y@sha256:1") == "x/y@sha256:1"


# ---------------------------------------------------------------- lock

def test_lock_blocks_a_second_driver_and_releases_only_its_own():
    s = Storage()
    b, lay = s.bucket("b"), driver.Layout("b", "root")
    gen = driver.acquire_lock(b, lay, "r1")
    with pytest.raises(driver.Locked):
        driver.acquire_lock(b, lay, "r2")
    driver.release_lock(b, lay, gen + 99)          # not ours: left in place
    assert lay.lock() in s.objs
    driver.release_lock(b, lay, gen)
    assert lay.lock() not in s.objs


def test_stale_lock_is_taken_over():
    s = Storage()
    b, lay = s.bucket("b"), driver.Layout("b", "root")
    driver.acquire_lock(b, lay, "old")
    later = lambda: T0 + dt.timedelta(hours=73)  # noqa: E731
    s.clock = later
    driver.acquire_lock(b, lay, "new", clock=later)
    assert json.loads(s.objs[lay.lock()]["data"])["run_id"] == "new"


def test_locked_driver_exits_75_and_submits_nothing(store, listing):
    driver.acquire_lock(store.bucket("b"), driver.Layout("b", "cellxgene-staging"), "other")
    rc, batch = run_driver(store, listing=listing)
    assert rc == driver.EXIT_LOCKED and not batch.jobs


# ---------------------------------------------------------------- worker task

@pytest.fixture(autouse=True)
def no_metadata(monkeypatch):
    monkeypatch.setattr(driver, "vm_name", lambda: "vm-1")


def task_env(store, item, run_id="r1"):
    lay = driver.Layout("b", "root")
    store.bucket("b").blob(lay.run(run_id, "jobs/j.items.json")).upload_from_string(json.dumps([item]))
    return {"RUN_ID": run_id, "ITEMS_URI": f"gs://b/{lay.run(run_id, 'jobs/j.items.json')}",
            "BUCKET": "b", "ROOT": "root", "BATCH_TASK_INDEX": "0"}, lay


ITEM = {"dataset_id": "d1", "dataset_version_id": "v1", "url": "https://x/v1.h5ad", "filesize": 10}


def fake_aggregate(tmp_path, reject=False):
    def run(ns):
        src = tmp_path / "in.h5ad"
        aggregate.download(ns.input, str(src))
        if reject:
            raise aggregate.InputRejected("float matrix")
        return {}
    return run


@pytest.fixture
def fake_download(monkeypatch):
    def dl(url, dest, *a, **kw):
        with open(dest, "wb") as fh:
            fh.write(b"h5ad-bytes")
        return {"url": url}
    monkeypatch.setattr(aggregate, "download", dl)


def test_task_uploads_raw_coldline_no_clobber(store, tmp_path, fake_download):
    env, lay = task_env(store, ITEM)
    assert driver.run_task(env, store, fake_aggregate(tmp_path)) == 0
    raw = store.objs[lay.raw("v1")]
    assert raw["storage_class"] == "COLDLINE" and raw["data"] == b"h5ad-bytes"
    gen = raw["generation"]
    out = store.json(lay.run("r1", "outcomes/v1.json"))
    assert (out["status"], out["raw"]) == ("ok", "uploaded")
    assert driver.run_task(env, store, fake_aggregate(tmp_path)) == 0    # spot retry
    assert store.objs[lay.raw("v1")]["generation"] == gen
    assert store.json(lay.run("r1", "outcomes/v1.json"))["raw"] == "present"
    assert aggregate.download.__name__ == "dl"     # the patch is undone


def test_raw_upload_race_is_present_not_an_error(store, tmp_path):
    lay = driver.Layout("b", "root")
    b = store.bucket("b")
    p = tmp_path / "f"
    p.write_bytes(b"x")
    blob = b.blob(lay.raw("v1"))
    blob.exists = lambda: False                     # a racing upload lands between check and write
    store.objs[lay.raw("v1")] = {"generation": 5}
    b.blob = lambda name: blob
    assert driver.upload_raw(b, lay, ITEM, str(p)) == "present"


def test_task_with_marker_is_skipped(store, tmp_path, fake_download):
    env, lay = task_env(store, ITEM)
    store.bucket("b").blob(lay.marker("v1")).upload_from_string("{}")
    assert driver.run_task(env, store, lambda ns: pytest.fail("aggregated a done dataset")) == 0
    assert store.json(lay.run("r1", "outcomes/v1.json"))["status"] == "skipped_done"


def test_rejected_input_exits_3_and_is_recorded(store, tmp_path, fake_download):
    env, lay = task_env(store, ITEM)
    assert driver.run_task(env, store, fake_aggregate(tmp_path, reject=True)) == driver.REJECTED_EXIT
    out = store.json(lay.run("r1", "outcomes/v1.json"))
    assert out["status"] == "rejected" and "float matrix" in out["reason"]
    marker = store.json(lay.rejected("v1"))
    assert marker["reason"] == "float matrix" and marker["dataset_version_id"] == "v1"
    assert lay.raw("v1") not in store.objs and out["raw"] == "deleted"
    # a rejected version is never re-aggregated, even if a task for it runs again
    assert driver.run_task(env, store, lambda ns: pytest.fail("re-aggregated")) == driver.REJECTED_EXIT


# ---------------------------------------------------------------- driver run

def test_breaker_stops_before_submitting(store, listing):
    rc, batch = run_driver(store, listing=listing, max_datasets=1)
    assert rc == driver.EXIT_BREAKER and not batch.jobs
    assert "cellxgene-staging/runs/r1/plan.json" in store.objs
    assert "cellxgene-staging/driver.lock" not in store.objs


def test_approve_large_passes_through(store, listing):
    rc, batch = run_driver(store, listing=listing, max_datasets=1, approve_large=True)
    assert rc == 0 and batch.jobs


def test_full_run_submits_publishes_and_writes_current_last(store, listing):
    rc, batch = run_driver(store, listing=listing)
    assert rc == 0
    n = sum(spec["taskGroups"][0]["taskCount"] for spec in batch.jobs.values())
    assert n == len(IN_SCOPE)
    puts = [n for op, n in store.log if op == "put"]
    assert puts[-1] == "cellxgene-staging/runs/r1/report.json"
    assert puts[-2] == "cellxgene-staging/releases/CURRENT"
    cur = store.json("cellxgene-staging/releases/CURRENT")
    assert cur["release"] == "2026-10-13" and cur["datasets"] == len(IN_SCOPE)


def test_rerun_is_a_no_op(store, listing):
    run_driver(store, listing=listing)
    gen = store.objs["cellxgene-staging/releases/CURRENT"]["generation"]
    rc, batch = run_driver(store, listing=listing, run_id="r2")
    assert rc == 0 and not batch.jobs
    assert store.objs["cellxgene-staging/releases/CURRENT"]["generation"] == gen
    assert "exactly these inputs" in store.json("cellxgene-staging/runs/r2/report.json")["publish"]["skipped"]


def test_failed_dataset_blocks_publish_unless_allow_partial(store, listing):
    bad = IN_SCOPE[0]["dataset_version_id"]
    outcome = lambda i: "failed" if i["dataset_version_id"] == bad else "ok"  # noqa: E731
    rc, _ = run_driver(store, outcome, listing=listing)
    assert rc == driver.EXIT_FAILED
    assert "cellxgene-staging/releases/CURRENT" not in store.objs
    rep = store.json("cellxgene-staging/runs/r1/report.json")
    assert rep["failed"] == [bad] and rep["datasets"][bad]["status"] == "failed"
    rc, batch = run_driver(store, outcome, listing=listing, run_id="r2", allow_partial=True)
    assert rc == 0 and sum(s["taskGroups"][0]["taskCount"] for s in batch.jobs.values()) == 1
    assert "cellxgene-staging/releases/CURRENT" in store.objs


def test_rejected_dataset_does_not_block_publish(store, listing):
    bad = IN_SCOPE[0]["dataset_version_id"]
    rc, _ = run_driver(store, lambda i: "rejected" if i["dataset_version_id"] == bad else "ok",
                       listing=listing)
    assert rc == 0
    assert store.json("cellxgene-staging/runs/r1/report.json")["rejected"] == [bad]
    rc, batch = run_driver(store, listing=listing, run_id="r2")
    assert rc == 0 and not batch.jobs
    rep = store.json("cellxgene-staging/runs/r2/report.json")
    assert [r["dataset_version_id"] for r in rep["previously_rejected"]] == [bad]
    assert rep["plan_totals"]["rejected"] == 1


def test_restricted_run_never_publishes(store, listing):
    rc, batch = run_driver(store, listing=listing, limit=2)
    assert rc == 0 and sum(s["taskGroups"][0]["taskCount"] for s in batch.jobs.values()) == 2
    assert "cellxgene-staging/releases/CURRENT" not in store.objs
    vid = IN_SCOPE[0]["dataset_version_id"]
    rc, batch = run_driver(store, listing=listing, only=vid, run_id="r2")
    submitted = [i["dataset_version_id"] for j in batch.jobs
                 for i in json.loads(store.objs[f"cellxgene-staging/runs/r2/jobs/{j}.items.json"]["data"])]
    assert set(submitted) <= {vid}


def test_absent_publish_fails_absent_load_skips_only_when_allowed(store, listing):
    rc, _ = run_driver(store, listing=listing, publish_cmd="/nonexistent/publish.py")
    assert rc == driver.EXIT_PUBLISH
    assert "cellxgene-staging/releases/CURRENT" not in store.objs
    rc, _ = run_driver(store, listing=listing, run_id="r2", load_cmd="/nonexistent/load_cellxgene.sh {release_uri}")
    assert rc == driver.EXIT_PUBLISH and "cellxgene-staging/releases/CURRENT" not in store.objs
    rc, _ = run_driver(store, listing=listing, run_id="r3", allow_missing_load=True,
                       load_cmd="/nonexistent/load_cellxgene.sh {release_uri}")
    assert rc == 0 and store.json("cellxgene-staging/releases/CURRENT")["bigquery_loaded"] is False


def test_superseded_raw_marked_only_after_new_version_succeeds(store, listing, tmp_path):
    new = IN_SCOPE[0]
    lay = driver.Layout("b", "cellxgene-staging")
    b = store.bucket("b")
    b.blob(lay.marker("old")).upload_from_string(json.dumps({"dataset_id": new["dataset_id"]}))
    b.blob(lay.raw("old")).upload_from_string("raw")
    b.blob(lay.raw(new["dataset_version_id"])).upload_from_string("raw")
    rc, _ = run_driver(store, lambda i: "failed" if i["dataset_version_id"] == new["dataset_version_id"]
                       else "ok", listing=listing)
    assert rc == driver.EXIT_FAILED
    assert store.objs[lay.raw("old")]["custom_time"] is None
    rep = store.json(lay.run("r1", "report.json"))
    assert rep["superseded"][0]["action"] == "waiting_for_new_version"
    rc, _ = run_driver(store, listing=listing, run_id="r2")
    marked = store.objs[lay.raw("old")]["custom_time"]
    assert abs((marked - dt.datetime.now(dt.timezone.utc)).total_seconds()) < 60
    assert store.objs[lay.raw(new["dataset_version_id"])]["custom_time"] is None
    rc, _ = run_driver(store, listing=listing, run_id="r3")
    assert store.json(lay.run("r3", "report.json"))["superseded"][0]["action"] == "already_marked"


def test_withdrawn_dataset_is_reported_not_touched(store, listing):
    lay = driver.Layout("b", "cellxgene-staging")
    b = store.bucket("b")
    b.blob(lay.marker("gone")).upload_from_string(json.dumps({"dataset_id": "no-longer-listed"}))
    b.blob(lay.raw("gone")).upload_from_string("raw")
    run_driver(store, listing=listing)
    assert store.json(lay.run("r1", "report.json"))["withdrawn_not_acted_on"] == ["gone"]
    assert store.objs[lay.raw("gone")]["custom_time"] is None
    assert lay.marker("gone") in store.objs


def test_publish_gets_the_listing_datasets_not_every_marker(store, listing, tmp_path):
    lay = driver.Layout("b", "cellxgene-staging")
    store.bucket("b").blob(lay.marker("gone")).upload_from_string(json.dumps({"dataset_id": "x"}))
    out = tmp_path / "datasets.txt"
    rc, _ = run_driver(store, listing=listing, publish_cmd=f"sh -c 'echo $0 > {out}' {{datasets}}")
    assert rc == 0
    assert sorted(out.read_text().strip().split(",")) == sorted(d["dataset_version_id"] for d in IN_SCOPE)


# ---------------------------------------------------------------- fix round

@pytest.mark.parametrize("kw", [{"only": ""}, {"only": " , "}, {"limit": 0}])
def test_empty_restriction_is_an_error(store, listing, kw):
    rc, batch = run_driver(store, listing=listing, **kw)
    assert rc == 2 and not batch.jobs and not store.objs


def test_sigterm_releases_the_lock(store, listing):
    def term():
        os.kill(os.getpid(), signal.SIGTERM)
        return "RUNNING"
    before = signal.getsignal(signal.SIGTERM)
    with pytest.raises(SystemExit) as e:
        run_driver(store, lambda i: "pending", state=term, listing=listing)
    assert e.value.code == 128 + signal.SIGTERM
    assert "cellxgene-staging/driver.lock" not in store.objs
    assert signal.getsignal(signal.SIGTERM) == before


def test_timeout_cancels_running_jobs_before_releasing_the_lock(store, listing):
    batch = FakeBatch(store, lambda i: "pending", "RUNNING")
    lock_held = []
    cancel = batch.cancel

    def cancel_spy(name):
        lock_held.append("cellxgene-staging/driver.lock" in store.objs)
        cancel(name)
    batch.cancel = cancel_spy
    rc = driver.drive(args(listing=listing, timeout_hours=0), client=store, batch=batch, sleep=lambda s: None)
    assert rc == driver.EXIT_FAILED
    assert sorted(batch.cancelled) == sorted(f"projects/p/locations/r/jobs/{j}" for j in batch.jobs)
    assert lock_held and all(lock_held)
    assert "cellxgene-staging/driver.lock" not in store.objs
    rep = store.json("cellxgene-staging/runs/r1/report.json")
    assert sorted(rep["cancelled_jobs"]) == sorted(batch.cancelled)
    # jobs that never reach a terminal state within the bound are named, and the lock released
    assert sorted(rep["cancel_unconfirmed"]) == sorted(batch.cancelled)


def test_lock_held_until_cancels_are_confirmed(store, listing):
    polls_after_cancel = []

    def state():
        if not batch.cancelled:
            return "RUNNING"
        polls_after_cancel.append("cellxgene-staging/driver.lock" in store.objs)
        return "RUNNING" if len(polls_after_cancel) < 3 else "CANCELLED"
    batch = FakeBatch(store, lambda i: "pending", state)
    rc = driver.drive(args(listing=listing, timeout_hours=0, cancel_wait_minutes=60), client=store,
                      batch=batch, sleep=lambda s: None)
    assert rc == driver.EXIT_FAILED
    assert len(polls_after_cancel) >= 3 and all(polls_after_cancel)
    assert "cellxgene-staging/driver.lock" not in store.objs
    rep = store.json("cellxgene-staging/runs/r1/report.json")
    assert rep["cancelled_jobs"] and rep["cancel_unconfirmed"] == []
    assert {j["state_after_cancel"] for j in rep["jobs"].values()} == {"CANCELLED"}


def test_release_name_is_unique_and_never_rewritten(store, listing):
    rc, _ = run_driver(store, listing=listing, release=None)
    cur = store.json("cellxgene-staging/releases/CURRENT")
    assert rc == 0 and re.fullmatch(r"\d{4}-\d\d-\d\d-\d{6}", cur["release"])
    store.objs.pop("cellxgene-staging/releases/CURRENT")
    store.bucket("b").blob("cellxgene-staging/releases/taken/manifest.json").upload_from_string("{}")
    rc, _ = run_driver(store, listing=listing, run_id="r2", release="taken")
    assert rc == driver.EXIT_PUBLISH and "cellxgene-staging/releases/CURRENT" not in store.objs


def test_shrink_guard_stops_before_publish(store, listing, tmp_path):
    store.bucket("b").blob("cellxgene-staging/releases/CURRENT").upload_from_string(
        json.dumps({"release": "old", "datasets": 10 * len(IN_SCOPE)}))
    ran = tmp_path / "ran"
    rc, _ = run_driver(store, listing=listing, publish_cmd=f"touch {ran}")
    assert rc == driver.EXIT_SHRINK and not ran.exists()
    assert store.json("cellxgene-staging/releases/CURRENT")["release"] == "old"
    rc, _ = run_driver(store, listing=listing, run_id="r2", publish_cmd=f"touch {ran}", allow_shrink=True)
    assert rc == 0 and ran.exists()
    assert store.json("cellxgene-staging/releases/CURRENT")["release"] == "2026-10-13"


def test_publish_key_covers_more_than_the_dataset_list(store, listing, tmp_path):
    run_driver(store, listing=listing)
    rc, batch = run_driver(store, listing=listing, run_id="r2", release="r2", image="r/cellxgene@sha256:new")
    assert rc == 0 and not batch.jobs
    assert store.json("cellxgene-staging/releases/CURRENT")["release"] == "r2"
    (tmp_path / "cell_type_rollup.tsv").write_text("a\tb\n")
    k1 = driver.publish_inputs(["v"], "img", "cmd", tmp_path)
    (tmp_path / "cell_type_rollup.tsv").write_text("a\tc\n")
    k2 = driver.publish_inputs(["v"], "img", "cmd", tmp_path)
    assert k1["rollup_sha256"] != k2["rollup_sha256"] and k1["datasets_sha256"] == k2["datasets_sha256"]
    assert driver.publish_inputs(["v"], "img", "cmd --min-cells 10", tmp_path) != k2
    cur = store.json("cellxgene-staging/releases/CURRENT")
    assert set(cur["publish_inputs"]) >= {"datasets_sha256", "rollup_sha256", "publish_sha256", "image",
                                          "aggregation_format_version", "publish_cmd"}
