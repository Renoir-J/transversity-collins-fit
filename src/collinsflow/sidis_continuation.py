"""Immutable same-attempt transport for frozen complete Julia jobs.

This module has only standard-library imports. The launcher imports the frozen
collinsflow package in an isolated process after verifying its exact Git blobs.
Python never evaluates a physics model or a derivative.
"""
from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import uuid

FROZEN_HEAD = "bec12db0d65ed67a74de80bca949e21b80a99757"
CONTRACT_SHA = "47059f3236115d0c171bb7b487929b13ca579d876e61839cb33632f0a62a2417"
POLICY_COMMIT = "e9f86df76ab5f1bedbea52290a0aad5f70be7a2e"
SCOPE = "same_attempt_reconstructed_full132_continuation_not_fitting"
PREFIX = 13
FLAGS = ("production_ready", "physical_extraction", "historical_fit_ready", "experimental_data_admitted",
         "fit_performed", "accepted_prediction_diagnostic", "accepted_fitting_gate", "G2PERF_pass",
         "interval_coverage_verified", "all132_coverage_verified")
POINT_PHASES = ("native123", "native4", "tailref123", "reference4")
GLOBAL_PHASES = ("native-select3", "native-selection", "source-math", "derivative-select3", "diagnostic")


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def atomic_json(path, value):
    path = Path(path)
    fd, temp = tempfile.mkstemp(prefix=path.name+".", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as stream:
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp, path)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)


def verify_git(root, expected=None):
    """Compare raw working bytes with every committed blob, in two Git calls."""
    root = Path(root).resolve()
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    if expected is not None and head != expected:
        raise ValueError("wrong committed source head")
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=root):
        raise ValueError("dirty source checkout")
    entries = subprocess.check_output(["git", "ls-tree", "-rz", "HEAD"], cwd=root).split(b"\0")
    hashes = {}
    for entry in filter(None, entries):
        meta, name = entry.split(b"\t", 1)
        _, kind, oid = meta.split()
        name = name.decode()
        if kind != b"blob" or "\n" in name or "\r" in name or '"' in name:
            raise ValueError("non-file or ambiguous Git source path")
        path = root/name
        if path.is_symlink() or not path.is_file():
            raise ValueError("missing or aliased committed source")
        data = path.read_bytes()
        actual = hashlib.sha1(f"blob {len(data)}\0".encode()+data).hexdigest()
        if actual != oid.decode():
            raise ValueError("working bytes differ from committed Git blob: "+name)
        hashes[name] = hashlib.sha256(data).hexdigest()
    return head, hashes


def contract(adapter):
    path = Path(adapter)/"config/sidis-continuation-contract.json"
    if file_hash(path) != CONTRACT_SHA:
        raise ValueError("changed prospective continuation contract")
    result = read_json(path)
    if result["scope"] != SCOPE or any(result.get(f) is not False for f in FLAGS):
        raise ValueError("self-admitting continuation contract")
    return result


def point_folder(manifest, phase, index=0):
    root = Path(manifest["run"])
    if phase in POINT_PHASES:
        if isinstance(index, bool) or not isinstance(index, int) or not 1 <= index <= 132:
            raise ValueError("exact source point index1..132 required")
        if phase == "native123" and index <= PREFIX:
            root = Path(manifest["original"])/"full"
        return root/phase/f"p{index:04d}"
    if phase not in GLOBAL_PHASES or type(index) is not int or index != 0:
        raise ValueError("declared global phase with index0 required")
    return root/phase


def require_new_folder(manifest, phase, index):
    folder = point_folder(manifest, phase, index)
    run = Path(manifest["run"])
    if not folder.is_relative_to(run) or folder.resolve() != folder:
        raise ValueError("never dispatch or write through an original/aliased folder")
    if folder.exists() and any(folder.iterdir()):
        raise ValueError("retain an occupied/failed/interrupted folder; fresh output required")
    return folder


def verify_inventory(original, inventory):
    original = Path(original).resolve()
    expected = {}
    for row in inventory:
        path = Path(row["path"])
        if (not path.is_absolute() or not path.is_relative_to(original) or
                path.resolve() != path or str(path) in expected):
            raise ValueError("invalid original inventory path")
        expected[str(path)] = row
    current = {str(p) for p in original.rglob("*") if p.is_file()}
    if current != set(expected):
        raise ValueError("original file inventory changed")
    for name, row in expected.items():
        path = Path(name)
        if path.stat().st_size != row["bytes"] or file_hash(path) != row["sha256"]:
            raise ValueError("original file changed: "+name)


def owner_absent():
    """An unrelated PID reuse is not the old owner; never signal anything."""
    if not Path("/proc/sys/kernel/random/boot_id").is_file():
        raise ValueError("continuation requires Linux process observations")
    matches = []
    for proc in Path("/proc").glob("[0-9]*"):
        try:
            args = (proc/"cmdline").read_bytes().decode(errors="replace").split("\0")
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        if "attempt-2" in args and any(a.endswith("/run-historical-sidis-jobs-real-data.py") for a in args):
            matches.append(proc.name)
    if matches:
        raise ValueError("the original attempt2 owner is live; do not duplicate it")
    return {"boot_id": Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
            "matching_original_owners": [], "original_processes_signaled": False}


def verify_manifest(path, frozen, adapter, api, *, expected_sha=None):
    path, frozen, adapter = map(lambda p: Path(p).resolve(), (path, frozen, adapter))
    if expected_sha is not None and file_hash(path) != expected_sha:
        raise ValueError("changed externally bound continuation manifest")
    m = read_json(path)
    c = contract(adapter)
    if (m["status"] != "PREPARED_IMMUTABLE_CONTINUATION_NOT_A_GATE" or m["scope"] != SCOPE or
            any(m.get(f) is not False for f in FLAGS) or m["contract_sha256"] != CONTRACT_SHA or
            m["frozen"] != str(frozen) or m["adapter"] != str(adapter)):
        raise ValueError("wrong or self-admitting continuation manifest")
    old_head, old_files = verify_git(frozen, FROZEN_HEAD)
    new_head, new_files = verify_git(adapter, m["adapter_commit"])
    if (len(old_files) != 335 or m["frozen_files"] != old_files or m["adapter_files"] != new_files or
            any(new_files.get(k) != h for k, h in old_files.items())):
        raise ValueError("source file bindings changed")
    if (old_head != c["frozen_numerical_head"] or m["frozen_head"] != old_head or new_head == POLICY_COMMIT or
            m["numerical_source_fingerprint"] != c["frozen_numerical_fingerprint"]):
        raise ValueError("implementation must follow the prospective contract")
    if api.ROOT != frozen or api.source_fingerprint(frozen) != c["frozen_numerical_fingerprint"]:
        raise ValueError("untrusted frozen numerical package")
    run, original = Path(m["run"]), Path(m["original"])
    if (run != path.parent or run.resolve() != run or not run.is_relative_to(frozen/"runs") or
            run.is_relative_to(original) or original.is_relative_to(run) or
            original != frozen/"runs/historical-sidis-jobs-real-data/attempt-2"):
        raise ValueError("fresh disjoint canonical ignored continuation root required")
    origin = c["supported_origin"]
    for key, sha in (("review", origin["independent_prefix_review_sha256"]),
                     ("integrity", origin["post_review_integrity_sha256"]),
                     ("loss", origin["original_loss_inventory_sha256"])):
        if file_hash(m[key]) != sha:
            raise ValueError("changed pinned original review/loss input")
    loss = read_json(m["loss"])
    workspace = frozen.parents[2]
    inventory = [dict(r, path=str(workspace/r["path"]))
                 for r in loss["interrupted_attempts"]["sidis132"]["retained_files"]]
    if len(inventory) != 141 or m["original_inventory"] != inventory:
        raise ValueError("the exact original141-file inventory is required")
    verify_inventory(original, inventory)
    if file_hash(original/"invocation.json") != origin["outer_invocation_sha256"]:
        raise ValueError("changed original outer invocation")
    deps = read_json(original/"invocation.json")["dependency_receipts"]
    if deps != m["dependency_receipts"] or {m["native"], m["parton"]} != set(deps):
        raise ValueError("changed original build receipt identities")
    if any(file_hash(p) != h for p, h in deps.items()):
        raise ValueError("changed frozen native/provider receipt")
    request_path = original/"full/request.json"
    request = api.verify_request(request_path)
    if str(request_path) != m["request"] or file_hash(request_path) != m["request_sha256"]:
        raise ValueError("changed exact original source request")
    for i in range(1, PREFIX+1):
        api.verify_metadata(request_path, m["native"], m["parton"], original/"full", "native123", i)
    if request["source_fingerprint"] != c["frozen_numerical_fingerprint"]:
        raise ValueError("wrong numerical source fingerprint")
    owner_absent()
    return m


def prepare(frozen, adapter, original, review, integrity, loss, run, native, parton, api):
    frozen, adapter, original, run = map(lambda p: Path(p).resolve(), (frozen, adapter, original, run))
    c = contract(adapter)
    old_head, old_files = verify_git(frozen, FROZEN_HEAD)
    new_head, new_files = verify_git(adapter)
    if subprocess.run(["git", "merge-base", "--is-ancestor", POLICY_COMMIT, new_head], cwd=adapter).returncode:
        raise ValueError("prospective policy is not an implementation ancestor")
    if run.exists() or not run.is_relative_to(frozen/"runs") or run.is_relative_to(original):
        raise ValueError("a fresh disjoint ignored run directory is required")
    workspace = frozen.parents[2]
    inventory = [dict(r, path=str(workspace/r["path"]))
                 for r in read_json(loss)["interrupted_attempts"]["sidis132"]["retained_files"]]
    m = dict(status="PREPARED_IMMUTABLE_CONTINUATION_NOT_A_GATE", scope=SCOPE,
             contract_sha256=CONTRACT_SHA, frozen=str(frozen), frozen_head=old_head,
             numerical_source_fingerprint=c["frozen_numerical_fingerprint"], frozen_files=old_files,
             adapter=str(adapter), adapter_commit=new_head, adapter_files=new_files,
             original=str(original), original_inventory=inventory, run=str(run),
             request=str(original/"full/request.json"), request_sha256=file_hash(original/"full/request.json"),
             native=str(Path(native).resolve()), parton=str(Path(parton).resolve()),
             dependency_receipts=read_json(original/"invocation.json")["dependency_receipts"],
             review=str(Path(review).resolve()), integrity=str(Path(integrity).resolve()),
             loss=str(Path(loss).resolve()), origin_process_observation=owner_absent(),
             **dict.fromkeys(FLAGS, False))
    run.mkdir(parents=True, exist_ok=False)
    path = run/"continuation.json"
    atomic_json(path, m)
    verify_manifest(path, frozen, adapter, api)
    return path


def provenance(m, path):
    return {"scope": SCOPE, "manifest_sha256": file_hash(path),
            "adapter_commit": m["adapter_commit"], "contract_sha256": CONTRACT_SHA}


def _selection(m, phase, api):
    return api._selection(m["run"], phase, api.read_json(m["request"]), m["numerical_source_fingerprint"])


def input_paths(m, path, phase, index, api):
    """Path-only equivalent of the frozen controller, with explicit transport bindings."""
    frozen = Path(m["frozen"])
    paths = [Path(m[k]) for k in ("request", "native", "parton")]
    paths += [frozen/"config/historical-sidis-job-policy.json", api.TIME_BINARY,
              frozen/"scripts/run_historical_sidis_job.py", Path(path)]
    paths += [Path(m["adapter"])/name for name in m["adapter_files"]]
    def points(name):
        paths.extend(point_folder(m, name, i)/"invocation.json" for i in range(1, 133))
    def glob(name, suffix="invocation.json"):
        paths.append(point_folder(m, name)/suffix)
    if phase in POINT_PHASES:
        if index > 1:
            paths.append(point_folder(m, phase, index-1)/"invocation.json")
        if phase == "native4":
            glob("native-select3", "report.json")
            paths.append(point_folder(m, "native123", index)/"report.json")
        if phase in ("tailref123", "reference4"):
            level = _selection(m, "native-selection", api)["selected_level"]
            glob("native-selection", "report.json")
            paths.append(point_folder(m, "native123" if level == 3 else "native4", index)/"report.json")
        if phase == "reference4":
            glob("derivative-select3", "report.json")
            paths.append(point_folder(m, "tailref123", index)/"report.json")
    elif phase != "source-math":
        points("native123")
        if phase != "native-select3":
            glob("native-select3")
            if _selection(m, "native-select3", api)["selected_level"] == 0:
                points("native4")
        if phase in ("derivative-select3", "diagnostic"):
            glob("native-selection")
            points("tailref123")
        if phase == "diagnostic":
            glob("derivative-select3")
            glob("source-math")
            if _selection(m, "derivative-select3", api)["selected_reference_level"] == 0:
                points("reference4")
    else:
        glob("native-selection")
    paths.append(Path(m["run"])/"prefix-review/invocation.json")
    return sorted({p.resolve() for p in paths})


def verify_job(m, path, phase, index, api):
    folder = point_folder(m, phase, index)
    if phase == "native123" and index <= PREFIX:
        return api.verify_metadata(m["request"], m["native"], m["parton"], Path(m["original"])/"full", phase, index)[3]
    inv = read_json(folder/"invocation.json")
    final = phase == "diagnostic"
    if (inv.get("status") != api.COMPLETED or inv.get("phase") != phase or
            type(inv.get("point_index")) is not int or inv["point_index"] != index or
            inv.get("continuation") != provenance(m, path) or
            inv.get("source_fingerprint") != m["numerical_source_fingerprint"] or
            inv.get("policy_sha256") != api.POLICY_SHA or not api._flags(inv, final) or
            inv.get("time_binary") != str(api.TIME_BINARY) or not inv.get("completed_at")):
        raise ValueError("missing, failed, stale or self-admitting continuation job")
    expected = {str(p): file_hash(p) for p in input_paths(m, path, phase, index, api)}
    if inv["input_hashes"] != expected or set(inv["outputs"]) != api._products(folder, final):
        raise ValueError("changed complete continuation dependencies/outputs")
    if any(file_hash(p) != h for p, h in inv["outputs"].items()):
        raise ValueError("changed complete continuation output")
    if (file_hash(folder/"initial-invocation.json") != inv["initial_record_sha256"] or
            read_json(folder/"initial-invocation.json")["record_role"] != "initial_sentinel_not_an_executed_failure" or
            [s["name"] for s in inv["steps"]] != ["producer", "audit"]):
        raise ValueError("changed sentinel or incomplete command ledger")
    budget = {"resources": api.load_policy()["resources"][phase if phase in POINT_PHASES else "global"]}
    for s in inv["steps"]:
        api._require_command_success(s, budget)
        if (file_hash(folder/(s["name"]+".log")) != s["log_sha256"] or
                file_hash(folder/(s["name"]+".resources.txt")) != s["resource_sha256"] or
                api.parse_resource_file(folder/(s["name"]+".resources.txt")) != s["resource"]):
            raise ValueError("changed complete-command measurements")
    if sum(s["wrapper_elapsed_seconds"] for s in inv["steps"]) > budget["resources"]["whole_Julia_jobs_time_budget_seconds"]:
        raise ValueError("whole-job resource budget exhausted")
    report, audit = api._report_audit(folder, phase, index, read_json(m["request"]),
                                     m["numerical_source_fingerprint"], m["request_sha256"])
    if report.get("continuation") != provenance(m, path) or audit.get("continuation") != provenance(m, path):
        raise ValueError("unbound continuation producer/reader")
    if final and (folder/"point-comparison.csv").read_text() != api.comparison_csv(report):
        raise ValueError("changed comparison CSV")
    return inv


def _julia_command(m, path, julia, program, *args):
    return [julia, "--startup-file=no", "--project="+str(Path(m["frozen"])/"numerics"),
            str(Path(m["adapter"])/"numerics/bin"/program), m["frozen"], str(path), file_hash(path), *map(str, args)]


@contextmanager
def job_lock(m, api):
    import fcntl
    # The existing shared checkout lock is opened read-only; no old file is created or touched.
    shared = Path(m["frozen"])/".runtime/sidis-pilot.lock"
    with shared.open("r") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError("another whole Julia job owns the frozen checkout") from error
        try:
            with api._pilot_lock(Path(m["run"])):
                yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def verify_prefix(m, path, api):
    folder = Path(m["run"])/"prefix-review"
    inv = read_json(folder/"invocation.json")
    if inv.get("status") != "ORIGINAL13_FRESHLY_REVIEWED_BY_TYPED_PATH_READER_NOT_DIAGNOSTIC_GATE":
        raise ValueError("fresh complete Julia path-resolver prefix review required")
    if inv.get("continuation") != provenance(m, path) or any(inv.get(f) is not False for f in FLAGS):
        raise ValueError("stale/self-admitting prefix review")
    budget = {"resources": api.load_policy()["resources"]["global"]}
    api._require_command_success(inv["step"], budget)
    expected = {str(folder/name) for name in ("audit.json", "prefix-audit.log", "prefix-audit.resources.txt")}
    step = inv["step"]
    if (set(inv["outputs"]) != expected or step["name"] != "prefix-audit" or
            file_hash(folder/"prefix-audit.log") != step["log_sha256"] or
            file_hash(folder/"prefix-audit.resources.txt") != step["resource_sha256"] or
            api.parse_resource_file(folder/"prefix-audit.resources.txt") != step["resource"] or
            step["wrapper_elapsed_seconds"] > budget["resources"]["whole_Julia_jobs_time_budget_seconds"]):
        raise ValueError("changed prefix command/resource inventory")
    for name, digest in inv["outputs"].items():
        if file_hash(name) != digest:
            raise ValueError("changed prefix review output")
    audit = read_json(folder/"audit.json")
    original_hashes = {str(point_folder(m, "native123", i)/"invocation.json"):
                       file_hash(point_folder(m, "native123", i)/"invocation.json") for i in range(1, 14)}
    if (audit["status"] != "ORIGINAL13_TYPED_PATH_ARITHMETIC_REVIEWED_NOT_DIAGNOSTIC_GATE" or
            audit["continuation"] != provenance(m, path) or audit["points_reviewed"] != list(range(1, 14)) or
            audit["original_invocations"] != original_hashes or
            audit["source_fingerprint"] != m["numerical_source_fingerprint"] or
            audit["request_sha256"] != m["request_sha256"] or any(audit.get(f) is not False for f in FLAGS)):
        raise ValueError("incomplete original prefix reader")
    return inv


def audit_prefix(m, path, api, julia):
    folder = Path(m["run"])/"prefix-review"
    folder.mkdir(exist_ok=False)
    rec = {"status": "FAIL", "continuation": provenance(m, path), **dict.fromkeys(FLAGS, False)}
    try:
        budget = {"resources": api.load_policy()["resources"]["global"]}
        with job_lock(m, api):
            row = api._timed_job("prefix-audit", _julia_command(m, path, julia,
                "audit_sidis_continuation_prefix.jl", folder/"audit.json"), Path(m["frozen"]), folder,
                api.TIME_BINARY, budget, api._deadline(budget["resources"]["whole_Julia_jobs_time_budget_seconds"], None))
        rec["step"] = row
        api._require_command_success(row, budget)
        verify_manifest(path, m["frozen"], m["adapter"], api)
        rec["outputs"] = {str(p): file_hash(p) for p in (folder/"audit.json", folder/"prefix-audit.log", folder/"prefix-audit.resources.txt")}
        rec["status"] = "ORIGINAL13_FRESHLY_REVIEWED_BY_TYPED_PATH_READER_NOT_DIAGNOSTIC_GATE"
        atomic_json(folder/"invocation.json", rec)
        verify_prefix(m, path, api)
    except BaseException as error:
        rec["status"] = "FAIL"
        rec["error"] = type(error).__name__+": "+str(error)
        raise
    finally:
        rec["completed_at"] = api.utc_now()
        atomic_json(folder/"invocation.json", rec)
    return verify_prefix(m, path, api)


def run_job(m, path, phase, index, api, julia, timeout_seconds=None):
    verify_manifest(path, m["frozen"], m["adapter"], api)
    verify_prefix(m, path, api)
    folder = require_new_folder(m, phase, index)
    if phase in POINT_PHASES and index > 1:
        verify_job(m, path, phase, index-1, api)
    if phase == "source-math":
        verify_job(m, path, "native-selection", 0, api)
    folder.mkdir(parents=True, exist_ok=True)
    inv = {"status": "FAIL", "scope": api.SCOPE, "phase": phase, "point_index": index,
           "started_at": api.utc_now(), "steps": [], "outputs": {}, "continuation": provenance(m, path),
           **dict.fromkeys(FLAGS, False)}
    atomic_json(folder/"initial-invocation.json", {**inv, "record_role": "initial_sentinel_not_an_executed_failure"})
    inv["initial_record_sha256"] = file_hash(folder/"initial-invocation.json")
    atomic_json(folder/"invocation.json", inv)
    try:
        paths = input_paths(m, path, phase, index, api)
        inv.update(source_fingerprint=m["numerical_source_fingerprint"], code_commit=FROZEN_HEAD,
                   policy_sha256=api.POLICY_SHA, time_binary=str(api.TIME_BINARY),
                   input_hashes={str(p): file_hash(p) for p in paths},
                   time_version=subprocess.check_output([str(api.TIME_BINARY), "--version"], text=True, timeout=10))
        budget = {"resources": api.load_policy()["resources"][phase if phase in POINT_PHASES else "global"]}
        program = "run_sidis_continuation_point.jl" if phase in POINT_PHASES else "run_sidis_continuation_global.jl"
        producer = _julia_command(m, path, julia, program, phase, index, folder/"report.json")
        auditor = _julia_command(m, path, julia, "audit_sidis_continuation_job.jl", phase, index, folder/"report.json", folder/"audit.json")
        with job_lock(m, api):
            deadline = api._deadline(budget["resources"]["whole_Julia_jobs_time_budget_seconds"], timeout_seconds)
            for name, command in (("producer", producer), ("audit", auditor)):
                row = api._timed_job(name, command, Path(m["frozen"]), folder, api.TIME_BINARY, budget, deadline)
                inv["steps"].append(row)
                atomic_json(folder/"invocation.json", inv)
                api._require_command_success(row, budget)
        api._report_audit(folder, phase, index, read_json(m["request"]), m["numerical_source_fingerprint"], m["request_sha256"])
        if any(read_json(folder/name).get("continuation") != provenance(m, path) for name in ("report.json", "audit.json")):
            raise ValueError("unbound continuation producer/reader")
        if phase == "diagnostic":
            api.render_comparison(read_json(folder/"report.json"), folder)
        verify_manifest(path, m["frozen"], m["adapter"], api)
        if any(file_hash(p) != h for p, h in inv["input_hashes"].items()):
            raise ValueError("job input changed during execution")
        inv["outputs"] = {str(p): file_hash(p) for p in api._products(folder, phase == "diagnostic")}
        inv["status"] = api.COMPLETED
        if phase == "diagnostic":
            inv.update(scope=api.SCIENTIFIC_SCOPE, accepted_prediction_diagnostic=True, all132_coverage_verified=True)
        inv["completed_at"] = api.utc_now()
        atomic_json(folder/"invocation.json", inv)
        verify_job(m, path, phase, index, api)
    except BaseException as error:
        inv["status"] = "FAIL"
        inv.update(dict.fromkeys(FLAGS, False))
        inv["error"] = type(error).__name__+": "+str(error)
        raise
    finally:
        inv["completed_at"] = api.utc_now()
        atomic_json(folder/"invocation.json", inv)
    return verify_job(m, path, phase, index, api)


def audit_cached(m, path, phase, index, api, julia):
    verify_manifest(path, m["frozen"], m["adapter"], api)
    original = verify_job(m, path, phase, index, api)
    review = Path(m["run"])/"cache-reviews"/(phase+f"-{index}-"+uuid.uuid4().hex)
    review.mkdir(parents=True, exist_ok=False)
    rec = {"status": "FAIL", "continuation": provenance(m, path), "phase": phase, "point_index": index,
           "original_invocation_sha256": file_hash(point_folder(m, phase, index)/"invocation.json"),
           **dict.fromkeys(FLAGS, False)}
    try:
        budget = {"resources": api.load_policy()["resources"]["global"]}
        command = _julia_command(m, path, julia, "audit_sidis_continuation_job.jl", phase, index,
                                 point_folder(m, phase, index)/"report.json", review/"audit.json")
        with job_lock(m, api):
            row = api._timed_job("cache-audit", command, Path(m["frozen"]), review, api.TIME_BINARY,
                                 budget, api._deadline(budget["resources"]["whole_Julia_jobs_time_budget_seconds"], None))
        rec["step"] = row
        api._require_command_success(row, budget)
        audit = read_json(review/"audit.json")
        if (audit.get("status") != api.AUDITED or audit.get("continuation") != provenance(m, path) or
                audit.get("report_sha256") != file_hash(point_folder(m, phase, index)/"report.json") or
                audit.get("source_fingerprint") != m["numerical_source_fingerprint"] or
                audit.get("request_sha256") != m["request_sha256"] or audit.get("policy_sha256") != api.POLICY_SHA or
                not api._flags(audit, phase == "diagnostic")):
            raise ValueError("failed fresh cached Julia review")
        verify_manifest(path, m["frozen"], m["adapter"], api)
        if original != verify_job(m, path, phase, index, api):
            raise ValueError("original invocation changed during cached review")
        rec["status"] = "CURRENT_CONTINUATION_ARITHMETIC_REVIEWED_NOT_FITTING"
    except BaseException as error:
        rec["error"] = type(error).__name__+": "+str(error)
        raise
    finally:
        rec["completed_at"] = api.utc_now()
        atomic_json(review/"review.json", rec)
    return rec


def run_all(m, path, api, julia):
    """Serial complete jobs; a checkpoint decision never selects a point subset."""
    def complete(phase, index=0):
        folder = point_folder(m, phase, index)
        if (folder/"invocation.json").is_file():
            audit_cached(m, path, phase, index, api, julia)
        else:
            run_job(m, path, phase, index, api, julia)
        print(json.dumps({"phase": phase, "index": index, "status": "COMPLETE_NOT_FITTING"}), flush=True)
    if (Path(m["run"])/"prefix-review").exists():
        verify_prefix(m, path, api)
    else:
        audit_prefix(m, path, api, julia)
    for i in range(PREFIX+1, 133):
        complete("native123", i)
    complete("native-select3")
    if _selection(m, "native-select3", api)["selected_level"] == 0:
        for i in range(1, 133):
            complete("native4", i)
    complete("native-selection")
    complete("source-math")
    for i in range(1, 133):
        complete("tailref123", i)
    complete("derivative-select3")
    if _selection(m, "derivative-select3", api)["selected_reference_level"] == 0:
        for i in range(1, 133):
            complete("reference4", i)
    complete("diagnostic")
    return audit_cached(m, path, "diagnostic", 0, api, julia)
