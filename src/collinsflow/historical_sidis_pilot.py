"""Whole Julia resource-pilot jobs; no Python predictions or fitting admission."""
from __future__ import annotations

from contextlib import contextmanager
import os
from pathlib import Path
import subprocess
import tempfile
import time

from .experimental import acquire
from .historical_points import finite_number
from .historical_sidis_points import FALSE_FLAGS as SOURCE_FALSE_FLAGS, verify_request
from .job_control import run_complete_process
from .provenance import atomic_json, file_hash, git_commit, read_json, source_fingerprint, utc_now

SCOPE = "historical_reconstructed_sidis_resource_pilot"
POLICY = "config/historical-sidis-pilot-policy.json"
POLICY_SHA = "395854b73f595fb361141f3993320ba5c127ff313dc0943e261346ae9ea19f4a"
ROOT = Path(__file__).resolve().parents[2]
FALSE_FLAGS = (*SOURCE_FALSE_FLAGS, "accepted_prediction_diagnostic", "accepted_fitting_gate", "all132_coverage_verified")
COMPLETED = "PILOT_COMPLETED_NOT_DIAGNOSTIC_GATE"
AUDITED = "PILOT_REVIEWED_NOT_DIAGNOSTIC_GATE"
TIME_FORMAT = "COLLINSFLOW_TIME_V1\t%M\t%e\t%U\t%S\t%x"


def load_policy(root):
    root = Path(root).resolve()
    if root != ROOT:
        raise ValueError("pilot execution requires the current trusted project source")
    if file_hash(root/POLICY) != POLICY_SHA:
        raise ValueError("changed frozen SIDIS pilot policy")
    policy = read_json(root/POLICY)
    if policy["scope"] != SCOPE or any(policy.get(k) is not False for k in FALSE_FLAGS):
        raise ValueError("pilot policy cannot approve a diagnostic or fitting")
    if file_hash(root/policy["parent_contract"]) != policy["parent_contract_sha256"]:
        raise ValueError("changed SIDIS pilot parent contract")
    return policy


def prepare_primary_cache(manifest_path, cache, output, *, offline=False):
    """Acquire the nine pinned convention files as a complete source job."""
    manifest_path, cache, output = Path(manifest_path).resolve(), Path(cache).resolve(), Path(output).resolve()
    root = manifest_path.parent.parent
    load_policy(root)
    parent = read_json(root/"config/historical-sidis-point-contract.json")
    if file_hash(manifest_path) != parent["primary_manifest_sha256"]:
        raise ValueError("changed convention manifest")
    if output.exists():
        raise ValueError("preserve earlier primary acquisition records; use a fresh receipt")
    record = {"status": "FAIL", "scope": "primary_convention_bytes_not_data_admission",
              **{k: False for k in FALSE_FLAGS}, "started_at": utc_now(), "inputs": {}}
    atomic_json(output, record)
    try:
        manifest = read_json(manifest_path)
        if len(manifest["sources"]) != 9:
            raise ValueError("nine locked primary files required")
        for source in manifest["sources"]:
            path = acquire(source, cache, offline=offline)
            record["inputs"][str(path)] = file_hash(path)
            atomic_json(output, record)
        record["status"] = "PINNED_PRIMARY_BYTES_VERIFIED"
        record["source_fingerprint"] = source_fingerprint(root)
    except BaseException as error:
        record["error"] = type(error).__name__+": "+str(error)
        raise
    finally:
        record["completed_at"] = utc_now()
        atomic_json(output, record)
    return record


def available_memory_bytes():
    if os.name != "posix" or not Path("/proc/meminfo").is_file():
        raise ValueError("SIDIS resource pilots require Linux memory observations")
    fields = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
    return int(fields["MemAvailable"].split()[0])*1024


@contextmanager
def _pilot_lock(root):
    import fcntl
    lock = Path(root)/".runtime/sidis-pilot.lock"
    lock.parent.mkdir(parents=True, exist_ok=True)
    with lock.open("a+") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError("another complete SIDIS pilot owns this checkout") from error
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def parse_resource_file(path):
    """Read the retained GNU time record, keeping its limited measurement scope."""
    if Path(path).stat().st_size > 20_000:
        raise ValueError("oversized complete-command resource record")
    lines = [line for line in Path(path).read_text(encoding="utf-8").splitlines()
             if line.startswith("COLLINSFLOW_TIME_V1\t")]
    if len(lines) != 1 or len(lines[0].split("\t")) != 6:
        raise ValueError("missing or ambiguous complete-command resource record")
    _, rss, elapsed, user, system, code = lines[0].split("\t")
    result = {"maximum_RSS_KiB": int(rss), "wall_seconds": float(elapsed),
              "user_seconds": float(user), "system_seconds": float(system), "child_exit_code": int(code)}
    if result["maximum_RSS_KiB"] <= 0 or result["child_exit_code"] < 0:
        raise ValueError("invalid complete-command resource record")
    for key in ("wall_seconds", "user_seconds", "system_seconds"):
        finite_number(result[key], key, nonnegative=True)
    result["scope"] = "GNU_time_direct_child_maximum_not_simultaneous_process_tree_peak"
    return result


def _case(request, policy, case_id):
    cases = [case for case in policy["selection"]["cases"] if case["id"] == case_id]
    if len(cases) != 1:
        raise ValueError("case is outside frozen pilot selection")
    case = cases[0]
    points = [p for p in request["points"] if p["id"] == case["point_id"]]
    if (len(request["points"]) != 132 or len(points) != 1 or
            any(points[0][k] != case[k] for k in ("publication", "target"))):
        raise ValueError("changed pilot case or incomplete full source request")
    return case


def _environment():
    return {**os.environ, **{k: "1" for k in ("JULIA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")}}


def _julia_command(julia, root, program, *arguments):
    return [julia, "--startup-file=no", "--project="+str(root/"numerics"),
            str(root/"numerics/bin"/program), *map(str, arguments)]


def _timed_job(name, command, root, folder, tool, policy, deadline):
    remaining = deadline-time.monotonic()
    if remaining <= 0:
        raise ValueError("whole Julia-jobs time budget exhausted before "+name)
    if available_memory_bytes() < policy["resources"]["minimum_available_memory_bytes"]:
        raise ValueError("insufficient available memory before complete pilot command")
    log, resource = folder/(name+".log"), folder/(name+".resources.txt")
    measured = [str(tool), "--output", str(resource), "--format", TIME_FORMAT, "--", *command]
    started = time.monotonic()
    with log.open("w", encoding="utf-8", newline="\n") as stream:
        result = run_complete_process(measured, cwd=root, env=_environment(), stdout=stream, timeout_seconds=remaining)
    row = {"name": name, "command": measured, **result._asdict(),
           "wrapper_elapsed_seconds": time.monotonic()-started, "log_sha256": file_hash(log)}
    if resource.exists():
        row["resource_sha256"] = file_hash(resource)
        try:
            row["resource"] = parse_resource_file(resource)
        except ValueError as error:
            row["resource_error"] = str(error)
    return row


def _require_command_success(row, policy):
    if (row["returncode"] != 0 or row["process_returncode"] != 0 or row["interruption"] is not None or
            row["termination_signal"] is not None):
        raise ValueError("interrupted or failed complete pilot "+row["name"])
    resource = row.get("resource")
    if resource is None or resource["child_exit_code"] != 0:
        raise ValueError("missing or failed measured Julia command")
    if resource["maximum_RSS_KiB"] > policy["resources"]["maximum_command_RSS_KiB"]:
        raise ValueError("complete-command RSS resource budget exceeded")


def _check_audit(request_path, report_path, audit_path, case_id, fingerprint):
    report, audit = read_json(report_path), read_json(audit_path)
    for item, status in ((report, COMPLETED), (audit, AUDITED)):
        if (item.get("status") != status or item.get("scope") != SCOPE or item.get("case_id") != case_id or
                item.get("source_fingerprint") != fingerprint or any(item.get(k) is not False for k in FALSE_FLAGS)):
            raise ValueError("incomplete, stale or self-admitting pilot result")
    if (audit["request_sha256"] != file_hash(request_path) or audit["report_sha256"] != file_hash(report_path) or
            audit["policy_sha256"] != POLICY_SHA or report["policy_sha256"] != POLICY_SHA or
            audit["saved_scalar_checks_recomputed"] != 4):
        raise ValueError("changed pilot audit bindings")
    return report, audit


def run_pilot_job(request_path, case_id, native_build, parton_build, output, *, julia="julia", time_binary="/usr/bin/time"):
    """Schedule two complete Julia jobs; all mathematical evaluation stays in Julia."""
    request_path, native_build, parton_build, output = map(
        lambda p: Path(p).resolve(), (request_path, native_build, parton_build, output))
    if output.exists() and any(output.iterdir()):
        raise ValueError("use a fresh pilot output directory; retain earlier artifacts")
    output.mkdir(parents=True, exist_ok=True)
    invocation = {"status": "FAIL", "scope": SCOPE, "case_id": case_id, **{k: False for k in FALSE_FLAGS},
                  "started_at": utc_now(), "steps": [], "outputs": {}}
    initial = output/"initial-invocation.json"
    atomic_json(initial, {**invocation, "record_role": "initial_sentinel_not_an_executed_failure"})
    invocation["initial_record_sha256"] = file_hash(initial)
    record_path = output/"invocation.json"
    atomic_json(record_path, invocation)
    try:
        request = verify_request(request_path)
        root = Path(request["source_preparation"]["manifest"]).parent.parent
        policy = load_policy(root)
        _case(request, policy, case_id)
        fingerprint = source_fingerprint(root)
        if request["source_fingerprint"] != fingerprint:
            raise ValueError("stale pilot source request")
        tool = Path(time_binary).resolve()
        if os.name != "posix" or not tool.is_file():
            raise ValueError("Linux GNU time required before dispatch")
        version = subprocess.run([str(tool), "--version"], capture_output=True, text=True, check=True, timeout=10).stdout
        inputs = {str(p): file_hash(p) for p in (request_path, native_build, parton_build, root/POLICY, tool)}
        invocation.update(source_fingerprint=fingerprint, code_commit=git_commit(root), input_hashes=inputs,
                          time_binary=str(tool), time_version=version, policy_sha256=POLICY_SHA)
        report_path, audit_path = output/"pilot.json", output/"audit.json"
        commands = [("producer", _julia_command(julia, root, "pilot_historical_sidis_point.jl",
                      request_path, case_id, native_build, parton_build, report_path)),
                    ("audit", _julia_command(julia, root, "audit_historical_sidis_pilot.jl",
                      request_path, case_id, report_path, audit_path))]
        with _pilot_lock(root):
            deadline = time.monotonic()+policy["resources"]["whole_Julia_jobs_time_budget_seconds"]
            for name, command in commands:
                row = _timed_job(name, command, root, output, tool, policy, deadline)
                invocation["steps"].append(row)
                atomic_json(record_path, invocation)
                _require_command_success(row, policy)
        _check_audit(request_path, report_path, audit_path, case_id, fingerprint)
        verify_request(request_path)
        if (source_fingerprint(root) != fingerprint or any(file_hash(p) != h for p, h in inputs.items()) or
                file_hash(initial) != invocation["initial_record_sha256"]):
            raise ValueError("pilot input/code/initial record changed during execution")
        products = [report_path, audit_path, *[output/(name+suffix) for name in ("producer", "audit")
                                             for suffix in (".log", ".resources.txt")]]
        invocation["outputs"] = {str(p): file_hash(p) for p in products}
        invocation["status"] = COMPLETED
    except BaseException as error:
        invocation["error"] = type(error).__name__+": "+str(error)
        raise
    finally:
        invocation["completed_at"] = utc_now()
        atomic_json(record_path, invocation)
    return {"status": COMPLETED, "scope": SCOPE, "case_id": case_id, **{k: False for k in FALSE_FLAGS}}


def verify_pilot_job(request_path, case_id, output, *, julia="julia", time_binary="/usr/bin/time"):
    """A marker alone cannot replace source replay, measured records and a fresh Julia audit."""
    request_path, output = Path(request_path).resolve(), Path(output).resolve()
    request = verify_request(request_path)
    root = Path(request["source_preparation"]["manifest"]).parent.parent
    policy = load_policy(root)
    _case(request, policy, case_id)
    fingerprint = source_fingerprint(root)
    invocation = read_json(output/"invocation.json")
    if (invocation.get("status") != COMPLETED or invocation.get("scope") != SCOPE or invocation.get("case_id") != case_id or
            any(invocation.get(k) is not False for k in FALSE_FLAGS) or invocation.get("source_fingerprint") != fingerprint):
        raise ValueError("stale, missing or self-admitting cached pilot")
    if invocation["input_hashes"].get(str(request_path)) != file_hash(request_path):
        raise ValueError("cached pilot request identity changed")
    report = read_json(output/"pilot.json")
    tool = Path(time_binary).resolve()
    expected_inputs = {str(p) for p in (request_path, Path(report["native_build_receipt"]).resolve(),
        Path(report["parton_build_receipt"]).resolve(), root/POLICY, tool)}
    if set(invocation["input_hashes"]) != expected_inputs or invocation["time_binary"] != str(tool):
        raise ValueError("incomplete input set or changed caller-selected timing tool")
    for p, h in invocation["input_hashes"].items():
        if file_hash(p) != h:
            raise ValueError("changed pilot input")
    if file_hash(output/"initial-invocation.json") != invocation["initial_record_sha256"]:
        raise ValueError("changed initial pilot record")
    expected = {str(output/name) for name in ("pilot.json", "audit.json", "producer.log", "producer.resources.txt", "audit.log", "audit.resources.txt")}
    if set(invocation["outputs"]) != expected or any(file_hash(p) != h for p, h in invocation["outputs"].items()):
        raise ValueError("incomplete or changed pilot outputs")
    if [s["name"] for s in invocation["steps"]] != ["producer", "audit"]:
        raise ValueError("incomplete pilot command ledger")
    for row in invocation["steps"]:
        _require_command_success(row, policy)
        if (file_hash(output/(row["name"]+".log")) != row["log_sha256"] or
                file_hash(output/(row["name"]+".resources.txt")) != row["resource_sha256"]):
            raise ValueError("changed pilot command log or resource bytes")
        if parse_resource_file(output/(row["name"]+".resources.txt")) != row["resource"]:
            raise ValueError("changed measured pilot resources")
    _check_audit(request_path, output/"pilot.json", output/"audit.json", case_id, fingerprint)
    with _pilot_lock(root), tempfile.TemporaryDirectory(prefix="sidis-pilot-audit-", dir=output.parent) as directory:
        folder = Path(directory)
        audit = folder/"audit.json"
        command = _julia_command(julia, root, "audit_historical_sidis_pilot.jl", request_path, case_id, output/"pilot.json", audit)
        deadline = time.monotonic()+policy["resources"]["whole_Julia_jobs_time_budget_seconds"]
        row = _timed_job("cache-audit", command, root, folder, tool, policy, deadline)
        try:
            _require_command_success(row, policy)
            _check_audit(request_path, output/"pilot.json", audit, case_id, fingerprint)
        except BaseException:
            failure = output/("cache-audit-failure-"+utc_now().replace(":", "-")+".json")
            atomic_json(failure, {"status": "FAIL", "scope": SCOPE, **{k: False for k in FALSE_FLAGS},
                "step": row, "audit": read_json(audit) if audit.exists() else None,
                "log": (folder/"cache-audit.log").read_text(encoding="utf-8")})
            raise
    return {"status": "CACHED_PILOT_REVIEWED_NOT_DIAGNOSTIC_GATE", "scope": SCOPE,
            "case_id": case_id, **{k: False for k in FALSE_FLAGS}}


def _measurement_ledger(request_path, run, *, julia="julia", time_binary="/usr/bin/time"):
    """Resource arithmetic and job identity only; all physics arithmetic stays in Julia."""
    request_path, run = Path(request_path).resolve(), Path(run).resolve()
    request = verify_request(request_path)
    policy = load_policy(ROOT)
    cases = []
    for case in policy["selection"]["cases"]:
        output = run/"cases"/case["id"]
        verify_pilot_job(request_path, case["id"], output, julia=julia, time_binary=time_binary)
        invocation = read_json(output/"invocation.json")
        report = read_json(output/"pilot.json")
        cases.append({"case_id": case["id"], "point_id": case["point_id"],
            "invocation_sha256": file_hash(output/"invocation.json"), "report_sha256": file_hash(output/"pilot.json"),
            "audit_sha256": file_hash(output/"audit.json"),
            "complete_command_resources": {s["name"]: s["resource"] for s in invocation["steps"]},
            "Julia_preparation_observations": report["resources"], "Julia_evaluation_observations": report["evaluation_times"]})
    payload = max(case["Julia_preparation_observations"]["row_coefficient_payload_bytes"] for case in cases)
    return {"status": COMPLETED, "scope": SCOPE, **{k: False for k in FALSE_FLAGS},
        "code_commit": git_commit(ROOT), "source_fingerprint": source_fingerprint(ROOT), "policy_sha256": POLICY_SHA,
        "request_sha256": file_hash(request_path), "source_request_points": len(request["points"]),
        "pilot_cases": cases, "pilot_case_coverage_verified": True,
        "engineering_projection": {"condition": "IF every full132 context had the same payload as the maximum of these three measured cases",
            "sampled_maximum_row_payload_bytes": payload, "conditional_all132_payload_bytes": 132*payload,
            "actual_all132_payload_measured": False, "actual_all132_RSS_measured": False,
            "production_performance_policy_calibrated_or_accepted": False}}


def build_pilot_ledger(request_path, run, output, *, julia="julia", time_binary="/usr/bin/time"):
    output = Path(output)
    if output.exists():
        raise ValueError("preserve earlier pilot ledgers; use a fresh output")
    ledger = _measurement_ledger(request_path, run, julia=julia, time_binary=time_binary)
    atomic_json(output, ledger)
    return ledger


def verify_pilot_ledger(request_path, run, output, *, julia="julia", time_binary="/usr/bin/time"):
    expected = _measurement_ledger(request_path, run, julia=julia, time_binary=time_binary)
    if read_json(output) != expected:
        raise ValueError("changed, incomplete or stale three-case pilot ledger")
    return {"status": "CACHED_THREE_CASE_PILOT_REVIEWED_NOT_DIAGNOSTIC_GATE", "scope": SCOPE,
            **{k: False for k in FALSE_FLAGS}}
