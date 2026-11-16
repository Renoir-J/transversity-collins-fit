"""Schedule complete Julia fixed-context storage jobs; never parameter callbacks."""
from __future__ import annotations

import math
from pathlib import Path
import shutil
import subprocess
import time
import uuid

from .historical_sia_points import verify_request as verify_sia
from .historical_sidis_points import verify_request as verify_sidis
from .historical_sidis_pilot import (
    TIME_FORMAT, _pilot_lock, _require_command_success, _timed_job, parse_resource_file,
)
from .provenance import atomic_json, file_hash, git_commit, read_json, source_fingerprint, utc_now

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/fixed-leptonic-storage-contract.json"
CONTRACT_SHA = "cf7de20446d7b80c57559de54d807ff8b73729240cfa0c3f6bbde464ed1c0996"
SCOPE = "reconstructed_fixed_leptonic_storage_roundtrip"
FLAGS = ("production_ready", "physical_extraction", "historical_fit_ready", "experimental_data_admitted",
         "fit_performed", "accepted_prediction_diagnostic", "accepted_fitting_gate", "G2PERF_pass", "interval_coverage_verified")
COMPLETE = "COMPLETE_STORAGE_CASE_RECORDED_NOT_PHYSICS_ADMISSION"
PRODUCED = "COLD_RELOADED_FIXED_CONTEXT_BITS_RECORDED_NOT_PHYSICS_ADMISSION"
AUDITED = "FRESH_LOADED_CONTEXT_BITS_AND_SOURCES_REVIEWED_NOT_PHYSICS_ADMISSION"
TIME = Path("/usr/bin/time")


def storage_policy():
    if file_hash(CONTRACT) != CONTRACT_SHA:
        raise ValueError("changed prospective storage contract")
    policy = read_json(CONTRACT)
    if policy["scope"] != "fixed_leptonic_context_storage_integrity" or any(policy.get(k) is not False for k in FLAGS):
        raise ValueError("storage contract cannot self-admit")
    if any(file_hash(ROOT / p) != h for p, h in policy["parent_sha256"].items()):
        raise ValueError("changed storage parent contract")
    return policy


def _budget(policy):
    r = policy["validation"]["physical_job_resources"]
    return {"resources": {
        "minimum_available_memory_bytes": r["minimum_available_memory_bytes"],
        "maximum_command_RSS_KiB": r["maximum_observed_direct_child_RSS_KiB"],
        "whole_Julia_jobs_time_budget_seconds": r["whole_producer_auditor_time_budget_seconds"],
    }}


def _deadline(policy, selected):
    limit = policy["validation"]["physical_job_resources"]["whole_producer_auditor_time_budget_seconds"]
    if selected is not None:
        if isinstance(selected, bool) or not isinstance(selected, (int, float)) or not math.isfinite(selected) or not 0 < selected <= limit:
            raise ValueError("caller budget may only tighten the frozen storage budget")
        limit = selected
    return time.monotonic() + limit


def _case(policy, index):
    if type(index) is not int or not 1 <= index <= 4:
        raise ValueError("storage case index must be1..4")
    return policy["validation"]["physical_cases"][index - 1]


def _sources(sidis, sia, native, parton):
    requests = (verify_sidis(sidis), verify_sia(sia))
    fingerprint = source_fingerprint(ROOT)
    paths = {str(p): file_hash(p) for p in (sidis, sia, native, parton)}
    for request in requests:
        if request["source_fingerprint"] != fingerprint or Path(request["source_preparation"]["manifest"]).parent.parent != ROOT:
            raise ValueError("current trusted source requests from this checkout required")
        for p, h in request["source_inputs"].items():
            if p in paths and paths[p] != h:
                raise ValueError("conflicting storage source input")
            paths[p] = h
    return fingerprint, paths


def _command(julia, program, *arguments):
    return [julia, "--startup-file=no", "--project=" + str(ROOT / "numerics"),
            str(ROOT / "numerics/bin" / program), *map(str, arguments)]


def _outputs(folder):
    return {str(folder / p) for p in (
        "report.json", "audit.json", "producer.log", "producer.resources.txt", "audit.log", "audit.resources.txt",
        "context/fixed-context-manifest.json", "context/fixed-context-payload.bin",
    )}


def _reports(folder, case_id, fingerprint):
    report, audit = read_json(folder / "report.json"), read_json(folder / "audit.json")
    for data, status in ((report, PRODUCED), (audit, AUDITED)):
        if (data.get("status") != status or data.get("scope") != SCOPE or data.get("case_id") != case_id or
                data.get("source_fingerprint") != fingerprint or data.get("storage_contract_sha256") != CONTRACT_SHA or
                data.get("input_kind") != "pinned_inputs" or any(data.get(k) is not False for k in FLAGS)):
            raise ValueError("failed or stale source-bound storage report")
    if audit["producer_report_sha256"] != file_hash(folder / "report.json"):
        raise ValueError("storage audit does not bind the producer")
    manifest = folder / "context/fixed-context-manifest.json"
    if report["manifest_sha256"] != audit["manifest_sha256"] or report["manifest_sha256"] != file_hash(manifest):
        raise ValueError("changed stored manifest")
    if report["payload_sha256"] != file_hash(folder / "context/fixed-context-payload.bin"):
        raise ValueError("changed stored payload")
    if (audit.get("full13_gradient_and_directional_HVP_reexecuted") is not True or
            audit.get("cold_native_preparation_reexecuted") is not False or
            audit.get("whole_command_resources_or_temporal_reuse_replayed") is not False or
            audit.get("loaded_observables_reexecuted") != (1 if "SIDIS" in case_id else 2)):
        raise ValueError("incomplete fresh loaded-context numerical review")
    return report, audit


def run_case(sidis, sia, native, parton, run, index, *, julia="julia", timeout_seconds=None):
    if type(index) is not int or not 1 <= index <= 4:
        raise ValueError("storage case index must be1..4")
    sidis, sia, native, parton, run = map(lambda p: Path(p).resolve(), (sidis, sia, native, parton, run))
    folder = run / f"case{index:02d}"
    if folder.exists() and any(folder.iterdir()):
        raise ValueError("preserve earlier storage attempts; use fresh paths")
    folder.mkdir(parents=True, exist_ok=True)
    record = {"status": "FAIL", "scope": SCOPE, "case_index": index, "started_at": utc_now(),
              "steps": [], "outputs": {}, **dict.fromkeys(FLAGS, False)}
    initial = folder / "initial-invocation.json"
    atomic_json(initial, {**record, "record_role": "initial_sentinel_not_an_executed_failure"})
    record["initial_record_sha256"] = file_hash(initial)
    atomic_json(folder / "invocation.json", record)
    try:
        policy = storage_policy()
        case_id = _case(policy, index)
        fingerprint, inputs = _sources(sidis, sia, native, parton)
        if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip():
            raise ValueError("freeze committed storage code before physical jobs")
        resolved_julia = shutil.which(julia)
        if resolved_julia is None:
            raise ValueError("selected Julia executable is unavailable")
        executable = Path(resolved_julia).absolute()
        record.update(case_id=case_id, source_fingerprint=fingerprint, code_commit=git_commit(ROOT),
                      contract_sha256=CONTRACT_SHA, input_hashes=inputs, time_binary=str(TIME),
                      julia_executable=str(executable), julia_binary_resolved_path=str(executable.resolve()),
                      julia_executable_sha256=file_hash(executable))
        arguments = (sidis, sia, native, parton, case_id)
        producer = _command(str(executable), "run_fixed_context_storage_case.jl", *arguments, folder / "report.json")
        auditor = _command(str(executable), "audit_fixed_context_storage_case.jl", *arguments, folder / "report.json", folder / "audit.json")
        budget = _budget(policy)
        with _pilot_lock(ROOT):
            deadline = _deadline(policy, timeout_seconds)
            for name, command in (("producer", producer), ("audit", auditor)):
                row = _timed_job(name, command, ROOT, folder, TIME, budget, deadline)
                record["steps"].append(row)
                atomic_json(folder / "invocation.json", record)
                _require_command_success(row, budget)
        _reports(folder, case_id, fingerprint)
        if source_fingerprint(ROOT) != fingerprint or any(file_hash(p) != h for p, h in inputs.items()):
            raise ValueError("storage source or inputs changed during complete case")
        record["outputs"] = {p: file_hash(p) for p in _outputs(folder)}
        record["status"] = COMPLETE
    except BaseException as error:
        record["error"] = type(error).__name__ + ": " + str(error)
        raise
    finally:
        record["completed_at"] = utc_now()
        atomic_json(folder / "invocation.json", record)
    return {k: record[k] for k in ("status", "case_index", *FLAGS)}


def verify_metadata(sidis, sia, native, parton, run, index):
    sidis, sia, native, parton, run = map(lambda p: Path(p).resolve(), (sidis, sia, native, parton, run))
    policy = storage_policy()
    case_id = _case(policy, index)
    fingerprint, inputs = _sources(sidis, sia, native, parton)
    folder = Path(run).resolve() / f"case{index:02d}"
    inv = read_json(folder / "invocation.json")
    if (inv.get("status") != COMPLETE or type(inv.get("case_index")) is not int or inv.get("case_index") != index or inv.get("case_id") != case_id or
            inv.get("source_fingerprint") != fingerprint or inv.get("contract_sha256") != CONTRACT_SHA or
            inv.get("code_commit") != git_commit(ROOT) or not inv.get("completed_at") or
            any(inv.get(k) is not False for k in FLAGS)):
        raise ValueError("incomplete or stale storage invocation")
    if inv["input_hashes"] != inputs or any(file_hash(p) != h for p, h in inputs.items()):
        raise ValueError("changed storage input binding")
    if set(inv["outputs"]) != _outputs(folder) or any(file_hash(p) != h for p, h in inv["outputs"].items()):
        raise ValueError("changed storage case output")
    if file_hash(folder / "initial-invocation.json") != inv["initial_record_sha256"]:
        raise ValueError("changed storage initial record")
    if read_json(folder / "initial-invocation.json")["record_role"] != "initial_sentinel_not_an_executed_failure":
        raise ValueError("wrong initial storage sentinel")
    if [s["name"] for s in inv["steps"]] != ["producer", "audit"] or inv["time_binary"] != str(TIME):
        raise ValueError("missing complete storage producer/auditor pair")
    executable = Path(inv["julia_executable"])
    if (not executable.is_file() or str(executable.resolve()) != inv["julia_binary_resolved_path"] or
            file_hash(executable) != inv["julia_executable_sha256"]):
        raise ValueError("changed storage Julia executable")
    args = (sidis, sia, native, parton, case_id)
    expected_commands = [
        _command(str(executable), "run_fixed_context_storage_case.jl", *args, folder / "report.json"),
        _command(str(executable), "audit_fixed_context_storage_case.jl", *args, folder / "report.json", folder / "audit.json"),
    ]
    for step, program in zip(inv["steps"], expected_commands):
        expected = [str(TIME), "--output", str(folder / (step["name"] + ".resources.txt")),
                    "--format", TIME_FORMAT, "--", *program]
        if step["command"] != expected:
            raise ValueError("changed complete storage command")
        if (file_hash(folder / (step["name"] + ".log")) != step["log_sha256"] or
                file_hash(folder / (step["name"] + ".resources.txt")) != step["resource_sha256"] or
                parse_resource_file(folder / (step["name"] + ".resources.txt")) != step["resource"]):
            raise ValueError("changed raw storage command resources")
        _require_command_success(step, _budget(policy))
    if sum(s["wrapper_elapsed_seconds"] for s in inv["steps"]) > _budget(policy)["resources"]["whole_Julia_jobs_time_budget_seconds"]:
        raise ValueError("storage pair exceeded its frozen whole-job budget")
    report, audit = _reports(folder, case_id, fingerprint)
    return policy, folder, inv, report, audit


def verify_case(sidis, sia, native, parton, run, index, *, julia="julia", timeout_seconds=None):
    sidis, sia, native, parton, run = map(lambda p: Path(p).resolve(), (sidis, sia, native, parton, run))
    policy, folder, inv, report, _ = verify_metadata(sidis, sia, native, parton, run, index)
    review = run / "cache-reviews" / (f"case{index:02d}-" + uuid.uuid4().hex)
    review.mkdir(parents=True, exist_ok=False)
    record = {"status": "FAIL", "scope": "fresh_loaded_storage_case_review_not_native_reexecution",
              "started_at": utc_now(), "original_invocation_sha256": file_hash(folder / "invocation.json"),
              **dict.fromkeys(FLAGS, False)}
    try:
        command = _command(julia, "audit_fixed_context_storage_case.jl", sidis, sia, native, parton,
                           inv["case_id"], folder / "report.json", review / "audit.json")
        with _pilot_lock(ROOT):
            row = _timed_job("cache-audit", command, ROOT, review, TIME, _budget(policy), _deadline(policy, timeout_seconds))
        record["step"] = row
        _require_command_success(row, _budget(policy))
        audit = read_json(review / "audit.json")
        if (audit["status"] != AUDITED or audit["producer_report_sha256"] != file_hash(folder / "report.json") or
                audit["manifest_sha256"] != report["manifest_sha256"] or audit.get("case_id") != inv["case_id"] or
                audit.get("source_fingerprint") != inv["source_fingerprint"] or audit.get("storage_contract_sha256") != CONTRACT_SHA or
                audit.get("input_kind") != "pinned_inputs" or audit.get("full13_gradient_and_directional_HVP_reexecuted") is not True or
                audit.get("cold_native_preparation_reexecuted") is not False or
                any(audit.get(k) is not False for k in FLAGS)):
            raise ValueError("fresh storage cache auditor rejected the case")
        verify_metadata(sidis, sia, native, parton, run, index)
        record["status"] = "CURRENT_LOADED_STORAGE_CASE_VERIFIED_NOT_PHYSICS_ADMISSION"
        record["audit_sha256"] = file_hash(review / "audit.json")
    except BaseException as error:
        record["error"] = type(error).__name__ + ": " + str(error)
        raise
    finally:
        record["completed_at"] = utc_now()
        atomic_json(review / "invocation.json", record)
    return {"case_id": inv["case_id"], "original_invocation_sha256": file_hash(folder / "invocation.json"),
            "manifest_sha256": report["manifest_sha256"], "payload_sha256": report["payload_sha256"],
            "logical_context_sha256": report["logical_context_sha256"], "row_alias_sha256": report["row_alias_sha256"]}


def build_or_verify_ledger(sidis, sia, native, parton, run, output, *, julia="julia", timeout_seconds=None):
    cases = [verify_case(sidis, sia, native, parton, run, i, julia=julia, timeout_seconds=timeout_seconds) for i in range(1, 5)]
    core = {"scope": SCOPE, "source_fingerprint": source_fingerprint(ROOT), "contract_sha256": CONTRACT_SHA,
            "cases": cases, "native_preparation_reexecuted_by_cache_reader": False, **dict.fromkeys(FLAGS, False)}
    output = Path(output)
    if output.exists():
        ledger = read_json(output)
        if any(ledger.get(k) != v for k, v in core.items()) or ledger.get("status") != "FOUR_STORAGE_CASES_RECORDED_NOT_PHYSICS_ADMISSION":
            raise ValueError("changed current storage ledger")
    else:
        atomic_json(output, {"status": "FOUR_STORAGE_CASES_RECORDED_NOT_PHYSICS_ADMISSION", "created_at": utc_now(), **core})
    return {"status": "CURRENT_FOUR_STORAGE_CASE_LEDGER_VERIFIED_NOT_PHYSICS_ADMISSION",
            "ledger_sha256": file_hash(output), **dict.fromkeys(FLAGS, False)}
