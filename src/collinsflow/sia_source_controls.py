"""Schedule complete synthetic Julia control jobs; no physical SIA dispatch."""
from __future__ import annotations

from datetime import datetime, timezone
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
import time
import uuid

ROOT = Path(__file__).resolve().parents[2]
SCOPE = "explicitly_synthetic_sia_source_and_stencil_controls"
COMPLETED = "COMPLETE_SYNTHETIC_SIA_SOURCE_CONTROL_PAIR_NOT_PHYSICAL_GATE"
COMPUTED = "SYNTHETIC_SIA_SOURCE_CONTROLS_COMPUTED_NOT_DIAGNOSTIC_GATE"
REVIEWED = "SYNTHETIC_SIA_SOURCE_CONTROL_ARITHMETIC_REVIEWED_NOT_DIAGNOSTIC_GATE"
FLAGS = ("production_ready", "physical_extraction", "historical_fit_ready", "experimental_data_admitted",
         "fit_performed", "G2PERF_pass", "interval_coverage_verified", "accepted_fitting_gate",
         "accepted_prediction_diagnostic", "numerical_dispatch_ready")
FORMAT = "COLLINSFLOW_TIME_V1\t%M\t%e\t%U\t%S\t%x"
TIME = Path("/usr/bin/time")


def atomic_json(path, value):
    path = Path(path)
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def planner():
    spec = importlib.util.spec_from_file_location("sia_source_plan", ROOT / "src/collinsflow/sia_geometry_jobs.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def resource_record(path):
    rows = [line.split("\t") for line in Path(path).read_text(encoding="utf-8").splitlines()
            if line.startswith("COLLINSFLOW_TIME_V1\t")]
    if len(rows) != 1 or len(rows[0]) != 6:
        raise ValueError("one complete GNU-time resource record required")
    _, rss, wall, user, system, code = rows[0]
    rss, code = int(rss), int(code)
    wall, user, system = float(wall), float(user), float(system)
    if rss <= 0 or not 0 <= code <= 255 or any(not math.isfinite(v) or v < 0 for v in (wall, user, system)):
        raise ValueError("invalid whole-command resource measurements")
    return dict(maximum_RSS_KiB=rss, wall_seconds=wall, user_seconds=user, system_seconds=system,
                child_exit_code=code, scope="GNU_time_direct_child_maximum_not_simultaneous_process_tree_peak")


def available_memory():
    values = {line.split(":", 1)[0]: line.split(":", 1)[1].split() for line in Path("/proc/meminfo").read_text().splitlines()}
    value, unit = values["MemAvailable"]
    if unit != "kB":
        raise ValueError("unknown available-memory units")
    return int(value) * 1024


def _environment():
    return {**os.environ, "JULIA_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}


def successful_command(step, budget):
    resource = step.get("resource", {})
    elapsed = step.get("wrapper_elapsed_seconds")
    rss = resource.get("maximum_RSS_KiB")
    return (type(step.get("returncode")) is int and step["returncode"] == 0 and
            type(step.get("process_returncode")) is int and step["process_returncode"] == 0 and
            step.get("interruption") is None and step.get("termination_signal") is None and
            type(resource.get("child_exit_code")) is int and resource["child_exit_code"] == 0 and
            type(rss) is int and 0 < rss <= budget["maximum_command_RSS_KiB"] and
            type(elapsed) in (float, int) and math.isfinite(elapsed) and elapsed >= 0)


def _supervisor(jobs, plan):
    jobs._trusted_api(Path(plan["frozen"]))
    module = importlib.import_module("collinsflow.job_control")
    if Path(module.__file__).resolve().parents[2] != Path(plan["frozen"]):
        raise ValueError("wrong frozen process supervisor")
    return module


def _command(plan, path, sha, folder, name, julia, report_sha=None):
    args = [plan["frozen"], str(path), sha]
    if name == "producer":
        program = "run_sia_source_controls.jl"
        args.append(str(folder / "report.json"))
    else:
        program = "audit_sia_source_controls.jl"
        args += [str(folder / "report.json"), report_sha, str(folder / "audit.json")]
    return [str(TIME), "--output", str(folder / (name + ".resources.txt")), "--format", FORMAT, "--",
            julia, "--startup-file=no", "--project=" + str(Path(plan["frozen"]) / "numerics"),
            str(ROOT / "numerics/bin" / program), *args]


def _inputs(jobs, plan, path, julia):
    paths = {Path(path), Path(plan["request"]), Path(plan["source_request"]), TIME, Path(julia)}
    paths.update(Path(p) for p in plan["source_inputs"])
    paths.update(Path(plan["adapter"]) / name for name in plan["adapter_files"])
    paths.update(Path(plan["frozen"]) / name for name in plan["frozen_files"])
    return {str(p): jobs.file_hash(p) for p in sorted(paths)}


def _reports(jobs, plan, sha, folder):
    report, audit = (jobs.read_json(folder / name) for name in ("report.json", "audit.json"))
    for data, status in ((report, COMPUTED), (audit, REVIEWED)):
        if (data.get("status") != status or data.get("scope") != SCOPE or data.get("theory_scope") != "synthetic_closure" or
                data.get("source_controls_success") is not True or any(data.get(key) is not False for key in FLAGS) or
                data.get("physical_SIA_evolution_or_integrals_or_AD_executed") is not False or not data.get("completed_at")):
            raise ValueError("source controls are incomplete or claim physical admission")
    binding = report["source_context"]
    if binding != audit["source_context"] or binding["source_plan_sha256"] != sha or binding["adapter_head"] != plan["adapter_head"]:
        raise ValueError("source-control binding changed")
    expected = dict(source_plan_sha256=sha, source_request_sha256=plan["request_sha256"], frozen_head=jobs.FROZEN_HEAD,
        numerical_source_fingerprint=jobs.FROZEN_FP, adapter_head=plan["adapter_head"], adapter_files=plan["adapter_files"],
        parent_contract_sha256=jobs.PARENT_SHA, job_policy_sha256=jobs.POLICY_SHA, source_inputs=plan["source_inputs"],
        raw_primary_replay="Python source replay attestation; Julia verifies current full source bindings",
        Julia_version="1.12.6", Julia_executable=binding["Julia_executable"],
        Julia_executable_sha256=jobs.file_hash(binding["Julia_executable"]))
    if binding != expected or audit.get("source_control_report_sha256") != jobs.file_hash(folder / "report.json"):
        raise ValueError("changed complete source or runtime provenance")
    if len(report["checks"]) != 433 or audit.get("stored_arithmetic_checks_recomputed") != 433:
        raise ValueError("incomplete synthetic comparison ledger")
    return report, audit


def run_pair(plan_path, expected_sha256, output, *, julia="julia"):
    jobs = planner()
    path, folder = Path(plan_path).resolve(), Path(output).absolute()
    if (jobs.file_hash(path) != expected_sha256 or folder.resolve() != folder or
            not folder.is_relative_to(path.parent) or not path.parent.is_relative_to(ROOT / "runs") or folder.exists()):
        raise ValueError("externally bound plan and fresh canonical control output required")
    julia = shutil.which(julia)
    if not julia or not Path(julia).is_absolute() or not TIME.is_file():
        raise ValueError("explicit Julia launcher and GNU time required")
    # Preserve Juliaup's launcher argv0; only hash the resolved file contents.
    folder.mkdir(parents=True, exist_ok=False)
    record = {"status": "FAIL", "record_role": "initial_sentinel_not_an_executed_failure", "scope": SCOPE,
              "theory_scope": "synthetic_closure", "started_at": datetime.now(timezone.utc).isoformat(),
              "source_plan": str(path), "source_plan_sha256": expected_sha256, "julia_launcher": julia,
              "steps": [], **dict.fromkeys(FLAGS, False)}
    initial = folder / "initial-invocation.json"
    atomic_json(initial, record)
    record["initial_record_sha256"] = jobs.file_hash(initial)
    record["record_role"] = "complete_source_control_pair_invocation"
    try:
        jobs.verify_plan(path, expected_sha256=expected_sha256)
        plan = jobs.read_json(path)
        supervisor = _supervisor(jobs, plan)
        budget = jobs.load_policy()["resources"]["global"]
        inputs = _inputs(jobs, plan, path, julia)
        record.update(source_plan=str(path), source_plan_sha256=expected_sha256, adapter_head=plan["adapter_head"],
                      julia_launcher=julia, input_hashes=inputs, resources=budget)
        atomic_json(folder / "invocation.json", record)
        started = time.monotonic()
        for name in ("producer", "audit"):
            memory = available_memory()
            remaining = budget["whole_Julia_jobs_time_budget_seconds"] - (time.monotonic() - started)
            if memory < budget["minimum_available_memory_bytes"] or remaining <= 0:
                raise ValueError("source-control whole-job resource prerequisite failed")
            report_sha = jobs.file_hash(folder / "report.json") if name == "audit" else None
            command = _command(plan, path, expected_sha256, folder, name, julia, report_sha)
            began = time.monotonic()
            with (folder / (name + ".log")).open("wb") as stream:
                result = supervisor.run_complete_process(command, cwd=ROOT, env=_environment(), stdout=stream, timeout_seconds=remaining)
            step = dict(name=name, command=command, returncode=result.returncode, process_returncode=result.process_returncode,
                interruption=result.interruption, termination_signal=result.termination_signal,
                wrapper_elapsed_seconds=time.monotonic() - began, minimum_available_memory_observed_before_command=memory,
                log_sha256=jobs.file_hash(folder / (name + ".log")))
            record["steps"].append(step)
            atomic_json(folder / "invocation.json", record)
            resource_path = folder / (name + ".resources.txt")
            if resource_path.is_file():
                step.update(resource=resource_record(resource_path), resource_sha256=jobs.file_hash(resource_path))
                atomic_json(folder / "invocation.json", record)
            if not successful_command(step, budget):
                raise ValueError("source-control command failed, was interrupted or exhausted resources")
        if (sum(s["resource"]["wall_seconds"] for s in record["steps"]) > budget["whole_Julia_jobs_time_budget_seconds"] or
                time.monotonic() - started > budget["whole_Julia_jobs_time_budget_seconds"]):
            raise ValueError("complete source-control pair exceeded its shared time budget")
        _reports(jobs, plan, expected_sha256, folder)
        jobs.verify_plan(path, expected_sha256=expected_sha256)
        if _inputs(jobs, plan, path, julia) != inputs:
            raise ValueError("source-control inputs changed during complete pair")
        record["outputs"] = {str(folder / name): jobs.file_hash(folder / name) for name in
            ("report.json", "audit.json", "producer.log", "audit.log", "producer.resources.txt", "audit.resources.txt")}
        record["status"] = COMPLETED
    except BaseException as error:
        record["error"] = type(error).__name__ + ": " + str(error)
        raise
    finally:
        record["completed_at"] = datetime.now(timezone.utc).isoformat()
        atomic_json(folder / "invocation.json", record)
    (folder / "job.completed").write_text(jobs.file_hash(folder / "invocation.json") + "\n", encoding="utf-8")
    return record


def _fresh_cached_review(jobs, plan, path, sha, folder, julia):
    review = path.parent / "source-control-cache-reviews" / uuid.uuid4().hex
    review.mkdir(parents=True, exist_ok=False)
    report_sha = jobs.file_hash(folder / "report.json")
    command = _command(plan, path, sha, folder, "audit", julia, report_sha)
    command[2], command[-1] = str(review / "audit.resources.txt"), str(review / "audit.json")
    budget = jobs.load_policy()["resources"]["global"]
    inputs = {**_inputs(jobs, plan, path, julia),
              **{str(folder / name): jobs.file_hash(folder / name) for name in ("report.json", "invocation.json", "job.completed")}}
    result = {"status": "FAIL", "scope": SCOPE, "theory_scope": "synthetic_closure", "command": command,
              "source_plan_sha256": sha, "source_control_report_sha256": report_sha, "input_hashes": inputs,
              "resources": budget, "record_role": "initial_sentinel_not_an_executed_failure", **dict.fromkeys(FLAGS, False)}
    atomic_json(review / "initial-invocation.json", result)
    result.update(initial_record_sha256=jobs.file_hash(review / "initial-invocation.json"),
                  record_role="cached_source_control_reader_invocation")
    atomic_json(review / "invocation.json", result)
    try:
        memory = available_memory()
        if memory < budget["minimum_available_memory_bytes"]:
            raise ValueError("cached source-control reader memory prerequisite failed")
        supervisor = _supervisor(jobs, plan)
        began = time.monotonic()
        with (review / "audit.log").open("wb") as stream:
            process = supervisor.run_complete_process(command, cwd=ROOT, env=_environment(), stdout=stream,
                timeout_seconds=budget["whole_Julia_jobs_time_budget_seconds"])
        result.update(returncode=process.returncode, process_returncode=process.process_returncode,
            interruption=process.interruption, termination_signal=process.termination_signal,
            wrapper_elapsed_seconds=time.monotonic() - began, minimum_available_memory_observed_before_command=memory,
            log_sha256=jobs.file_hash(review / "audit.log"))
        if (review / "audit.resources.txt").is_file():
            result.update(resource=resource_record(review / "audit.resources.txt"),
                          resource_sha256=jobs.file_hash(review / "audit.resources.txt"))
        if (not successful_command(result, budget) or
                result["resource"]["wall_seconds"] > budget["whole_Julia_jobs_time_budget_seconds"] or
                result["wrapper_elapsed_seconds"] > budget["whole_Julia_jobs_time_budget_seconds"]):
            raise ValueError("fresh cached source-control arithmetic reader failed or exceeded resources")
        audit = jobs.read_json(review / "audit.json")
        if (audit.get("status") != REVIEWED or audit.get("scope") != SCOPE or audit.get("theory_scope") != "synthetic_closure" or
                audit.get("source_controls_success") is not True or any(audit.get(key) is not False for key in FLAGS) or
                audit.get("physical_SIA_evolution_or_integrals_or_AD_executed") is not False or
                audit.get("raw_primary_tables_independently_reparsed_in_Julia") is not False or
                audit.get("Fourier_quadrature_reexecuted_by_reader") is not False or
                audit.get("closed_mathematical_functions_reevaluated") is not True or
                audit.get("source_control_report_sha256") != report_sha or
                audit.get("source_context") != jobs.read_json(folder / "report.json")["source_context"] or
                audit.get("stored_arithmetic_checks_recomputed") != 433 or jobs.file_hash(folder / "report.json") != report_sha or
                any(jobs.file_hash(p) != value for p, value in inputs.items())):
            raise ValueError("fresh cached source-control audit is stale or incomplete")
        result.update(status="FRESH_CACHED_SYNTHETIC_SOURCE_CONTROL_REVIEW_COMPLETE_NOT_PHYSICAL_GATE",
                      audit_sha256=jobs.file_hash(review / "audit.json"))
    except BaseException as error:
        result["error"] = type(error).__name__ + ": " + str(error)
        raise
    finally:
        result["completed_at"] = datetime.now(timezone.utc).isoformat()
        atomic_json(review / "invocation.json", result)
    return dict(path=str(review / "invocation.json"), sha256=jobs.file_hash(review / "invocation.json"), status=result["status"])


def verify_pair(output, plan_path, expected_sha256, *, julia="julia"):
    jobs = planner()
    folder, path = Path(output).absolute(), Path(plan_path).absolute()
    if folder.resolve() != folder or path.resolve() != path or not folder.is_relative_to(path.parent):
        raise ValueError("canonical source-control run and externally bound plan required")
    julia = shutil.which(julia)
    if not julia or not Path(julia).is_absolute():
        raise ValueError("explicit Julia launcher required")
    invocation_sha = jobs.file_hash(folder / "invocation.json")
    record = jobs.read_json(folder / "invocation.json")
    if (record.get("status") != COMPLETED or record.get("scope") != SCOPE or record.get("theory_scope") != "synthetic_closure" or
            record.get("record_role") != "complete_source_control_pair_invocation" or
            not record.get("completed_at") or any(record.get(key) is not False for key in FLAGS) or
            jobs.file_hash(folder / "invocation.json") != invocation_sha or
            (folder / "job.completed").read_text().strip() != invocation_sha):
        raise ValueError("complete source-control provenance required; a marker alone is insufficient")
    sha = expected_sha256
    if record["source_plan"] != str(path) or record["source_plan_sha256"] != sha or record["julia_launcher"] != julia:
        raise ValueError("source-control run differs from external plan or launcher")
    jobs.verify_plan(path, expected_sha256=sha)
    plan = jobs.read_json(path)
    if (record["adapter_head"] != plan["adapter_head"] or record["input_hashes"] != _inputs(jobs, plan, path, julia) or
            record["initial_record_sha256"] != jobs.file_hash(folder / "initial-invocation.json") or
            jobs.read_json(folder / "initial-invocation.json")["record_role"] != "initial_sentinel_not_an_executed_failure"):
        raise ValueError("changed source-control input or initial sentinel")
    products = {str(folder / name) for name in
        ("report.json", "audit.json", "producer.log", "audit.log", "producer.resources.txt", "audit.resources.txt")}
    if set(record["outputs"]) != products or any(jobs.file_hash(p) != sha for p, sha in record["outputs"].items()):
        raise ValueError("changed complete source-control output inventory")
    budget = jobs.load_policy()["resources"]["global"]
    if record["resources"] != budget or [s["name"] for s in record["steps"]] != ["producer", "audit"]:
        raise ValueError("source-control budget or command pair changed")
    for step in record["steps"]:
        name = step["name"]
        expected = _command(plan, path, sha, folder, name, record["julia_launcher"], jobs.file_hash(folder / "report.json"))
        if (step["command"] != expected or not successful_command(step, budget) or
                step["log_sha256"] != jobs.file_hash(folder / (name + ".log")) or
                step["resource_sha256"] != jobs.file_hash(folder / (name + ".resources.txt")) or
                step["resource"] != resource_record(folder / (name + ".resources.txt")) or
                step["minimum_available_memory_observed_before_command"] < budget["minimum_available_memory_bytes"]):
            raise ValueError("stale source-control command or complete resource record")
    if (sum(s["resource"]["wall_seconds"] for s in record["steps"]) > budget["whole_Julia_jobs_time_budget_seconds"] or
            sum(s["wrapper_elapsed_seconds"] for s in record["steps"]) > budget["whole_Julia_jobs_time_budget_seconds"]):
        raise ValueError("source-control whole-job budget exceeded")
    _reports(jobs, plan, sha, folder)
    fresh = _fresh_cached_review(jobs, plan, path, sha, folder, julia)
    jobs.verify_plan(path, expected_sha256=sha)
    if (record["input_hashes"] != _inputs(jobs, plan, path, julia) or
            any(jobs.file_hash(p) != value for p, value in record["outputs"].items()) or
            jobs.file_hash(folder / "invocation.json") != invocation_sha or
            (folder / "job.completed").read_text().strip() != invocation_sha):
        raise ValueError("source-control inputs or outputs changed during cached review")
    return {**record, "fresh_cache_review": fresh}
