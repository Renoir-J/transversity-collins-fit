"""Schedule complete SIA Julia jobs; all selection and numerical checks stay in Julia."""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import importlib.util
import json
import os
from pathlib import Path
import shutil
import time
import uuid


def sibling(name):
    spec = importlib.util.spec_from_file_location(name + "_execution_adapter", Path(__file__).with_name(name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


jobs, controls, host = (sibling(name) for name in ("sia_geometry_jobs", "sia_source_controls", "sia_geometry_host"))
ROOT = jobs.ROOT
FLAGS = controls.FLAGS
PREPARED = "PREPARED_COMPLETE_SIA_JOB_RUN_NOT_PHYSICS_OR_FITTING_ADMISSION"
COMPLETE = "COMPLETE_SIA_JULIA_JOB_PAIR_NOT_FITTING_GATE"
PAIR_SCOPE = "source_bound_complete_sia_Julia_job_pair"
GLOBALS = ("native-select3", "native-selection", "tail-selection", "diagnostic")
PHYSICAL = ("native123", "native4", "tailref", "derivative")
PRODUCTS = ("report.json", "audit.json", "producer.log", "audit.log", "producer.resources.txt", "audit.resources.txt")
file_hash, read_json, atomic_json = jobs.file_hash, jobs.read_json, controls.atomic_json


def now():
    return datetime.now(timezone.utc).isoformat()


def canonical(path):
    path = Path(path).absolute()
    if path.resolve() != path:
        raise ValueError("canonical non-aliased run paths required")
    return path


def check_hashes(inventory):
    if not isinstance(inventory, dict) or not inventory or any(file_hash(path) != sha for path, sha in inventory.items()):
        raise ValueError("missing or changed complete input/output inventory")


def build_inventory(native, parton):
    result = {}
    for path in (native, parton):
        receipt = read_json(path)
        for key in ("source_hashes", "output_hashes"):
            values = receipt[key]
            if not values or any(name in result and result[name] != sha for name, sha in values.items()):
                raise ValueError("inconsistent native/parton source or binary inventory")
            result.update(values)
        result[str(path)] = file_hash(path)
    check_hashes(result)
    return result


def base_inputs(manifest, path):
    return {**manifest["execution_inputs"], str(path): file_hash(path)}


def execution_inputs(manifest):
    plan = read_json(manifest["source_plan"])
    result = controls._inputs(jobs, plan, Path(manifest["source_plan"]), manifest["julia_launcher"])
    result.update(manifest["build_inputs"])
    result.update(manifest["retained_files"])
    result.update(manifest["source_control_files"])
    result[manifest["sidis_owner_file"]] = manifest["sidis_owner_file_sha256"]
    return result


def prepare(plan_path, plan_sha, native, parton, source_controls, sidis_owner, sidis_roots,
            retained_inventory, output, *, julia="julia"):
    """Prepare a new immutable run; preparation never launches physical jobs."""
    plan_path, output = canonical(plan_path), canonical(output)
    if output.exists() or output.parent != plan_path.parent or output.name != "execution.json":
        raise ValueError("use a fresh execution.json beside the current source plan")
    jobs.verify_plan(plan_path, expected_sha256=plan_sha)
    plan = read_json(plan_path)
    julia = shutil.which(julia)
    if not julia or not Path(julia).is_absolute():
        raise ValueError("explicit Julia launcher required; preserve its argv0")
    pair = controls.verify_pair(source_controls, plan_path, plan_sha, julia=julia)
    owner_path = canonical(sidis_owner)
    owner = read_json(owner_path)
    host._owner(owner)
    roots = [host._run_root(str(canonical(root))) for root in sidis_roots]
    if not roots or not any(owner_path.is_relative_to(Path(root)) for root in roots):
        raise ValueError("bind the actual current SIDIS owner file and its run root")
    retained = read_json(retained_inventory)
    # This is an externally reviewed inventory, never a source of executable code.
    if retained.get("scope") != "preserved_original_SIA_files_not_reusable_jobs" or len(retained.get("files", {})) != 9:
        raise ValueError("retain the reviewed original9 SIA file inventory explicitly")
    check_hashes(retained["files"])
    folder = canonical(source_controls)
    control_files = {str(folder / name): file_hash(folder / name) for name in (*PRODUCTS, "invocation.json", "initial-invocation.json", "job.completed")}
    fresh = pair["fresh_cache_review"]
    # Include the complete fresh-review record, logs and resources, not only its status.
    review_folder = canonical(fresh["path"]).parent
    for name in ("invocation.json", "initial-invocation.json", "audit.json", "audit.log", "audit.resources.txt"):
        control_files[str(review_folder / name)] = file_hash(review_folder / name)
    control_files[str(canonical(retained_inventory))] = file_hash(retained_inventory)
    manifest = dict(schema_version=1, status=PREPARED, scope=jobs.SCOPE, theory_scope="reconstructed_theory",
        **dict.fromkeys(FLAGS, False), source_plan=str(plan_path), source_plan_sha256=plan_sha,
        run=str(output.parent), manifest=str(output), adapter_head=plan["adapter_head"],
        native=str(canonical(native)), parton=str(canonical(parton)),
        build_inputs=build_inventory(canonical(native), canonical(parton)),
        source_controls=str(folder), source_control_files=control_files,
        source_control_fresh_review=fresh, sidis_owner_file=str(owner_path), sidis_owner_file_sha256=file_hash(owner_path),
        sidis_owner=owner, sidis_run_roots=roots, retained_files=retained["files"],
        julia_launcher=julia, resources=jobs.load_policy()["resources"], prepared_at=now(),
        physical_dispatch_started=False, limitations=plan["limitations"])
    manifest["execution_inputs"] = execution_inputs(manifest)
    atomic_json(output, manifest)
    return manifest


def load(path, sha, *, replay=False):
    path = canonical(path)
    if file_hash(path) != sha:
        raise ValueError("changed externally bound execution manifest")
    m = read_json(path)
    if (m.get("schema_version") != 1 or type(m.get("schema_version")) is not int or m.get("status") != PREPARED or
            m.get("scope") != jobs.SCOPE or m.get("theory_scope") != "reconstructed_theory" or
            any(m.get(key) is not False for key in (*FLAGS, "physical_dispatch_started")) or
            m["manifest"] != str(path) or m["run"] != str(path.parent) or path.name != "execution.json"):
        raise ValueError("wrong complete-job manifest or attempted admission")
    plan = read_json(m["source_plan"])
    if (canonical(m["source_plan"]).parent != path.parent or not path.is_relative_to(ROOT / "runs") or
            file_hash(m["source_plan"]) != m["source_plan_sha256"] or m["adapter_head"] != plan["adapter_head"] or
            m["resources"] != jobs.load_policy()["resources"] or m["limitations"] != plan["limitations"]):
        raise ValueError("changed run, source plan, resources or limitations")
    jobs.verify_checkout(ROOT, head=plan["adapter_head"], ancestor=jobs.POLICY_HEAD)
    jobs.verify_checkout(Path(plan["frozen"]), head=jobs.FROZEN_HEAD)
    if replay:
        jobs.verify_plan(m["source_plan"], expected_sha256=m["source_plan_sha256"])
    if build_inventory(Path(m["native"]), Path(m["parton"])) != m["build_inputs"]:
        raise ValueError("changed build receipt or binaries")
    host._owner(m["sidis_owner"])
    if read_json(m["sidis_owner_file"]) != m["sidis_owner"]:
        raise ValueError("SIDIS ownership changed; require a new reviewed run")
    if execution_inputs(m) != m["execution_inputs"]:
        raise ValueError("changed complete execution inputs, launcher or GNU-time executable")
    check_hashes(read_json(Path(m["source_controls"]) / "invocation.json")["input_hashes"])
    check_hashes(base_inputs(m, path))
    return m


def location(m, phase, index=0):
    if phase in GLOBALS and type(index) is int and index == 0:
        path = Path(m["run"]) / "global" / phase
    elif phase in PHYSICAL and type(index) is int and index in (range(1, 37) if phase != "derivative" else (1, 15, 36)):
        path = Path(m["run"]) / phase / f"g{index:04d}"
    else:
        raise ValueError("undeclared complete job phase/index")
    return canonical(path)


def needs_native4(m):
    report = read_json(location(m, "native-select3") / "report.json")
    return report["selection"]["status"] == "NATIVE4_REQUIRED_NOT_DIAGNOSTIC_GATE"


def dependencies(m, phase):
    first = [("native123", i) for i in range(1, 37)]
    if phase == "native123":
        return []
    if phase == "native-select3":
        return first
    if phase == "native4":
        return [("native-select3", 0)]
    native = first + [("native-select3", 0)]
    if needs_native4(m):
        native += [("native4", i) for i in range(1, 37)]
    if phase == "native-selection":
        return native
    native += [("native-selection", 0)]
    if phase == "tailref":
        return native
    tails = native + [("tailref", i) for i in range(1, 37)]
    if phase == "tail-selection":
        return tails
    tails += [("tail-selection", 0)]
    if phase == "derivative":
        return tails
    if phase == "diagnostic":
        return tails + [("derivative", i) for i in (1, 15, 36)]
    raise ValueError("undeclared dependency phase")


def command(m, path, sha, phase, index, folder, name, report_sha=None):
    arguments = [str(path), sha, phase, str(index)]
    if name == "producer":
        program = "run_sia_geometry_job.jl"
    else:
        program = "audit_sia_geometry_job.jl"
        arguments += [str(location(m, phase, index) / "report.json"), report_sha]
    arguments += [str(folder / ("report.json" if name == "producer" else "audit.json"))]
    return [str(controls.TIME), "--output", str(folder / (name + ".resources.txt")), "--format", controls.FORMAT, "--",
        m["julia_launcher"], "--startup-file=no", "--project=" + str(Path(read_json(m["source_plan"])["frozen"]) / "numerics"),
        str(ROOT / "numerics/bin" / program), *arguments]


def report_context(m, sha):
    binding = read_json(Path(m["source_controls"]) / "report.json")["source_context"]
    return {**binding, "execution_manifest_sha256": sha,
        "native_receipt_sha256": file_hash(m["native"]), "parton_receipt_sha256": file_hash(m["parton"])}


def verify_report_header(m, sha, phase, index, data, status):
    if (type(data.get("schema_version")) is not int or data["schema_version"] != 1 or data.get("status") != status or
            data.get("scope") != jobs.SCOPE or data.get("theory_scope") != "reconstructed_theory" or
            data.get("manifest_sha256") != sha or data.get("phase") != phase or type(data.get("index")) is not int or
            data["index"] != index or not isinstance(data.get("completed_at"), str) or not data["completed_at"] or
            data.get("source_context") != report_context(m, sha) or any(data.get(key) is not False for key in FLAGS)):
        raise ValueError("incomplete or stale complete report identity and source context")


def verify_pair(m, path, sha, phase, index=0, *, require_marker=True, base=None):
    """Metadata review only; a fresh Julia reader still recomputes saved arithmetic."""
    folder = location(m, phase, index)
    invpath = folder / "invocation.json"
    invsha, record = file_hash(invpath), read_json(invpath)
    expected_files = {str(folder / name) for name in PRODUCTS}
    if (record.get("status") != COMPLETE or record.get("scope") != PAIR_SCOPE or
            not isinstance(record.get("completed_at"), str) or not record["completed_at"] or
            record.get("theory_scope") != "reconstructed_theory" or record.get("record_role") != "complete_job_execution_invocation" or
            record.get("manifest_sha256") != sha or record.get("phase") != phase or record.get("index") != index or
            type(record.get("index")) is not int or any(record.get(key) is not False for key in FLAGS) or
            (require_marker and (folder / "job.completed").read_text().strip() != invsha) or
            set(record.get("outputs", {})) != expected_files or record["inputs"] != job_inputs(m, path, sha, phase, index, base=base)):
        raise ValueError("incomplete, substituted or stale complete job ledger")
    check_hashes(record["inputs"])
    check_hashes(record["outputs"])
    initial = read_json(folder / "initial-invocation.json")
    budget = m["resources"][phase if phase in PHYSICAL else "global"]
    if (record["initial_record_sha256"] != file_hash(folder / "initial-invocation.json") or initial["status"] != "FAIL" or
            initial["record_role"] != "initial_sentinel_not_an_executed_failure" or initial["steps"] != [] or
            initial["completed_at"] is not None or initial["manifest_sha256"] != sha or initial.get("scope") != PAIR_SCOPE or
            initial.get("theory_scope") != "reconstructed_theory" or initial.get("phase") != phase or
            type(initial.get("index")) is not int or initial["index"] != index or initial.get("resources") != budget or
            not isinstance(initial.get("started_at"), str) or not initial["started_at"] or
            any(initial.get(key) is not False for key in FLAGS)):
        raise ValueError("changed initial nonterminal sentinel")
    if record["resources"] != budget or [s["name"] for s in record["steps"]] != ["producer", "audit"]:
        raise ValueError("changed whole-job budget or command pair")
    for step in record["steps"]:
        name = step["name"]
        if (step["command"] != command(m, path, sha, phase, index, folder, name, file_hash(folder / "report.json")) or
                not controls.successful_command(step, budget) or
                step["resource"] != controls.resource_record(folder / (name + ".resources.txt")) or
                step["log_sha256"] != file_hash(folder / (name + ".log")) or
                step["resource_sha256"] != file_hash(folder / (name + ".resources.txt")) or
                step["minimum_available_memory_observed_before_command"] < budget["minimum_available_memory_bytes"]):
            raise ValueError("changed measured command, raw/effective exit, log or resource operands")
    if any(sum(s[key] if key == "wrapper_elapsed_seconds" else s["resource"]["wall_seconds"] for s in record["steps"]) >
           budget["whole_Julia_jobs_time_budget_seconds"] for key in ("wrapper_elapsed_seconds", "wall_seconds")):
        raise ValueError("whole producer/auditor budget exceeded")
    producer, audit = (read_json(folder / name) for name in ("report.json", "audit.json"))
    verify_report_header(m, sha, phase, index, producer, "COMPLETE_SIA_JOB_VALUES_COMPUTED_NOT_GATE")
    verify_report_header(m, sha, phase, index, audit, "COMPLETE_SIA_SAVED_JOB_ARITHMETIC_REVIEWED_NOT_GATE")
    if (audit.get("report_sha256") != file_hash(folder / "report.json") or
            type(audit.get("complete_prerequisite_pair_ledgers_verified")) is not int or
            audit["complete_prerequisite_pair_ledgers_verified"] != len(dependencies(m, phase)) or
            any(audit.get(key) is not False for key in ("native_evolution_reexecuted_by_reader",
                "cold_AD_or_integrals_reexecuted_by_reader", "raw_primary_tables_replayed_in_Julia")) or
            (phase in GLOBALS and audit.get("complete72_saved_arithmetic_verified") is not (phase == "diagnostic"))):
        raise ValueError("missing independent complete saved-record audit")
    return record


def job_inputs(m, path, sha, phase, index, *, base=None):
    result = dict(base if base is not None else base_inputs(m, path))
    for dep_phase, dep_index in dependencies(m, phase):
        folder = location(m, dep_phase, dep_index)
        for name in (*PRODUCTS, "initial-invocation.json", "invocation.json", "job.completed", "lease.json"):
            result[str(folder / name)] = file_hash(folder / name)
    own = location(m, phase, index) / "lease.json"
    result[str(own)] = file_hash(own)
    return result


@contextmanager
def driver_lease():
    path = Path("/tmp") / f"collinsflow-sia-family-{os.getuid()}.lock"
    with host._lock_at(path) as lease:
        yield lease


def require_clear(m):
    result = host.observe_host(m["sidis_owner"], m["sidis_run_roots"])
    if result["status"] != "HOST_RESOURCE_GUARD_CLEAR_NOT_SCIENTIFIC_ADMISSION" or not result["live_host_observed"]:
        raise ValueError("current SIDIS owner, orphan or unknown host state blocks SIA dispatch")
    return result


def measured(m, argv, folder, name, budget, remaining):
    available = controls.available_memory()
    if available < budget["minimum_available_memory_bytes"] or remaining <= 0:
        raise ValueError("insufficient declared available memory or whole-job budget")
    supervisor = controls._supervisor(jobs, read_json(m["source_plan"]))
    started = time.monotonic()
    with (folder / (name + ".log")).open("xb") as stream:
        result = supervisor.run_complete_process(argv, cwd=ROOT, env=controls._environment(), stdout=stream, timeout_seconds=remaining)
    step = dict(name=name, command=argv, **result._asdict(), wrapper_elapsed_seconds=time.monotonic() - started,
        minimum_available_memory_observed_before_command=available,
        log_sha256=file_hash(folder / (name + ".log")))
    try:
        step.update(resource_sha256=file_hash(folder / (name + ".resources.txt")),
                    resource=controls.resource_record(folder / (name + ".resources.txt")))
    except (OSError, ValueError) as error:
        step["resource_error"] = type(error).__name__ + ": " + str(error)
    return step


def run_pair(m, path, sha, phase, index, family):
    folder = location(m, phase, index)
    if folder.exists():
        return verify_pair(m, path, sha, phase, index)
    load(path, sha)
    base = base_inputs(m, path)
    for dep in dependencies(m, phase):
        verify_pair(m, path, sha, *dep, base=base)
    with host.geometry_lock() as geometry:
        observation = require_clear(m)
        folder.mkdir(parents=True, exist_ok=False)
        atomic_json(folder / "lease.json", dict(family=family, geometry=geometry, host_observation=observation))
        budget = m["resources"][phase if phase in PHYSICAL else "global"]
        initial = dict(status="FAIL", scope=PAIR_SCOPE, theory_scope="reconstructed_theory", **dict.fromkeys(FLAGS, False),
            record_role="initial_sentinel_not_an_executed_failure", manifest_sha256=sha, phase=phase, index=index,
            started_at=now(), completed_at=None, steps=[], resources=budget)
        atomic_json(folder / "initial-invocation.json", initial)
        record = {**initial, "record_role": "complete_job_execution_invocation", "inputs": job_inputs(m, path, sha, phase, index),
            "initial_record_sha256": file_hash(folder / "initial-invocation.json")}
        atomic_json(folder / "invocation.json", record)
        try:
            for name in ("producer", "audit"):
                load(path, sha)
                require_clear(m)
                report_sha = file_hash(folder / "report.json") if name == "audit" else None
                argv = command(m, path, sha, phase, index, folder, name, report_sha)
                remaining = budget["whole_Julia_jobs_time_budget_seconds"] - sum(s["wrapper_elapsed_seconds"] for s in record["steps"])
                record["steps"].append(dict(name=name, command=argv, state="COMMAND_PREPARED_NOT_COMPLETED"))
                atomic_json(folder / "invocation.json", record)
                step = measured(m, argv, folder, name, budget, remaining)
                record["steps"][-1] = step
                atomic_json(folder / "invocation.json", record)
                if not controls.successful_command(step, budget):
                    raise ValueError("complete Julia command failed; retain all partial outputs")
            load(path, sha)
            check_hashes(record["inputs"])
            record["outputs"] = {str(folder / name): file_hash(folder / name) for name in PRODUCTS}
            record.update(status=COMPLETE, completed_at=now())
            atomic_json(folder / "invocation.json", record)
            verify_pair(m, path, sha, phase, index, require_marker=False, base=base)
            with (folder / "job.completed").open("x", encoding="utf-8") as stream:
                stream.write(file_hash(folder / "invocation.json") + "\n")
            verify_pair(m, path, sha, phase, index)
            return record
        except BaseException as error:
            record.update(status="FAIL", completed_at=now(), error=type(error).__name__ + ": " + str(error))
            atomic_json(folder / "invocation.json", record)
            raise


def run_all(path, sha):
    path = canonical(path)
    m = load(path, sha, replay=True)
    attempt = Path(m["run"]) / "driver" / uuid.uuid4().hex
    attempt.mkdir(parents=True, exist_ok=False)
    record = dict(status="FAIL", scope="complete_sia_run_execution_not_fitting", **dict.fromkeys(FLAGS, False),
        manifest_sha256=sha, started_at=now(), completed_at=None, completed_pairs=[], dispatch_attempted=False)
    atomic_json(attempt / "initial-invocation.json", {**record, "record_role": "initial_sentinel_not_an_executed_failure"})
    atomic_json(attempt / "invocation.json", record)
    try:
        with driver_lease() as family:
            owner = {key: family[key] for key in ("pid", "start_ticks", "boot_id")}
            owner["command"] = Path("/proc/self/cmdline").read_bytes().decode().rstrip("\0").split("\0")
            atomic_json(attempt / "owner.json", owner)
            record["host_observation"] = require_clear(m)
            schedule = [("native123", i) for i in range(1, 37)] + [("native-select3", 0)]
            def complete(phase, index):
                record.update(status="IN_PROGRESS_COMPLETE_SIA_JOBS_NOT_GATE", current_phase=phase, current_index=index)
                record["dispatch_attempted"] = True
                atomic_json(attempt / "invocation.json", record)
                run_pair(m, path, sha, phase, index, family)
                record["completed_pairs"].append(dict(phase=phase, index=index,
                    invocation_sha256=file_hash(location(m, phase, index) / "invocation.json")))
                atomic_json(attempt / "invocation.json", record)
            for pair in schedule:
                complete(*pair)
            if needs_native4(m):
                for i in range(1, 37):
                    complete("native4", i)
            complete("native-selection", 0)
            for i in range(1, 37):
                complete("tailref", i)
            complete("tail-selection", 0)
            for i in (1, 15, 36):
                complete("derivative", i)
            complete("diagnostic", 0)
        result = cached_products(path, sha)
        record.update(status="COMPLETE_SIA_CONDITIONAL_DIAGNOSTIC_NOT_FITTING", completed_at=now(),
            product_receipt=result)
        atomic_json(attempt / "invocation.json", record)
        return result
    except BaseException as error:
        record.update(status="FAIL", completed_at=now(), error=type(error).__name__ + ": " + str(error))
        atomic_json(attempt / "invocation.json", record)
        raise


def cached_products(path, sha):
    """Always run a fresh Julia complete-ledger reader before rendering products."""
    path = canonical(path)
    m = load(path, sha, replay=True)
    original = location(m, "diagnostic")
    verify_pair(m, path, sha, "diagnostic")
    folder = Path(m["run"]) / "cached-reviews" / uuid.uuid4().hex
    folder.mkdir(parents=True, exist_ok=False)
    inputs = job_inputs(m, path, sha, "diagnostic", 0)
    inputs.update({str(original / name): file_hash(original / name) for name in (*PRODUCTS, "invocation.json", "job.completed")})
    record = dict(status="FAIL", scope="fresh_complete72_saved_arithmetic_review_not_fitting", theory_scope="reconstructed_theory",
        **dict.fromkeys(FLAGS, False), manifest_sha256=sha, started_at=now(), completed_at=None,
        inputs=inputs, steps=[], record_role="initial_sentinel_not_an_executed_failure")
    atomic_json(folder / "initial-invocation.json", record)
    record.update(record_role="fresh_cached_complete_ledger_reader", initial_record_sha256=file_hash(folder / "initial-invocation.json"))
    atomic_json(folder / "invocation.json", record)
    try:
        budget = m["resources"]["global"]
        argv = command(m, path, sha, "diagnostic", 0, folder, "audit", file_hash(original / "report.json"))
        record["steps"] = [dict(name="audit", command=argv, state="COMMAND_PREPARED_NOT_COMPLETED")]
        atomic_json(folder / "invocation.json", record)
        step = measured(m, argv, folder, "audit", budget, budget["whole_Julia_jobs_time_budget_seconds"])
        record["steps"] = [step]
        atomic_json(folder / "invocation.json", record)
        if (not controls.successful_command(step, budget) or step["wrapper_elapsed_seconds"] > budget["whole_Julia_jobs_time_budget_seconds"] or
                step["resource"]["wall_seconds"] > budget["whole_Julia_jobs_time_budget_seconds"]):
            raise ValueError("fresh complete-ledger reader failed or exceeded declared resources")
        audit = read_json(folder / "audit.json")
        if (audit.get("status") != "COMPLETE_SIA_SAVED_JOB_ARITHMETIC_REVIEWED_NOT_GATE" or
                audit.get("complete72_saved_arithmetic_verified") is not True or audit.get("manifest_sha256") != sha or
                audit.get("report_sha256") != file_hash(original / "report.json") or any(audit.get(key) is not False for key in FLAGS)):
            raise ValueError("fresh complete72 independent saved-arithmetic review is absent or stale")
        load(path, sha)
        check_hashes(inputs)
        record.update(status="FRESH_COMPLETE72_CONDITIONAL_SAVED_ARITHMETIC_REVIEW_NOT_FITTING", completed_at=now(),
            audit_sha256=file_hash(folder / "audit.json"))
        atomic_json(folder / "invocation.json", record)
        products = render_products(m, read_json(original / "report.json"), folder / "products")
        check_hashes(inputs)
        receipt = dict(status="COMPLETE72_CONDITIONAL_DIAGNOSTIC_PRODUCTS_NOT_FITTING", theory_scope="reconstructed_theory",
            **dict.fromkeys(FLAGS, False), manifest_sha256=sha,
            fresh_cached_invocation=str(folder / "invocation.json"), fresh_cached_invocation_sha256=file_hash(folder / "invocation.json"),
            complete72_saved_arithmetic_verified=True, cold_evolution_or_AD_reexecuted_by_reader=False,
            outputs={str(p): file_hash(p) for p in products}, limitations=m["limitations"], completed_at=now())
        receipt["accepted_prediction_diagnostic"] = True
        atomic_json(folder / "products/receipt.json", receipt)
        return receipt
    except BaseException as error:
        record.update(status="FAIL", completed_at=now(), error=type(error).__name__ + ": " + str(error))
        atomic_json(folder / "invocation.json", record)
        raise


def render_products(m, report, output):
    """Render Julia values and separate measured errors; no likelihood or fitting."""
    import csv
    request = read_json(read_json(m["source_plan"])["request"])
    if (request.get("theory_scope") != "reconstructed_theory" or report.get("theory_scope") != "reconstructed_theory" or
            "fixture_scope" in request or "fixture_scope" in report):
        raise ValueError("the reconstructed-data renderer cannot label synthetic or undeclared results as physical source comparisons")
    source, predicted = request["points"], report["points"]
    if (len(predicted) != 72 or [p["id"] for p in predicted] != [p["id"] for p in source] or
            any(report.get(key) is not False for key in FLAGS)):
        raise ValueError("complete ordered source-bound72 predictions required")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    output = canonical(output)
    output.mkdir(parents=True, exist_ok=False)
    table = output / "comparison72.csv"
    with table.open("x", encoding="utf-8", newline="") as stream:
        stream.write("# reconstructed_theory; conditional fixed rounded parameters; no fit or likelihood\n")
        writer = csv.writer(stream)
        writer.writerow(["point_id", "ratio", "z1", "z2", "Q2_GeV2", "reported_pT_mean_GeV", "angular_mean_of_ratio",
            "qT_lower_GeV", "qT_upper_GeV", "P_lower_GeV", "P_upper_GeV", "measured", "statistical_error",
            "reported_systematic", "systematic_components_json", "fixed_theory"])
        for p, prediction in zip(source, predicted):
            k, measured = p["kinematics"], p["measured"]
            writer.writerow([p["id"], p["ratio"], k["z1"], k["z2"], k["Q2_GeV2"], p["reported_pT_mean_GeV"],
                p["source_angular"]["value"], *prediction["qT_bounds_GeV"], *prediction["mapped_P_bounds_GeV"],
                measured["central_value"], measured["statistical_error"], measured["systematic_components"]["reported_systematic"],
                json.dumps(measured["systematic_components"], sort_keys=True), prediction["amplitude"]])
    products = [table]
    bins = sorted({tuple(p["primary_bins"]["z1"]) for p in source})
    for ratio in ("UL", "UC"):
        figure, axes = plt.subplots(len(bins), 1, figsize=(7, 2.5 * len(bins)), squeeze=False)
        for axis, interval in zip(axes[:, 0], bins):
            indices = [i for i, p in enumerate(source) if p["ratio"] == ratio and tuple(p["primary_bins"]["z1"]) == interval]
            xs = [source[i]["kinematics"]["z2"] for i in indices]
            axis.errorbar(xs, [source[i]["measured"]["central_value"] for i in indices],
                yerr=[source[i]["measured"]["statistical_error"] for i in indices], fmt="o", capsize=3, label="Measured; statistical errors")
            axis.scatter(xs, [predicted[i]["amplitude"] for i in indices], marker="D", facecolors="none", edgecolors="#c25b12",
                         label="Fixed conditional integral")
            axis.set(title=f"{ratio}; source z1 bin [{interval[0]}, {interval[1]}]", ylabel="Cosine amplitude (fraction)")
            axis.grid(alpha=.2)
        axes[0, 0].legend(fontsize=8)
        axes[-1, 0].set_xlabel("Reported mean z2")
        figure.suptitle("BaBar RF0: reconstructed fixed-parameter diagnostic")
        figure.text(.5, .01, "qT [0,3.5] GeV; fixed z/angular means and uniform momentum acceptance.\n"
            "Systematic errors remain separate in the table; covariance and joint acceptance unresolved. No fit.", ha="center", fontsize=8)
        figure.tight_layout(rect=(0, .045, 1, .97))
        image = output / (ratio + "-comparison.png")
        figure.savefig(image, dpi=160)
        plt.close(figure)
        products.append(image)
    return products
