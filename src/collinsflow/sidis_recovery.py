"""Schedule whole Julia jobs across immutable SIDIS recovery namespaces."""
from __future__ import annotations

from contextlib import contextmanager
import copy
import importlib
import importlib.util
import math
import os
from pathlib import Path
import subprocess
import uuid


def sibling(name):
    spec = importlib.util.spec_from_file_location(name + "_whole_job_adapter", Path(__file__).with_name(name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


lineage = sibling("sidis_recovery_lineage")
host = sibling("sidis_recovery_host")
locations = sibling("sidis_recovery_locations")
predecessor_records = sibling("sidis_predecessor_records")
old, paths = lineage.old, lineage.paths
FLAGS, POINT_PHASES, GLOBAL_PHASES = paths.FLAGS, paths.POINT_PHASES, paths.GLOBAL_PHASES
SCOPE, PREPARED = lineage.SCOPE, lineage.PREPARED
file_hash, read_json, atomic_json = old.file_hash, old.read_json, old.atomic_json
REVIEWED = "ALL_RETAINED_PAIRS_FRESHLY_REVIEWED_NOT_DIAGNOSTIC_GATE"
TIME_FORMAT = "COLLINSFLOW_TIME_V1\t%M\t%e\t%U\t%S\t%x"


def provenance(m, path):
    return dict(scope=SCOPE, manifest_sha256=file_hash(path), adapter_commit=m["adapter_commit"], contract_sha256=paths.POLICY_SHA)


class RecoveryContext(dict):
    """One owned manifest and location snapshot; never a cross-job physics cache."""

    __slots__ = ("_locations", "_origins", "_predecessor_files")

    def __init__(self, data, checked_view, origins=(), predecessor_files=()):
        if data["run"] != str(checked_view.root) or data["retained_jobs"] != checked_view.entries:
            raise ValueError("location snapshot does not belong to this verified manifest")
        dict.__init__(self, copy.deepcopy(data))
        object.__setattr__(self, "_locations", locations.snapshot(checked_view, paths))
        object.__setattr__(self, "_origins", tuple(origins))
        object.__setattr__(self, "_predecessor_files", tuple(sorted(dict(predecessor_files).items())))

    @property
    def locations(self):
        return self._locations

    def origin_context(self, path):
        return next((context for name, context in self._origins if name == str(path)), None)

    def _immutable(self, *args, **kwargs):
        raise TypeError("owned recovery context cannot be reassigned")

    __setitem__ = __delitem__ = clear = pop = popitem = setdefault = update = __ior__ = _immutable
    __setattr__ = __delattr__ = _immutable


def view(m):
    return m.locations if isinstance(m, RecoveryContext) else paths.RecoveryRunPaths(m["run"], m["retained_jobs"])


def point_folder(m, phase, index=0):
    return view(m).folder(phase, index)


def runtime_binding(julia):
    """Keep Juliaup argv[0] spelling; bind its bytes and the actual runtime separately."""
    path = Path(julia)
    if not path.is_absolute() or str(path) != julia or not path.is_file():
        raise ValueError("explicit absolute Julia launcher required")
    env = {**os.environ, **dict.fromkeys(("JULIA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"), "1")}
    result = subprocess.check_output([julia, "--startup-file=no", "-e", 'println(VERSION); println(realpath("/proc/self/exe"))'],
                                    env=env, text=True, timeout=60).splitlines()
    if len(result) != 2 or result[0] != "1.12.6":
        raise ValueError("the verified Julia1.12.6 runtime is required")
    executable = paths.canonical(result[1], directory=False)
    return dict(launcher=julia, launcher_sha256=file_hash(path), executable=str(executable),
                executable_sha256=file_hash(executable), version=result[0])


def verify_runtime(m):
    runtime = m["runtime"]
    if (runtime["version"] != "1.12.6" or not Path(runtime["launcher"]).is_absolute() or
            file_hash(runtime["launcher"]) != runtime["launcher_sha256"]):
        raise ValueError("changed Julia launcher or declared runtime")
    lineage.sha_file(runtime["executable"], runtime["executable_sha256"])


def source_replay(m):
    """A fresh explicit output avoids the legacy verifier's sibling temporary files."""
    module = importlib.import_module("collinsflow.historical_sidis_points")
    if Path(module.__file__).resolve().parents[2] != Path(m["frozen"]):
        raise ValueError("source preparation must use the explicit frozen package")
    folder = Path(m["run"]) / "source-replay"
    folder.mkdir(exist_ok=False)
    original = read_json(m["request"])
    src = original["source_preparation"]
    rec = dict(status="FAIL", scope="reconstructed_source_replay_only_not_numerical_admission", **dict.fromkeys(FLAGS, False))
    try:
        output = folder / "request.json"
        module.prepare_request(src["manifest"], src["cache"], src["candidates"], src["primary_cache"], output)
        if read_json(output) != original or file_hash(m["request"]) != m["request_sha256"]:
            raise ValueError("raw source replay disagrees with the exact original request")
        rec.update(status="CURRENT_RAW_SOURCE_REPLAY_VERIFIED_NEW_OUTPUT_ONLY_NOT_NUMERICAL_GATE",
                   original_request_sha256=m["request_sha256"], fresh_request_sha256=file_hash(output))
    except BaseException as error:
        rec["error"] = type(error).__name__ + ": " + str(error)
        raise
    finally:
        atomic_json(folder / "replay.json", rec)
    return dict(request=str(output), request_sha256=file_hash(output), receipt=str(folder / "replay.json"),
                receipt_sha256=file_hash(folder / "replay.json"))


def verify_source_replay(m):
    replay = m["source_replay"]
    root = Path(m["run"]) / "source-replay"
    if replay["request"] != str(root / "request.json") or replay["receipt"] != str(root / "replay.json"):
        raise ValueError("source replay must belong to the current namespace")
    lineage.sha_file(replay["request"], replay["request_sha256"])
    lineage.sha_file(replay["receipt"], replay["receipt_sha256"])
    receipt = read_json(replay["receipt"])
    if (read_json(replay["request"]) != read_json(m["request"]) or
            receipt.get("status") != "CURRENT_RAW_SOURCE_REPLAY_VERIFIED_NEW_OUTPUT_ONLY_NOT_NUMERICAL_GATE" or
            receipt.get("original_request_sha256") != m["request_sha256"] or
            receipt.get("fresh_request_sha256") != replay["request_sha256"] or
            any(receipt.get(key) is not False for key in FLAGS)):
        raise ValueError("changed, failed or self-admitting source replay")


def verify_manifest(path, frozen, predecessor, adapter, api, expected_sha):
    """Read-only verification; never call an ancestor's writing manifest verifier."""
    path = lineage.sha_file(path, expected_sha)
    m = read_json(path)
    paths.policy()
    if (path.name != "recovery.json" or m.get("status") != PREPARED or m.get("scope") != SCOPE or
            m.get("schema_version") != 2 or m.get("contract_sha256") != paths.POLICY_SHA or
            m.get("run") != str(path.parent) or any(m.get(key) is not False for key in FLAGS)):
        raise ValueError("wrong or self-admitting recovery manifest")
    source_binding = lineage.source_bindings(frozen, predecessor, adapter, m)
    if not path.parent.is_relative_to(Path(frozen) / "runs"):
        raise ValueError("current namespace must be inside the explicit frozen checkout's runs")
    state = lineage.validate_ancestors(m["ancestors"], frozen, predecessor, path.parent)
    if m["retained_jobs"] != state["entries"] or any(m.get(key) != state["first"][key] for key in lineage.ANCHOR_FIELDS):
        raise ValueError("changed original family or retained-job map")
    lineage.verify_original(state["first"], frozen, api)
    verify_runtime(m)
    verify_source_replay(m)
    origins = []
    for node, parent in zip(m["ancestors"], state["nodes"], strict=True):
        if node["protocol"] == "recovery_v2":
            parent_view = paths.RecoveryRunPaths(parent["run"], parent["retained_jobs"])
            origins.append((node["manifest"], RecoveryContext(parent, parent_view, origins, source_binding["predecessor_files"])))
    m = RecoveryContext(m, state["view"], origins, source_binding["predecessor_files"])
    for row in m["retained_jobs"]:
        verify_job(m, path, row["phase"], row["index"], api, state=state)
    return m


def _owner_from_record(record):
    owner = dict(pid=record["pid"], start_ticks=record.get("start_ticks", record.get("process_start_ticks")),
                 boot_id=record["boot_id"], command=record["command"])
    host._owner(owner)
    return owner


def guard(m):
    roots = [m["original"]] + [str(Path(node["manifest"]).parent) for node in m["ancestors"]]
    first = Path(roots[1]) / "launches/20261117T065752.754097Z/launch.json"
    lineage.sha_file(first, lineage.LAUNCH_SHA)
    owners = [_owner_from_record(read_json(first))]
    owners += [_owner_from_record(read_json(Path(root) / "driver/owner.json")) for root in roots[2:]]
    observations = [host.observe_host(owner, roots) for owner in owners]
    if any(row["status"] != "HOST_RESOURCE_GUARD_CLEAR_NOT_SCIENTIFIC_ADMISSION" or
           row["live_host_observed"] is not True or row["expected_owner_absent"] is not True for row in observations):
        raise ValueError("old owner, physical worker or ambiguous host observation blocks dispatch")
    return observations


def prepare(frozen, predecessor, adapter, parent_manifest, parent_sha, inventory_file, inventory_sha, output, api, julia):
    bound = lineage.source_bindings(frozen, predecessor, adapter)
    parent_manifest = lineage.sha_file(parent_manifest, parent_sha)
    lineage.sha_file(inventory_file, inventory_sha)
    parent = read_json(parent_manifest)
    output = Path(output)
    if (not output.is_absolute() or output.exists() or output.is_symlink() or output.resolve() != output or
            not output.is_relative_to(Path(frozen) / "runs")):
        raise ValueError("fresh canonical ignored recovery namespace required")
    if parent_manifest.name == "continuation.json":
        nodes = []
        protocol = "continuation_v1"
    elif parent_manifest.name == "recovery.json":
        nodes = parent["ancestors"]
        protocol = "recovery_v2"
    else:
        raise ValueError("declared predecessor manifest required")
    nodes = [*nodes, dict(protocol=protocol, manifest=str(parent_manifest), manifest_sha256=parent_sha,
                         inventory=str(inventory_file), inventory_sha256=inventory_sha)]
    # Create only the new root. Any preparation failure remains visible there.
    output.mkdir(parents=True, exist_ok=False)
    state = lineage.validate_ancestors(nodes, frozen, predecessor, output)
    lineage.verify_original(state["first"], frozen, api)
    m = {**{key: state["first"][key] for key in lineage.ANCHOR_FIELDS}, **bound,
         "schema_version": 2, "status": PREPARED, "scope": SCOPE, "contract_sha256": paths.POLICY_SHA,
         "run": str(output), "ancestors": nodes, "retained_jobs": state["entries"],
         "runtime": runtime_binding(julia), **dict.fromkeys(FLAGS, False)}
    m["origin_process_observations"] = guard(m)
    m["source_replay"] = source_replay(m)
    path = output / "recovery.json"
    atomic_json(path, m)
    verify_manifest(path, frozen, predecessor, adapter, api, file_hash(path))
    return path


def selection(m, phase, api):
    return api._selection(m["run"], phase, read_json(m["request"]), m["numerical_source_fingerprint"])


def input_paths(m, path, phase, index, api):
    """Complete invocation inputs; retained point reports keep their actual paths."""
    runview = view(m)
    frozen = Path(m["frozen"])
    inputs = [Path(m[key]) for key in ("request", "native", "parton")]
    inputs += [frozen / "config/historical-sidis-job-policy.json", api.TIME_BINARY,
               frozen / "scripts/run_historical_sidis_job.py", Path(path)]
    inputs += [Path(m["adapter"]) / name for name in m["adapter_files"]]
    inputs += [Path(m["runtime"][key]) for key in ("launcher", "executable")]
    inputs += [Path(m["source_replay"][key]) for key in ("request", "receipt")]
    for node in m["ancestors"]:
        inputs += [Path(node["manifest"]), Path(node["inventory"])]
    def points(name):
        inputs.extend(runview.folder(name, i) / "invocation.json" for i in range(1, 133))
    def glob(name, suffix="invocation.json"):
        inputs.append(runview.folder(name) / suffix)
    if phase in POINT_PHASES:
        if index > 1:
            inputs.append(runview.folder(phase, index - 1) / "invocation.json")
        if phase == "native4":
            glob("native-select3", "report.json")
            inputs.append(runview.folder("native123", index) / "report.json")
        if phase in ("tailref123", "reference4"):
            level = selection(m, "native-selection", api)["selected_level"]
            glob("native-selection", "report.json")
            inputs.append(runview.folder("native123" if level == 3 else "native4", index) / "report.json")
        if phase == "reference4":
            glob("derivative-select3", "report.json")
            inputs.append(runview.folder("tailref123", index) / "report.json")
    elif phase != "source-math":
        points("native123")
        if phase != "native-select3":
            glob("native-select3")
            if selection(m, "native-select3", api)["selected_level"] == 0:
                points("native4")
        if phase in ("derivative-select3", "diagnostic"):
            glob("native-selection")
            points("tailref123")
        if phase == "diagnostic":
            glob("derivative-select3")
            glob("source-math")
            if selection(m, "derivative-select3", api)["selected_reference_level"] == 0:
                points("reference4")
    else:
        glob("native-selection")
    inputs.append(Path(m["run"]) / "retained-review/invocation.json")
    return sorted(set(inputs))


def julia_command(m, path, program, *args):
    return [m["runtime"]["launcher"], "--startup-file=no", "--project=" + str(Path(m["frozen"]) / "numerics"),
            str(Path(m["adapter"]) / "numerics/bin" / program), m["frozen"], m["predecessor"], str(path),
            file_hash(path), *map(str, args)]


def measured_command(folder, name, command):
    return ["/usr/bin/time", "--output", str(folder / (name + ".resources.txt")), "--format", TIME_FORMAT, "--", *command]


def verify_step(step, folder, command, budget, api):
    if (any(type(step.get(key)) is not int or step[key] != 0 for key in ("returncode", "process_returncode")) or
            step.get("command") != measured_command(folder, step["name"], command)):
        raise ValueError("inexact or wrong complete-command identity/outcome")
    elapsed = step.get("wrapper_elapsed_seconds")
    if isinstance(elapsed, bool) or not isinstance(elapsed, (int, float)) or not math.isfinite(elapsed) or elapsed < 0:
        raise ValueError("finite nonnegative complete-command elapsed time required")
    api._require_command_success(step, budget)
    name = step["name"]
    if (file_hash(folder / (name + ".log")) != step["log_sha256"] or
            file_hash(folder / (name + ".resources.txt")) != step["resource_sha256"] or
            api.parse_resource_file(folder / (name + ".resources.txt")) != step["resource"]):
        raise ValueError("changed complete-command log or resource evidence")


def _new_job_metadata(m, path, phase, index, api):
    folder = point_folder(m, phase, index)
    inv = read_json(folder / "invocation.json")
    final = phase == "diagnostic"
    if (inv.get("status") != api.COMPLETED or inv.get("phase") != phase or type(inv.get("point_index")) is not int or
            inv["point_index"] != index or inv.get("recovery") != provenance(m, path) or not api._flags(inv, final) or
            inv.get("source_fingerprint") != m["numerical_source_fingerprint"] or inv.get("code_commit") != old.FROZEN_HEAD or
            inv.get("policy_sha256") != api.POLICY_SHA or inv.get("time_binary") != str(api.TIME_BINARY) or not inv.get("completed_at")):
        raise ValueError("missing, failed, stale or self-admitting recovery job")
    expected = {str(p): file_hash(p) for p in input_paths(m, path, phase, index, api)}
    if (inv["input_hashes"] != expected or set(inv["outputs"]) != api._products(folder, final) or
            any(file_hash(p) != sha for p, sha in inv["outputs"].items()) or
            file_hash(folder / "initial-invocation.json") != inv["initial_record_sha256"] or
            read_json(folder / "initial-invocation.json")["record_role"] != "initial_sentinel_not_an_executed_failure" or
            [step["name"] for step in inv["steps"]] != ["producer", "audit"]):
        raise ValueError("changed or incomplete inputs, outputs, sentinel or command ledger")
    budget = {"resources": api.load_policy()["resources"][phase if phase in POINT_PHASES else "global"]}
    program = "run_sidis_recovery_point.jl" if phase in POINT_PHASES else "run_sidis_recovery_global.jl"
    commands = [julia_command(m, path, program, phase, index, folder / "report.json"),
                julia_command(m, path, "audit_sidis_recovery_job.jl", phase, index, folder / "report.json", folder / "audit.json")]
    for step, command in zip(inv["steps"], commands, strict=True):
        verify_step(step, folder, command, budget, api)
    if sum(step["wrapper_elapsed_seconds"] for step in inv["steps"]) > budget["resources"]["whole_Julia_jobs_time_budget_seconds"]:
        raise ValueError("whole-job budget exhausted")
    report, audit = api._report_audit(folder, phase, index, read_json(m["request"]), m["numerical_source_fingerprint"], m["request_sha256"])
    if any(item.get("recovery") != provenance(m, path) for item in (report, audit)):
        raise ValueError("unbound recovery producer/reader")
    if final and (folder / "point-comparison.csv").read_text() != api.comparison_csv(report):
        raise ValueError("changed final comparison")
    return inv


def verify_job(m, path, phase, index, api, *, state=None):
    row = view(m).retained(phase, index) if phase in POINT_PHASES else None
    if row is None:
        inv = _new_job_metadata(m, path, phase, index, api)
    elif row["protocol"] == "frozen_original":
        inv = old._verify_original_job(m, index, read_json(m["request"]), api)
    else:
        origin_path = Path(row["producer_manifest"])
        origin = m.origin_context(origin_path) if isinstance(m, RecoveryContext) else None
        if origin is None:
            origin = read_json(origin_path)
        if row["protocol"] == "continuation_v1":
            if isinstance(m, RecoveryContext) and m._predecessor_files:
                inv = predecessor_records.verify(origin, origin_path, index, m._predecessor_files, old, api)
            else:
                inv = old.verify_job(origin, origin_path, phase, index, api)
        else:
            verify_runtime(origin)
            verify_source_replay(origin)
            inv = _new_job_metadata(origin, origin_path, phase, index, api)
            verify_retained_review(origin, origin_path, api)
    paths.canonical(point_folder(m, phase, index) / "job.completed", directory=False)
    return inv


def retained_hashes(m):
    return {str(Path(row["folder"]) / "invocation.json"): file_hash(Path(row["folder"]) / "invocation.json") for row in m["retained_jobs"]}


def review_inputs(m, path, api):
    inputs = [Path(path), api.TIME_BINARY, *[Path(m[key]) for key in ("request", "native", "parton")]]
    inputs += [Path(m["adapter"]) / name for name in m["adapter_files"]]
    inputs += [Path(m["runtime"][key]) for key in ("launcher", "executable")]
    inputs += [Path(m["source_replay"][key]) for key in ("request", "receipt")]
    for node in m["ancestors"]:
        inputs += [Path(node["manifest"]), Path(node["inventory"])]
    inputs += [Path(path) for path in retained_hashes(m)]
    return {str(path): file_hash(path) for path in sorted(set(inputs))}


def verify_retained_review(m, path, api):
    folder = Path(m["run"]) / "retained-review"
    inv = read_json(folder / "invocation.json")
    if (inv.get("status") != REVIEWED or inv.get("recovery") != provenance(m, path) or
            inv.get("input_hashes") != review_inputs(m, path, api) or any(inv.get(key) is not False for key in FLAGS)):
        raise ValueError("fresh complete independent retained-job reader required")
    budget = {"resources": api.load_policy()["resources"]["global"]}
    command = julia_command(m, path, "audit_sidis_recovery_retained.jl", folder / "audit.json")
    if inv["step"]["name"] != "retained-audit":
        raise ValueError("wrong retained-reader command")
    verify_step(inv["step"], folder, command, budget, api)
    if inv["step"]["wrapper_elapsed_seconds"] > budget["resources"]["whole_Julia_jobs_time_budget_seconds"]:
        raise ValueError("retained-reader budget exhausted")
    expected = {str(folder / name) for name in ("audit.json", "retained-audit.log", "retained-audit.resources.txt")}
    if set(inv["outputs"]) != expected or any(file_hash(name) != sha for name, sha in inv["outputs"].items()):
        raise ValueError("changed retained-reader outputs")
    audit = read_json(folder / "audit.json")
    if (audit.get("status") != "RETAINED_POINT_ARITHMETIC_REVIEWED_NOT_DIAGNOSTIC_GATE" or
            audit.get("recovery") != provenance(m, path) or audit.get("retained_invocations") != retained_hashes(m) or
            audit.get("jobs_reviewed") != [[row["phase"], row["index"]] for row in m["retained_jobs"]] or
            audit.get("source_fingerprint") != m["numerical_source_fingerprint"] or audit.get("request_sha256") != m["request_sha256"] or
            audit.get("current_global_prerequisites_accepted") is not False or any(audit.get(key) is not False for key in FLAGS)):
        raise ValueError("incomplete or self-admitting independent retained reader")
    return inv


@contextmanager
def pair_lock(m, api, family_lease):
    """The same physical resource lock used by SIA spans producer and auditor."""
    host.validate_owned_family_lease(family_lease)
    guard(m)
    with host.geometry_lock() as lease, old.job_lock(m, api):
        guard(m)
        yield lease


def recheck(m, path, api):
    if read_json(path) != m:
        raise ValueError("prepared recovery manifest changed during its owned lifetime")
    if isinstance(m, RecoveryContext):
        expected = tuple(tuple(row[key] for key in paths.FIELDS) for row in m["retained_jobs"])
        if m.locations.root != Path(m["run"]) or m.locations.rows != expected:
            raise ValueError("owned location snapshot changed during its whole-job context")
    return verify_manifest(path, m["frozen"], m["predecessor"], m["adapter"], api, file_hash(path))


def audit_retained(m, path, api, family_lease):
    recheck(m, path, api)
    folder = Path(m["run"]) / "retained-review"
    folder.mkdir(exist_ok=False)
    rec = dict(status="FAIL", recovery=provenance(m, path), family_resource_lease=family_lease,
               input_hashes=review_inputs(m, path, api), **dict.fromkeys(FLAGS, False))
    try:
        budget = {"resources": api.load_policy()["resources"]["global"]}
        with pair_lock(m, api, family_lease) as lease:
            rec["geometry_resource_lease"] = lease
            rec["step"] = api._timed_job("retained-audit", julia_command(m, path, "audit_sidis_recovery_retained.jl", folder / "audit.json"),
                Path(m["frozen"]), folder, api.TIME_BINARY, budget, api._deadline(budget["resources"]["whole_Julia_jobs_time_budget_seconds"], None))
        api._require_command_success(rec["step"], budget)
        recheck(m, path, api)
        rec["outputs"] = {str(folder / name): file_hash(folder / name) for name in ("audit.json", "retained-audit.log", "retained-audit.resources.txt")}
        rec["status"] = REVIEWED
        atomic_json(folder / "invocation.json", rec)
        verify_retained_review(m, path, api)
    except BaseException as error:
        rec.update(status="FAIL", error=type(error).__name__ + ": " + str(error))
        raise
    finally:
        rec["completed_at"] = api.utc_now()
        atomic_json(folder / "invocation.json", rec)
    return verify_retained_review(m, path, api)


def run_job(m, path, phase, index, api, family_lease, timeout_seconds=None):
    recheck(m, path, api)
    verify_retained_review(m, path, api)
    folder = view(m).fresh_folder(phase, index)
    if phase in POINT_PHASES and index > 1:
        verify_job(m, path, phase, index - 1, api)
    if phase == "source-math":
        verify_job(m, path, "native-selection", 0, api)
    folder.mkdir(parents=True, exist_ok=False)
    inv = dict(status="FAIL", scope=api.SCOPE, phase=phase, point_index=index, started_at=api.utc_now(), steps=[], outputs={},
               recovery=provenance(m, path), family_resource_lease=family_lease, **dict.fromkeys(FLAGS, False))
    atomic_json(folder / "initial-invocation.json", {**inv, "record_role": "initial_sentinel_not_an_executed_failure"})
    inv["initial_record_sha256"] = file_hash(folder / "initial-invocation.json")
    atomic_json(folder / "invocation.json", inv)
    try:
        inv.update(source_fingerprint=m["numerical_source_fingerprint"], code_commit=old.FROZEN_HEAD, policy_sha256=api.POLICY_SHA,
                   time_binary=str(api.TIME_BINARY), input_hashes={str(p): file_hash(p) for p in input_paths(m, path, phase, index, api)})
        budget = {"resources": api.load_policy()["resources"][phase if phase in POINT_PHASES else "global"]}
        program = "run_sidis_recovery_point.jl" if phase in POINT_PHASES else "run_sidis_recovery_global.jl"
        commands = (("producer", julia_command(m, path, program, phase, index, folder / "report.json")),
                    ("audit", julia_command(m, path, "audit_sidis_recovery_job.jl", phase, index, folder / "report.json", folder / "audit.json")))
        with pair_lock(m, api, family_lease) as lease:
            inv["geometry_resource_lease"] = lease
            deadline = api._deadline(budget["resources"]["whole_Julia_jobs_time_budget_seconds"], timeout_seconds)
            for name, command in commands:
                row = api._timed_job(name, command, Path(m["frozen"]), folder, api.TIME_BINARY, budget, deadline)
                inv["steps"].append(row)
                atomic_json(folder / "invocation.json", inv)
                api._require_command_success(row, budget)
        if phase == "diagnostic":
            api.render_comparison(read_json(folder / "report.json"), folder)
        recheck(m, path, api)
        if any(file_hash(p) != sha for p, sha in inv["input_hashes"].items()):
            raise ValueError("complete-job input changed during execution")
        inv["outputs"] = {str(p): file_hash(p) for p in api._products(folder, phase == "diagnostic")}
        inv["status"] = api.COMPLETED
        if phase == "diagnostic":
            inv.update(scope=api.SCIENTIFIC_SCOPE, accepted_prediction_diagnostic=True, all132_coverage_verified=True)
        inv["completed_at"] = api.utc_now()
        atomic_json(folder / "invocation.json", inv)
        _new_job_metadata(m, path, phase, index, api)
        (folder / "job.completed").touch(exist_ok=False)
    except BaseException as error:
        inv.update(status="FAIL", error=type(error).__name__ + ": " + str(error), **dict.fromkeys(FLAGS, False))
        raise
    finally:
        inv.setdefault("completed_at", api.utc_now())
        atomic_json(folder / "invocation.json", inv)
    return verify_job(m, path, phase, index, api)


def audit_cached(m, path, phase, index, api, family_lease):
    recheck(m, path, api)
    original = verify_job(m, path, phase, index, api)
    folder = Path(m["run"]) / "cache-reviews" / (phase + f"-{index}-" + uuid.uuid4().hex)
    folder.mkdir(parents=True, exist_ok=False)
    rec = dict(status="FAIL", recovery=provenance(m, path), phase=phase, point_index=index, family_resource_lease=family_lease,
               original_invocation_sha256=file_hash(point_folder(m, phase, index) / "invocation.json"), **dict.fromkeys(FLAGS, False))
    try:
        budget = {"resources": api.load_policy()["resources"]["global"]}
        command = julia_command(m, path, "audit_sidis_recovery_job.jl", phase, index, point_folder(m, phase, index) / "report.json", folder / "audit.json")
        with pair_lock(m, api, family_lease) as lease:
            rec["geometry_resource_lease"] = lease
            rec["step"] = api._timed_job("cache-audit", command, Path(m["frozen"]), folder, api.TIME_BINARY,
                budget, api._deadline(budget["resources"]["whole_Julia_jobs_time_budget_seconds"], None))
        verify_step(rec["step"], folder, command, budget, api)
        audit = read_json(folder / "audit.json")
        if (audit.get("status") != api.AUDITED or audit.get("recovery") != provenance(m, path) or
                audit.get("report_sha256") != file_hash(point_folder(m, phase, index) / "report.json") or
                audit.get("source_fingerprint") != m["numerical_source_fingerprint"] or audit.get("request_sha256") != m["request_sha256"] or
                audit.get("policy_sha256") != api.POLICY_SHA or not api._flags(audit, phase == "diagnostic")):
            raise ValueError("failed fresh current-prerequisite cached Julia review")
        recheck(m, path, api)
        if original != verify_job(m, path, phase, index, api):
            raise ValueError("original invocation changed during cached review")
        rec["outputs"] = {str(folder / name): file_hash(folder / name) for name in ("audit.json", "cache-audit.log", "cache-audit.resources.txt")}
        rec["status"] = "CURRENT_RECOVERY_ARITHMETIC_REVIEWED_NOT_FITTING"
    except BaseException as error:
        rec["error"] = type(error).__name__ + ": " + str(error)
        raise
    finally:
        rec["completed_at"] = api.utc_now()
        atomic_json(folder / "review.json", rec)
    return rec


def run_all(m, path, api, family_lease):
    """All132 common selections stay global; retained jobs never shorten the request."""
    def complete(phase, index=0):
        if (point_folder(m, phase, index) / "invocation.json").is_file():
            audit_cached(m, path, phase, index, api, family_lease)
        else:
            run_job(m, path, phase, index, api, family_lease)
        print(__import__("json").dumps(dict(phase=phase, index=index, status="COMPLETE_NOT_FITTING")), flush=True)
    if (Path(m["run"]) / "retained-review").exists():
        verify_retained_review(m, path, api)
    else:
        audit_retained(m, path, api, family_lease)
    for index in range(1, 133):
        # One fresh reader already reviewed every retained initial point.
        if view(m).retained("native123", index) is None:
            complete("native123", index)
    complete("native-select3")
    if selection(m, "native-select3", api)["selected_level"] == 0:
        for index in range(1, 133):
            complete("native4", index)
    complete("native-selection")
    complete("source-math")
    for index in range(1, 133):
        complete("tailref123", index)
    complete("derivative-select3")
    if selection(m, "derivative-select3", api)["selected_reference_level"] == 0:
        for index in range(1, 133):
            complete("reference4", index)
    complete("diagnostic")
    return audit_cached(m, path, "diagnostic", 0, api, family_lease)


@contextmanager
def driver(m, path, api, action):
    """One immutable owner identity for the complete scheduling lifetime."""
    with host.family_lock() as lease:
        guard(m)
        folder = Path(m["run"]) / "driver"
        folder.mkdir(exist_ok=False)
        identity = host.process_stat(Path("/proc/self/stat").read_text(encoding="utf-8"))
        owner = dict(pid=identity["pid"], start_ticks=identity["start_ticks"],
            boot_id=Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
            command=Path("/proc/self/cmdline").read_bytes().decode().rstrip("\0").split("\0"))
        atomic_json(folder / "owner.json", owner)
        rec = dict(status="IN_PROGRESS_IMMUTABLE_RECOVERY_NOT_A_GATE", action=action, recovery=provenance(m, path),
                   started_at=api.utc_now(), family_resource_lease=lease, **dict.fromkeys(FLAGS, False))
        atomic_json(folder / "invocation.json", rec)
        try:
            yield lease
            rec["status"] = "REQUESTED_RECOVERY_ACTION_COMPLETED_NOT_FITTING"
        except BaseException as error:
            rec.update(status="FAIL", error=type(error).__name__ + ": " + str(error))
            raise
        finally:
            rec["completed_at"] = api.utc_now()
            atomic_json(folder / "invocation.json", rec)


def snapshot(parent_manifest, parent_sha, output):
    """External stopped-run inventory; never rewrite an unknown historical exit."""
    parent_manifest = lineage.sha_file(parent_manifest, parent_sha)
    parent = read_json(parent_manifest)
    root = paths.canonical(parent_manifest.parent, directory=True)
    target = Path(output)
    if (parent.get("status") != PREPARED or target.exists() or target.resolve() != target or
            target.is_relative_to(root) or not target.is_absolute()):
        raise ValueError("fresh external snapshot for a declared recovery namespace required")
    owner = _owner_from_record(read_json(root / "driver/owner.json"))
    observation = host.observe_host(owner, [str(root)])
    if observation["status"] != "HOST_RESOURCE_GUARD_CLEAR_NOT_SCIENTIFIC_ADMISSION" or observation["live_host_observed"] is not True:
        raise ValueError("owner or orphan observation does not release this namespace")
    rows = lineage.inventory(root)
    rec = dict(status="STOPPED_SIDIS_NAMESPACE_INVENTORIED_NOT_A_GATE", manifest_sha256=parent_sha, run=str(root),
               retained_files=rows, owner=owner, old_owner_absent_verified=True, host_observation=observation,
               numerical_exit_codes_not_invented=True, **dict.fromkeys(FLAGS, False))
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("x", encoding="utf-8", newline="\n") as stream:
        __import__("json").dump(rec, stream, indent=2, allow_nan=False)
        stream.write("\n")
    lineage.verify_inventory(root, rows, root.parent)
    return rec
