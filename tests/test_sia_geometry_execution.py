"""Synthetic command/resource/provenance controls; no numerical or experimental jobs."""
import copy
from types import SimpleNamespace

import pytest

from collinsflow import sia_geometry_execution as execution
from collinsflow.provenance import atomic_json, read_json
from test_sia_geometry_host import fake_proc, add_process, OWNER


def fixture_manifest(tmp_path):
    plan = tmp_path / "plan.json"
    atomic_json(plan, dict(frozen="/explicitly-synthetic/frozen"))
    controls = tmp_path / "explicit-synthetic-source-controls"
    controls.mkdir()
    atomic_json(controls / "report.json", dict(source_context=dict(fixture_scope="explicitly_synthetic_metadata_not_dispatch")))
    native, parton = tmp_path / "synthetic-native.json", tmp_path / "synthetic-parton.json"
    for path in (native, parton):
        atomic_json(path, dict(scope="explicitly_synthetic_metadata_not_a_native_build"))
    return dict(run=str(tmp_path), source_plan=str(plan), julia_launcher="/explicitly-synthetic/julia",
        source_controls=str(controls), native=str(native), parton=str(parton),
        execution_inputs={str(plan): execution.file_hash(plan)},
        resources=execution.jobs.load_policy()["resources"],
        sidis_owner=copy.deepcopy(OWNER), sidis_run_roots=["/explicitly-synthetic/sidis-run"])


def test_exact_complete_job_shapes_reject_subsets_and_alternate_paths(tmp_path):
    m = fixture_manifest(tmp_path)
    for phase in execution.PHYSICAL:
        indices = range(1, 37) if phase != "derivative" else (1, 15, 36)
        assert len({execution.location(m, phase, i) for i in indices}) == len(indices)
    for phase in execution.GLOBALS:
        assert execution.location(m, phase).parent == tmp_path / "global"
    for phase, index in (("native123", True), ("native123", 0), ("native123", 37), ("derivative", 2),
                         ("native-selection", 1), ("native123/../../other", 1)):
        with pytest.raises(ValueError):
            execution.location(m, phase, index)
    assert execution.dependencies(m, "native123") == []
    assert execution.dependencies(m, "native-select3") == [("native123", i) for i in range(1, 37)]


@pytest.mark.parametrize("needs4", [False, True])
def test_global_dependencies_cover_every_required_geometry(tmp_path, needs4):
    m = fixture_manifest(tmp_path)
    path = execution.location(m, "native-select3") / "report.json"
    path.parent.mkdir(parents=True)
    atomic_json(path, dict(selection=dict(status="NATIVE4_REQUIRED_NOT_DIAGNOSTIC_GATE" if needs4 else "COMMON_NATIVE_SELECTION_OBSERVATION_NOT_GATE")))
    values = execution.dependencies(m, "diagnostic")
    assert len(values) == (114 if needs4 else 78)
    assert len(set(values)) == len(values)
    assert [(p, i) for p, i in values if p == "derivative"] == [("derivative", i) for i in (1, 15, 36)]
    assert len([i for p, i in values if p == "tailref"]) == 36
    assert len([i for p, i in values if p == "native4"]) == (36 if needs4 else 0)


def synthetic_pair(tmp_path):
    m = fixture_manifest(tmp_path)
    manifest = tmp_path / "execution.json"
    atomic_json(manifest, dict(scope="explicitly_synthetic_metadata_fixture_not_physics", **m))
    sha = execution.file_hash(manifest)
    folder = execution.location(m, "native123", 1)
    folder.mkdir(parents=True)
    atomic_json(folder / "lease.json", dict(scope="explicitly_synthetic_saved_lease_not_live_host"))
    header = dict(schema_version=1, scope=execution.jobs.SCOPE, theory_scope="reconstructed_theory",
        source_context=execution.report_context(m, sha), manifest_sha256=sha, phase="native123", index=1,
        completed_at="explicitly-synthetic-time-not-a-measurement", fixture_scope="explicitly_synthetic_protocol_shape_not_dispatch",
        **dict.fromkeys(execution.FLAGS, False))
    atomic_json(folder / "report.json", dict(status="COMPLETE_SIA_JOB_VALUES_COMPUTED_NOT_GATE", **header))
    steps = []
    for name in ("producer", "audit"):
        (folder / (name + ".log")).write_text("explicitly synthetic command metadata; no Julia execution\n")
        (folder / (name + ".resources.txt")).write_text("COLLINSFLOW_TIME_V1\t100\t1.0\t0.5\t0.1\t0\n")
        step = dict(name=name, command=execution.command(m, manifest, sha, "native123", 1, folder, name, execution.file_hash(folder / "report.json")),
            returncode=0, process_returncode=0, interruption=None, termination_signal=None, wrapper_elapsed_seconds=1.1,
            minimum_available_memory_observed_before_command=2**40, resource=execution.controls.resource_record(folder / (name + ".resources.txt")),
            log_sha256=execution.file_hash(folder / (name + ".log")), resource_sha256=execution.file_hash(folder / (name + ".resources.txt")))
        steps.append(step)
    audit = dict(status="COMPLETE_SIA_SAVED_JOB_ARITHMETIC_REVIEWED_NOT_GATE", **header,
        report_sha256=execution.file_hash(folder / "report.json"), complete_prerequisite_pair_ledgers_verified=0,
        native_evolution_reexecuted_by_reader=False, cold_AD_or_integrals_reexecuted_by_reader=False,
        raw_primary_tables_replayed_in_Julia=False)
    atomic_json(folder / "audit.json", audit)
    initial = dict(status="FAIL", scope=execution.PAIR_SCOPE, theory_scope="reconstructed_theory", phase="native123", index=1,
        resources=m["resources"]["native123"], started_at="explicit-synthetic-time-not-measurement", **dict.fromkeys(execution.FLAGS, False),
        record_role="initial_sentinel_not_an_executed_failure", steps=[], completed_at=None, manifest_sha256=sha)
    atomic_json(folder / "initial-invocation.json", initial)
    record = dict(status=execution.COMPLETE, scope=execution.PAIR_SCOPE, completed_at="synthetic-time-not-measurement",
        theory_scope="reconstructed_theory", record_role="complete_job_execution_invocation",
        **dict.fromkeys(execution.FLAGS, False), manifest_sha256=sha, phase="native123", index=1,
        resources=m["resources"]["native123"], steps=steps, initial_record_sha256=execution.file_hash(folder / "initial-invocation.json"),
        inputs=execution.job_inputs(m, manifest, sha, "native123", 1),
        outputs={str(folder / name): execution.file_hash(folder / name) for name in execution.PRODUCTS})
    atomic_json(folder / "invocation.json", record)
    (folder / "job.completed").write_text(execution.file_hash(folder / "invocation.json") + "\n")
    return m, manifest, sha, folder, record


def test_metadata_review_explicitly_does_not_claim_numerical_arithmetic(tmp_path):
    m, path, sha, _, record = synthetic_pair(tmp_path)
    assert execution.verify_pair(m, path, sha, "native123", 1) == record
    assert record["accepted_prediction_diagnostic"] is False
    assert record["accepted_fitting_gate"] is False
    with pytest.raises(ValueError):
        execution.load(path, sha)
    # The placeholder producer cannot pass the Julia saved-arithmetic reader;
    # this test covers only the scheduling metadata boundary.


@pytest.mark.parametrize("field,value", [("source_context", {}), ("theory_scope", "synthetic_closure"), ("scope", "production_fit"),
    ("schema_version", True), ("completed_at", ""), ("index", True), ("complete_prerequisite_pair_ledgers_verified", False),
    ("native_evolution_reexecuted_by_reader", True), ("cold_AD_or_integrals_reexecuted_by_reader", True),
    ("raw_primary_tables_replayed_in_Julia", True)])
def test_rehashed_audit_requires_complete_identity_and_reader_scope(tmp_path, field, value):
    m, path, sha, folder, record = synthetic_pair(tmp_path)
    audit = read_json(folder / "audit.json")
    audit[field] = value
    atomic_json(folder / "audit.json", audit)
    record["outputs"][str(folder / "audit.json")] = execution.file_hash(folder / "audit.json")
    atomic_json(folder / "invocation.json", record)
    (folder / "job.completed").write_text(execution.file_hash(folder / "invocation.json") + "\n")
    with pytest.raises(ValueError):
        execution.verify_pair(m, path, sha, "native123", 1)


@pytest.mark.parametrize("request_scope,report_scope,fixture", [("synthetic_closure", "synthetic_closure", True),
    ("reconstructed_theory", "synthetic_closure", False), ("synthetic_closure", "reconstructed_theory", False),
    ("reconstructed_theory", "reconstructed_theory", True), (None, "reconstructed_theory", False)])
def test_reconstructed_renderer_rejects_synthetic_or_undeclared_artifacts_before_writing(tmp_path, request_scope, report_scope, fixture):
    request = dict(theory_scope=request_scope, points=[dict(id=str(i)) for i in range(72)])
    report = dict(theory_scope=report_scope, points=request["points"], **dict.fromkeys(execution.FLAGS, False))
    if fixture:
        request["fixture_scope"] = "explicitly_synthetic_renderer_mutation_not_primary_data"
    request_path, plan_path = tmp_path / "request.json", tmp_path / "plan.json"
    atomic_json(request_path, request)
    atomic_json(plan_path, dict(request=str(request_path)))
    output = tmp_path / "must-remain-absent"
    with pytest.raises(ValueError, match="cannot label synthetic"):
        execution.render_products(dict(source_plan=str(plan_path)), report, output)
    assert not output.exists()


@pytest.mark.parametrize("mutation", ["effective-exit", "raw-exit", "signal", "interruption", "rss", "wall", "duration", "memory",
    "command", "resource-bytes", "output", "inputs", "sentinel", "missing-step", "scope", "gate", "marker"])
def test_rehashed_marker_cannot_replace_complete_pair_proof(tmp_path, mutation):
    m, path, sha, folder, record = synthetic_pair(tmp_path)
    step = record["steps"][0]
    if mutation == "effective-exit":
        step["returncode"] = 124
    elif mutation == "raw-exit":
        step["process_returncode"] = 2
    elif mutation == "signal":
        step["termination_signal"] = 2
    elif mutation == "interruption":
        step["interruption"] = "parent_signal"
    elif mutation == "rss":
        step["resource"]["maximum_RSS_KiB"] = 2**40
    elif mutation == "wall":
        step["resource"]["wall_seconds"] = 1e9
    elif mutation == "duration":
        step["wrapper_elapsed_seconds"] = 1e9
    elif mutation == "memory":
        step["minimum_available_memory_observed_before_command"] = 1
    elif mutation == "command":
        step["command"][-2] = "substituted-output"
    elif mutation == "resource-bytes":
        (folder / "producer.resources.txt").write_text("missing measurement\n")
    elif mutation == "output":
        (folder / "report.json").write_text("{}\n")
    elif mutation == "inputs":
        record["inputs"].clear()
    elif mutation == "sentinel":
        initial = read_json(folder / "initial-invocation.json")
        initial["completed_at"] = "invented-terminal-time"
        atomic_json(folder / "initial-invocation.json", initial)
        record["initial_record_sha256"] = execution.file_hash(folder / "initial-invocation.json")
    elif mutation == "missing-step":
        record["steps"].pop()
    elif mutation == "scope":
        record["scope"] = "production_fit"
    elif mutation == "gate":
        record["fit_performed"] = True
    elif mutation == "marker":
        record["status"] = "FAIL"
    atomic_json(folder / "invocation.json", record)
    (folder / "job.completed").write_text(execution.file_hash(folder / "invocation.json") + "\n")
    with pytest.raises((ValueError, OSError)):
        execution.verify_pair(m, path, sha, "native123", 1)


def test_live_or_ambiguous_sid_is_not_a_dispatch_receipt(tmp_path, monkeypatch):
    m = fixture_manifest(tmp_path)
    proc = fake_proc(tmp_path)
    add_process(proc, OWNER["pid"], OWNER["start_ticks"], OWNER["command"])
    real = execution.host.observe_host
    monkeypatch.setattr(execution.host, "observe_host", lambda owner, roots: real(owner, roots, proc=proc))
    with pytest.raises(ValueError, match="blocks SIA dispatch"):
        execution.require_clear(m)
    (proc / str(OWNER["pid"]) / "stat").unlink()
    with pytest.raises(ValueError):
        execution.require_clear(m)


def test_failed_process_keeps_raw_effective_codes_even_without_time_record(tmp_path, monkeypatch):
    m = fixture_manifest(tmp_path)
    result = SimpleNamespace(_asdict=lambda: dict(returncode=124, process_returncode=0, interruption="resource_timeout", termination_signal=None))
    def fake_run(command, *, stdout, **kwargs):
        stdout.write(b"explicitly synthetic process result; no numerical worker\n")
        return result
    monkeypatch.setattr(execution.controls, "available_memory", lambda: 2**40)
    monkeypatch.setattr(execution.controls, "_supervisor", lambda *args: SimpleNamespace(run_complete_process=fake_run))
    measured = execution.measured(m, ["synthetic-not-executed"], tmp_path, "audit", m["resources"]["global"], 10)
    assert measured["returncode"] == 124 and measured["process_returncode"] == 0
    assert measured["interruption"] == "resource_timeout" and "resource_error" in measured
    assert not execution.controls.successful_command(measured, m["resources"]["global"])


@pytest.mark.parametrize("program", ["run_sidis_recovery.py", "run_sidis_recovery_point.jl", "audit_sidis_recovery_job.jl",
    "run_sidis_continuation.py", "run_sidis_continuation_global.jl"])
def test_current_recovery_and_continuation_orphans_block_sia(tmp_path, program):
    proc = fake_proc(tmp_path)
    add_process(proc, 456, 700, ["synthetic-interpreter", "/synthetic/" + program])
    result = execution.host.observe_host(OWNER, proc=proc)
    assert result["status"] == "HOST_RESOURCE_GUARD_BLOCKED"
    assert result["expected_owner_absent"] is True
    assert result["blockers"][0]["reason"] == "physical_producer_or_pair_reader_present"
