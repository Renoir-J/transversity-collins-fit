"""Synthetic complete-job metadata and resource controls; no Julia physical work."""
import hashlib
import copy
import json
import sys

import pytest

from collinsflow import sia_source_controls as controls


def test_gnu_time_keeps_a_nonzero_child_exit_and_scope(tmp_path):
    path = tmp_path / "resources.txt"
    path.write_text("Command exited with non-zero status7\nCOLLINSFLOW_TIME_V1\t12345\t1.25\t0.9\t0.1\t7\n", encoding="utf-8")
    result = controls.resource_record(path)
    assert result["child_exit_code"] == 7
    assert result["maximum_RSS_KiB"] == 12345
    assert result["wall_seconds"] == 1.25
    assert "not_simultaneous_process_tree_peak" in result["scope"]


@pytest.mark.parametrize("contents", ["", "COLLINSFLOW_TIME_V1\t1\t1\t1\t1\n",
    "COLLINSFLOW_TIME_V1\t0\t1\t1\t1\t0\n", "COLLINSFLOW_TIME_V1\t1\tnan\t1\t1\t0\n",
    "COLLINSFLOW_TIME_V1\t1\t1\t-1\t1\t0\n", "COLLINSFLOW_TIME_V1\t1\t1\t1\t1\t-2\n",
    "COLLINSFLOW_TIME_V1\t1\t1\t1\t1\t0\nCOLLINSFLOW_TIME_V1\t2\t2\t2\t2\t0\n"])
def test_incomplete_or_nonfinite_resource_records_reject(tmp_path, contents):
    path = tmp_path / "resources.txt"
    path.write_text(contents, encoding="utf-8")
    with pytest.raises(ValueError):
        controls.resource_record(path)


def test_worker_command_keeps_the_juliaup_launcher_and_fixed_program(tmp_path):
    launcher = "/explicitly-synthetic/.juliaup/bin/julia"
    plan = dict(frozen="/explicitly-synthetic/frozen")
    command = controls._command(plan, tmp_path / "plan.json", "bound-plan", tmp_path, "producer", launcher)
    assert command[6] == launcher
    assert command[8] == "--project=/explicitly-synthetic/frozen/numerics"
    assert command[9] == str(controls.ROOT / "numerics/bin/run_sia_source_controls.jl")
    assert command[-1] == str(tmp_path / "report.json")
    audit = controls._command(plan, tmp_path / "plan.json", "bound-plan", tmp_path, "audit", launcher, "bound-report")
    assert audit[-2:] == ["bound-report", str(tmp_path / "audit.json")]


@pytest.mark.parametrize("field,value", [("returncode", False), ("process_returncode", False),
    ("returncode", 124), ("interruption", "resource_timeout"), ("termination_signal", 2),
    ("wrapper_elapsed_seconds", -1.), ("wrapper_elapsed_seconds", float("nan"))])
def test_raw_zero_or_invalid_metadata_does_not_promote_interruption(field, value):
    step = dict(returncode=0, process_returncode=0, interruption=None, termination_signal=None,
        wrapper_elapsed_seconds=1., resource=dict(child_exit_code=0, maximum_RSS_KiB=123))
    budget = dict(maximum_command_RSS_KiB=1000)
    assert controls.successful_command(step, budget)
    changed = copy.deepcopy(step)
    changed[field] = value
    assert not controls.successful_command(changed, budget)


def synthetic_invocation():
    return dict(status=controls.COMPLETED, scope=controls.SCOPE, theory_scope="synthetic_closure",
        record_role="complete_source_control_pair_invocation", completed_at="explicitly-synthetic-time",
        **dict.fromkeys(controls.FLAGS, False))


@pytest.mark.parametrize("flag", controls.FLAGS)
def test_forged_self_hashed_marker_cannot_claim_physical_admission(tmp_path, flag):
    folder = tmp_path / "pair"
    folder.mkdir()
    record = synthetic_invocation()
    record[flag] = True
    path = folder / "invocation.json"
    path.write_text(json.dumps(record), encoding="utf-8")
    (folder / "job.completed").write_text(hashlib.sha256(path.read_bytes()).hexdigest() + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="marker alone"):
        controls.verify_pair(folder, tmp_path / "plan.json", "external-plan", julia=sys.executable)


def test_failed_invocation_and_wrong_marker_are_not_complete(tmp_path):
    folder = tmp_path / "pair"
    folder.mkdir()
    record = synthetic_invocation()
    record["status"] = "FAIL"
    path = folder / "invocation.json"
    path.write_text(json.dumps(record), encoding="utf-8")
    (folder / "job.completed").write_text(hashlib.sha256(path.read_bytes()).hexdigest(), encoding="utf-8")
    with pytest.raises(ValueError, match="marker alone"):
        controls.verify_pair(folder, tmp_path / "plan.json", "external-plan", julia=sys.executable)
    path.write_text(json.dumps(synthetic_invocation()), encoding="utf-8")
    (folder / "job.completed").write_text("0" * 64, encoding="utf-8")
    with pytest.raises(ValueError, match="marker alone"):
        controls.verify_pair(folder, tmp_path / "plan.json", "external-plan", julia=sys.executable)


def test_occupied_source_control_output_is_preserved(tmp_path, monkeypatch):
    planner = controls.planner()
    monkeypatch.setattr(controls, "planner", lambda: planner)
    monkeypatch.setattr(controls, "ROOT", tmp_path)
    plan = tmp_path / "runs/source/plan.json"
    plan.parent.mkdir(parents=True)
    plan.write_text("{}\n", encoding="utf-8")
    folder = plan.parent / "controls"
    folder.mkdir()
    sentinel = folder / "old.json"
    sentinel.write_text('{"status":"FAIL"}\n', encoding="utf-8")
    original = sentinel.read_bytes()
    with pytest.raises(ValueError, match="fresh canonical"):
        controls.run_pair(plan, planner.file_hash(plan), folder, julia=sys.executable)
    assert sentinel.read_bytes() == original
    assert list(folder.iterdir()) == [sentinel]
