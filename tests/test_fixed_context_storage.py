"""Synthetic controller failures never become physical storage receipts."""
import math
import os
from pathlib import Path
import signal
import subprocess

import pytest

from collinsflow import fixed_context_storage as storage
from collinsflow import historical_sidis_pilot as pilot
from collinsflow.provenance import file_hash, read_json


def test_frozen_storage_contract_and_budget():
    policy = storage.storage_policy()
    assert len(policy["validation"]["physical_cases"]) == 4
    assert all(policy[k] is False for k in storage.FLAGS)
    assert storage._budget(policy)["resources"] == {
        "minimum_available_memory_bytes": 8 * 1024 ** 3,
        "maximum_command_RSS_KiB": 3 * 1024 ** 2,
        "whole_Julia_jobs_time_budget_seconds": 3600,
    }
    for value in (True, 0, -1, math.nan, math.inf, 3601):
        with pytest.raises(ValueError):
            storage._deadline(policy, value)
    for value in (True, 0, 5, "1"):
        with pytest.raises(ValueError):
            storage._case(policy, value)


@pytest.fixture
def synthetic_failure_inputs(tmp_path, monkeypatch):
    paths = []
    for name in ("sidis", "sia", "native", "parton"):
        path = tmp_path / (name + ".json")
        path.write_text('{"scope":"explicitly_synthetic_controller_failure_control"}\n')
        paths.append(path)
    inputs = {str(p): file_hash(p) for p in paths}
    monkeypatch.setattr(storage, "_sources", lambda *args: ("1" * 64, inputs))
    monkeypatch.setattr(storage, "git_commit", lambda root: "synthetic-failure-control-not-a-physical-commit")
    original = subprocess.check_output

    def output(command, **kwargs):
        if command == ["git", "status", "--porcelain"]:
            return ""
        return original(command, **kwargs)

    monkeypatch.setattr(storage.subprocess, "check_output", output)
    monkeypatch.setattr(pilot, "available_memory_bytes", lambda: 16 * 1024 ** 3)
    return paths


@pytest.mark.skipif(os.name != "posix" or not Path("/usr/bin/time").is_file(), reason="Linux GNU-time controller integration")
@pytest.mark.parametrize("mode", ["fail", "timeout", "zero_without_accepted_reports"])
def test_complete_case_retains_failed_or_interrupted_files(tmp_path, synthetic_failure_inputs, mode):
    worker = tmp_path / "synthetic-julia-failure-control"
    worker.write_text("#!/usr/bin/env python3\n"
        "import pathlib,signal,sys,time\n"
        "out=pathlib.Path(sys.argv[-1]);out.parent.mkdir(parents=True,exist_ok=True)\n"
        "out.write_text('{\"scope\":\"synthetic_failure_control\",\"status\":\"FAIL\"}')\n"
        "(out.parent/'synthetic-partial.bin').write_bytes(b'preserve-synthetic-partial')\n"
        "print('explicitly synthetic complete-job failure control',flush=True)\n"
        + ("sys.exit(2)\n" if mode == "fail" else
           "sys.exit(0)\n" if mode == "zero_without_accepted_reports" else
           "signal.signal(signal.SIGINT,lambda *args:sys.exit(0))\nwhile True:time.sleep(.01)\n"))
    worker.chmod(0o755)
    run = tmp_path / "run"
    with pytest.raises((ValueError, FileNotFoundError)):
        storage.run_case(*synthetic_failure_inputs, run, 1, julia=str(worker),
                         timeout_seconds=.5 if mode == "timeout" else 10)
    folder = run / "case01"
    inv = read_json(folder / "invocation.json")
    assert inv["status"] == "FAIL" and inv["completed_at"]
    assert all(inv[k] is False for k in storage.FLAGS)
    assert not (folder / "job.completed").exists()
    assert (folder / "synthetic-partial.bin").read_bytes() == b"preserve-synthetic-partial"
    assert (folder / "initial-invocation.json").is_file()
    assert (folder / "producer.log").is_file()
    assert (folder / "producer.resources.txt").is_file()
    if mode == "timeout":
        assert inv["steps"][0]["returncode"] == 124
        assert inv["steps"][0]["process_returncode"] == 0
        assert inv["steps"][0]["interruption"] == "resource_timeout"
    elif mode == "zero_without_accepted_reports":
        assert [s["returncode"] for s in inv["steps"]] == [0, 0]
    preserved = {p: file_hash(p) for p in folder.iterdir() if p.is_file()}
    with pytest.raises(ValueError, match="preserve earlier"):
        storage.run_case(*synthetic_failure_inputs, run, 1, julia=str(worker))
    assert all(file_hash(p) == h for p, h in preserved.items())


def test_forged_four_markers_do_not_dispatch_or_accept(tmp_path, synthetic_failure_inputs, monkeypatch):
    run = tmp_path / "forged"
    for i in range(1, 5):
        folder = run / f"case{i:02d}"
        folder.mkdir(parents=True)
        (folder / "job.completed").touch()

    def forbidden(*args, **kwargs):
        raise AssertionError("a forged marker caused numerical dispatch")

    monkeypatch.setattr(storage, "_timed_job", forbidden)
    with pytest.raises(FileNotFoundError):
        storage.build_or_verify_ledger(*synthetic_failure_inputs, run, run / "ledger.json")
    assert not (run / "ledger.json").exists()


def test_parent_signal_with_zero_child_is_retained(tmp_path, synthetic_failure_inputs, monkeypatch):
    """Inject only a synthetic wrapper outcome; real signals are tested by job_control."""
    def interrupted(name, command, root, folder, tool, policy, deadline):
        log, resource = folder / (name + ".log"), folder / (name + ".resources.txt")
        log.write_text("explicitly synthetic parent-signal retention control\n")
        resource.write_text("COLLINSFLOW_TIME_V1\t1\t0.01\t0\t0\t0\n")
        return {"name": name, "command": command, "returncode": 143, "process_returncode": 0,
                "interruption": "parent_signal", "termination_signal": signal.SIGTERM,
                "wrapper_elapsed_seconds": .01, "log_sha256": file_hash(log),
                "resource_sha256": file_hash(resource), "resource": pilot.parse_resource_file(resource)}

    monkeypatch.setattr(storage, "_timed_job", interrupted)
    with pytest.raises(ValueError, match="interrupted"):
        storage.run_case(*synthetic_failure_inputs, tmp_path / "signal", 1, julia=os.sys.executable)
    folder = tmp_path / "signal/case01"
    inv = read_json(folder / "invocation.json")
    assert inv["status"] == "FAIL" and inv["completed_at"]
    assert inv["steps"][0]["returncode"] == 143
    assert inv["steps"][0]["process_returncode"] == 0
    assert not (folder / "job.completed").exists()


@pytest.mark.skipif(os.name != "posix", reason="Linux memory admission")
def test_missing_memory_does_not_dispatch(tmp_path, synthetic_failure_inputs, monkeypatch):
    monkeypatch.setattr(pilot, "available_memory_bytes", lambda: 0)
    with pytest.raises(ValueError, match="insufficient available memory"):
        storage.run_case(*synthetic_failure_inputs, tmp_path / "memory", 1, julia=os.sys.executable)
    folder = tmp_path / "memory/case01"
    assert read_json(folder / "invocation.json")["status"] == "FAIL"
    assert not (folder / "producer.log").exists()
    assert not (folder / "job.completed").exists()
