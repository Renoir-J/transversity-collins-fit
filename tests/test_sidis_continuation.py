"""Synthetic transport/corruption controls; no physical gate is created."""
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from collinsflow import sidis_continuation as c


@pytest.fixture
def namespace(tmp_path):
    original = tmp_path/"original"
    run = tmp_path/"continuation"
    original.mkdir()
    run.mkdir()
    return dict(original=str(original), run=str(run))


def test_only_original_prefix_resolves_to_old_canonical_paths(namespace):
    for phase in c.POINT_PHASES:
        for i in range(1, 133):
            path = c.point_folder(namespace, phase, i)
            expected = Path(namespace["original"])/"full" if phase == "native123" and i <= 13 else Path(namespace["run"])
            assert path == expected/phase/f"p{i:04d}"
    for phase in c.GLOBAL_PHASES:
        assert c.point_folder(namespace, phase) == Path(namespace["run"])/phase


@pytest.mark.parametrize("phase,index", [("native123", True), ("native4", 0), ("native4", 133),
    ("native123", 1.0), ("native123", "1"), ("diagnostic", True), ("diagnostic", 1), ("../native123", 1)])
def test_no_ambiguous_or_traversing_source_key(namespace, phase, index):
    with pytest.raises(ValueError):
        c.point_folder(namespace, phase, index)


def test_original_complete_and_partial_outputs_cannot_be_written(namespace):
    old = Path(namespace["original"])/"full/native123/p0014"
    old.mkdir(parents=True)
    (old/"report.json.progress.json").write_text("ORIGINAL_PARTIAL")
    before = c.file_hash(old/"report.json.progress.json")
    for i in range(1, 14):
        with pytest.raises(ValueError, match="original/aliased"):
            c.require_new_folder(namespace, "native123", i)
    new = c.require_new_folder(namespace, "native123", 14)
    new.mkdir(parents=True)
    (new/"report.json.progress.json").write_text("NEW_INTERRUPTED")
    with pytest.raises(ValueError, match="occupied"):
        c.require_new_folder(namespace, "native123", 14)
    assert c.file_hash(old/"report.json.progress.json") == before


def test_alias_to_original_is_not_a_fresh_new_folder(namespace):
    original = Path(namespace["original"])
    dest = Path(namespace["run"])/"native123"
    dest.symlink_to(original, target_is_directory=True)
    with pytest.raises(ValueError, match="original/aliased"):
        c.require_new_folder(namespace, "native123", 14)


@pytest.mark.parametrize("mutation", ("bytes", "size", "added", "missing", "duplicate", "escape", "relative", "alias"))
def test_exact_inventory_rejects_changed_originals(tmp_path, mutation):
    original = tmp_path/"original"
    original.mkdir()
    path = original/"a.json"
    path.write_text("ORIGINAL")
    rows = [dict(path=str(path), bytes=path.stat().st_size, sha256=c.file_hash(path))]
    c.verify_inventory(original, rows)
    if mutation == "bytes":
        path.write_text("CORRUPT!")
    elif mutation == "size":
        rows[0]["bytes"] += 1
    elif mutation == "added":
        (original/"extra").write_text("EXTRA")
    elif mutation == "missing":
        path.unlink()
    elif mutation == "duplicate":
        rows.append(dict(rows[0]))
    elif mutation == "escape":
        rows[0]["path"] = str(tmp_path/"other")
    elif mutation == "relative":
        rows[0]["path"] = "a.json"
    else:
        alias = original/"alias"
        alias.symlink_to(path)
        rows.append(dict(rows[0], path=str(alias)))
    with pytest.raises(ValueError):
        c.verify_inventory(original, rows)


def test_raw_git_bytes_are_checked_even_with_assume_unchanged(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    path = tmp_path/"fixture.txt"
    path.write_text("EXPLICITLY SYNTHETIC\n")
    subprocess.run(["git", "add", "fixture.txt"], cwd=tmp_path, check=True)
    subprocess.run(["git", "-c", "user.name=Synthetic Test", "-c", "user.email=synthetic@example.invalid",
                    "commit", "-qm", "synthetic fixture"], cwd=tmp_path, check=True)
    head, hashes = c.verify_git(tmp_path)
    assert set(hashes) == {"fixture.txt"}
    assert c.verify_git(tmp_path, head) == (head, hashes)
    subprocess.run(["git", "update-index", "--assume-unchanged", "fixture.txt"], cwd=tmp_path, check=True)
    path.write_text("CORRUPTED SYNTHETIC\n")
    with pytest.raises(ValueError, match="working bytes"):
        c.verify_git(tmp_path, head)


def test_contract_has_no_automatic_physics_or_production_admission():
    root = Path(__file__).resolve().parents[1]
    policy = c.contract(root)
    assert all(policy[f] is False for f in c.FLAGS)
    assert policy["supported_origin"]["completed_native123_prefix"] == 13
    assert policy["supported_origin"]["prior_attempt1_records_accepted"] is False
    assert "all132" in policy["continuation"]["all_phases_required"]


def test_parent_input_is_original_point13_not_a_relocated_receipt(namespace, tmp_path):
    namespace.update(frozen=str(tmp_path/"kernel"), adapter=str(tmp_path/"adapter"),
        request=str(tmp_path/"request.json"), native=str(tmp_path/"native.json"),
        parton=str(tmp_path/"parton.json"), adapter_files={})
    api = SimpleNamespace(TIME_BINARY=Path("/usr/bin/time"))
    paths = c.input_paths(namespace, tmp_path/"manifest.json", "native123", 14, api)
    assert Path(namespace["original"])/"full/native123/p0013/invocation.json" in paths
    assert Path(namespace["run"])/"native123/p0013/invocation.json" not in paths
    assert Path(namespace["run"])/"prefix-review/invocation.json" in paths


def test_forged_complete_marker_does_not_make_a_failed_job_complete(namespace, tmp_path):
    namespace.update(numerical_source_fingerprint="synthetic", adapter_commit="synthetic")
    path = tmp_path/"manifest.json"
    path.write_text(json.dumps({"scope": "explicitly_synthetic"}))
    folder = c.point_folder(namespace, "native123", 14)
    folder.mkdir(parents=True)
    (folder/"job.completed").touch()
    c.atomic_json(folder/"invocation.json", dict(status="FAIL"))
    api = SimpleNamespace(COMPLETED="COMPLETE_JOB_RECORDED_NOT_PRODUCTION_GATE")
    with pytest.raises(ValueError, match="missing, failed"):
        c.verify_job(namespace, path, "native123", 14, api)
    assert c.read_json(folder/"invocation.json")["status"] == "FAIL"


@pytest.mark.parametrize("outcome", ("exit7", "cooperative_timeout", "exit0_invalid_reports"))
def test_real_complete_command_failures_retain_new_files_and_original_bytes(namespace, tmp_path, monkeypatch, outcome):
    """Isolate controller lifetime using synthetic inputs and real GNU-time jobs."""
    from collinsflow import historical_sidis_jobs as jobs
    m = dict(namespace, frozen=str(tmp_path/"synthetic-kernel"), adapter=str(tmp_path/"synthetic-adapter"),
        numerical_source_fingerprint="EXPLICITLY_SYNTHETIC", adapter_commit="EXPLICITLY_SYNTHETIC", adapter_files={})
    path = Path(m["run"])/"continuation.json"
    path.write_text('{"scope":"explicitly_synthetic_controller_test"}')
    m["request"] = str(tmp_path/"synthetic-request.json")
    m["request_sha256"] = "EXPLICITLY_SYNTHETIC"
    m["native"] = str(tmp_path/"synthetic-native.json")
    m["parton"] = str(tmp_path/"synthetic-parton.json")
    budget = dict(whole_Julia_jobs_time_budget_seconds=5, maximum_command_RSS_KiB=3*1024**2,
                  minimum_available_memory_bytes=0)
    api = SimpleNamespace(**{name: getattr(jobs, name) for name in ("SCOPE", "POLICY_SHA", "TIME_BINARY",
        "_pilot_lock", "_timed_job", "_require_command_success", "_deadline", "_report_audit", "utc_now")})
    api.load_policy = lambda: {"resources": {"native123": budget}}
    for item in c.input_paths(m, path, "native123", 14, api):
        if item == jobs.TIME_BINARY or item == path:
            continue
        item.parent.mkdir(parents=True, exist_ok=True)
        item.write_text('{"scope":"EXPLICITLY_SYNTHETIC"}')
    shared = Path(m["frozen"])/".runtime/sidis-pilot.lock"
    shared.parent.mkdir()
    shared.write_text("SYNTHETIC_ORIGINAL_LOCK_BYTES")
    original = c.point_folder(m, "native123", 13)/"invocation.json"
    old_hash, lock_hash = c.file_hash(original), c.file_hash(shared)
    # Actual source/physics preflight is deliberately absent in this lifetime
    # fixture. It can never be used by prepare(), prefix(), or a physical job.
    monkeypatch.setattr(c, "verify_manifest", lambda *a: m)
    monkeypatch.setattr(c, "verify_prefix", lambda *a: {})
    monkeypatch.setattr(c, "verify_job", lambda *a: {})
    folder = c.point_folder(m, "native123", 14)
    def command(m, path, julia, program, *args):
        dest = folder/("audit.json" if "audit_" in program else "report.json")
        progress = folder/"report.json.progress.json"
        code = (f"from pathlib import Path; Path({str(dest)!r}).write_text('{{\"status\":\"FAIL\",\"scope\":\"explicitly_synthetic\"}}');"
                f"Path({str(progress)!r}).write_text('SYNTHETIC_PARTIAL');")
        if outcome == "exit7":
            code += "raise SystemExit(7)"
        elif outcome == "cooperative_timeout":
            code += ("import signal,time,sys; signal.signal(signal.SIGTERM,lambda *args:sys.exit(0));"
                     "signal.signal(signal.SIGINT,lambda *args:sys.exit(0)); time.sleep(20)")
        return [sys.executable, "-c", code]
    monkeypatch.setattr(c, "_julia_command", command)
    with pytest.raises(ValueError):
        c.run_job(m, path, "native123", 14, api, sys.executable,
                  .4 if outcome == "cooperative_timeout" else None)
    inv = c.read_json(folder/"invocation.json")
    assert inv["status"] == "FAIL" and inv["completed_at"]
    assert all(inv[f] is False for f in c.FLAGS)
    assert (folder/"initial-invocation.json").is_file()
    assert (folder/"report.json.progress.json").is_file()
    assert (folder/"producer.log").is_file() and (folder/"producer.resources.txt").is_file()
    assert not (folder/"job.completed").exists()
    assert c.file_hash(original) == old_hash and c.file_hash(shared) == lock_hash
    if outcome == "cooperative_timeout":
        assert inv["steps"][0]["returncode"] == 124
        assert inv["steps"][0]["process_returncode"] == 0
    elif outcome == "exit7":
        assert inv["steps"][0]["returncode"] == 7
    else:
        assert len(inv["steps"]) == 2 and all(s["returncode"] == 0 for s in inv["steps"])
