import copy
import json
import os

import pytest

from collinsflow import sia_geometry_host as host

BOOT = "11111111-2222-4333-8444-555555555555"
OWNER = dict(pid=123, start_ticks=500, boot_id=BOOT, command=["/synthetic/python", "/synthetic/driver.py"])


def fake_proc(tmp_path):
    proc = tmp_path / "proc"
    boot = proc / "sys/kernel/random/boot_id"
    boot.parent.mkdir(parents=True)
    boot.write_text(BOOT + "\n")
    return proc


def add_process(proc, pid, start, command, *, state="S", parent=1):
    path = proc / str(pid)
    path.mkdir()
    # The comm field may contain both spaces and a closing parenthesis.
    fields = [state, str(parent), *["0"] * 17, str(start), *["0"] * 8]
    (path / "stat").write_text(f"{pid} (synthetic ) fixture) " + " ".join(fields) + "\n")
    (path / "cmdline").write_bytes("\0".join(command).encode() + (b"\0" if command else b""))
    return path


def test_synthetic_owner_is_not_release_and_never_live_host(tmp_path):
    proc = fake_proc(tmp_path)
    add_process(proc, OWNER["pid"], OWNER["start_ticks"], OWNER["command"])
    result = host.observe_host(OWNER, ["/synthetic/sidis-run"], proc=proc)
    assert result["status"] == "HOST_RESOURCE_GUARD_BLOCKED"
    assert result["expected_owner_live"] is True and result["expected_owner_absent"] is False
    assert result["live_host_observed"] is False
    assert result["input_scope"] == "explicitly_synthetic_proc_fixture_not_live_host"
    assert result["jobs_started_stopped_or_restarted"] is False
    assert all(result[key] is False for key in ("fit_performed", "production_ready", "accepted_fitting_gate"))
    assert json.loads(json.dumps(result))["blockers"][0]["start_ticks"] == 500


@pytest.mark.parametrize("kind", ["reused-pid", "zombie", "old-boot", "empty-kernel-command"])
def test_absent_owner_is_only_resource_observation(tmp_path, kind):
    proc = fake_proc(tmp_path)
    if kind == "reused-pid":
        add_process(proc, 123, 501, ["/unrelated/python"])
    elif kind == "zombie":
        add_process(proc, 123, 500, OWNER["command"], state="Z")
    elif kind == "old-boot":
        (proc / "sys/kernel/random/boot_id").write_text("99999999-8888-4777-8666-555555555555\n")
        add_process(proc, 123, 500, ["/unrelated/python"])
    else:
        add_process(proc, 124, 600, [], parent=2)
    result = host.observe_host(OWNER, proc=proc)
    assert result["status"] == "HOST_RESOURCE_GUARD_CLEAR_NOT_SCIENTIFIC_ADMISSION"
    assert result["expected_owner_absent"] is True and result["expected_owner_live"] is False
    assert result["production_ready"] is False and result["live_host_observed"] is False


@pytest.mark.parametrize("command,reason", [
    (OWNER["command"], "replacement_owner_command_present_not_old_owner_completion"),
    (["julia", "/unrelated/run_sia_geometry_job.jl"], "physical_producer_or_pair_reader_present"),
    (["time", "julia", "audit_sia_geometry_job.jl"], "physical_producer_or_pair_reader_present"),
    (["julia", "worker.jl", "--run=/synthetic/sidis-run/p3/report.json"], "bound_run_process_present_including_orphans"),
])
def test_replacement_or_orphan_cannot_release_host(tmp_path, command, reason):
    proc = fake_proc(tmp_path)
    add_process(proc, 456, 600, command, parent=1)
    result = host.observe_host(OWNER, ["/synthetic/sidis-run"], proc=proc)
    assert result["status"] == "HOST_RESOURCE_GUARD_BLOCKED"
    assert result["expected_owner_absent"] is True
    assert result["blockers"][0]["reason"] == reason
    assert "command" not in result["blockers"][0]


def test_owner_command_change_and_bad_observation_are_not_release(tmp_path):
    proc = fake_proc(tmp_path)
    path = add_process(proc, 123, 500, ["/changed/command"])
    result = host.observe_host(OWNER, proc=proc)
    assert result["status"] == "HOST_RESOURCE_GUARD_BLOCKED"
    assert result["expected_owner_live"] is False and result["expected_owner_absent"] is False
    (path / "stat").write_text("not a process identity")
    result = host.observe_host(OWNER, proc=proc)
    assert result["status"] == "HOST_RESOURCE_OBSERVATION_UNKNOWN_NOT_RELEASED"
    assert result["observation_errors"]


@pytest.mark.parametrize("key,value", [("pid", True), ("start_ticks", False), ("boot_id", "unknown"), ("command", [])])
def test_invalid_owner_tuple_rejected(tmp_path, key, value):
    owner = copy.deepcopy(OWNER)
    owner[key] = value
    with pytest.raises(ValueError):
        host.observe_host(owner, proc=fake_proc(tmp_path))


@pytest.mark.parametrize("root", ["/", "relative/run", "/synthetic//run", "/synthetic/run/", "/synthetic/../run", "//synthetic/run"])
def test_noncanonical_run_roots_rejected(tmp_path, root):
    with pytest.raises(ValueError):
        host.observe_host(OWNER, [root], proc=fake_proc(tmp_path))


def test_zero_tick_kernel_process_and_sibling_run_are_not_false_blockers(tmp_path):
    proc = fake_proc(tmp_path)
    add_process(proc, 2, 0, [], parent=0)
    add_process(proc, 999, 900, ["metadata-reader", "/synthetic/sidis-run-copy/report.json"])
    result = host.observe_host(OWNER, ["/synthetic/sidis-run"], proc=proc)
    assert result["status"] == "HOST_RESOURCE_GUARD_CLEAR_NOT_SCIENTIFIC_ADMISSION"


def test_missing_boot_is_unknown_not_release(tmp_path):
    result = host.observe_host(OWNER, proc=tmp_path / "missing-proc")
    assert result["status"] == "HOST_RESOURCE_OBSERVATION_UNKNOWN_NOT_RELEASED"
    assert result["expected_owner_absent"] is False and result["live_host_observed"] is False


@pytest.mark.parametrize("state,start", [("RS", 600), ("S", -1)])
def test_malformed_process_state_or_time_cannot_clear_host(tmp_path, state, start):
    proc = fake_proc(tmp_path)
    add_process(proc, 456, start, ["unrelated"], state=state)
    assert host.observe_host(OWNER, proc=proc)["status"] == "HOST_RESOURCE_OBSERVATION_UNKNOWN_NOT_RELEASED"


def test_empty_userspace_command_is_unknown(tmp_path):
    proc = fake_proc(tmp_path)
    add_process(proc, 456, 600, [], parent=1)
    result = host.observe_host(OWNER, proc=proc)
    assert result["status"] == "HOST_RESOURCE_OBSERVATION_UNKNOWN_NOT_RELEASED"


def test_owner_exiting_before_listing_cannot_create_false_release(tmp_path, monkeypatch):
    proc = fake_proc(tmp_path)
    path = add_process(proc, 123, 500, OWNER["command"])
    original = type(proc).iterdir

    def exit_at_listing(directory):
        if directory == proc:
            (path / "stat").unlink()
            (path / "cmdline").unlink()
            path.rmdir()
        return original(directory)

    monkeypatch.setattr(type(proc), "iterdir", exit_at_listing)
    result = host.observe_host(OWNER, proc=proc)
    assert result["status"] == "HOST_RESOURCE_GUARD_BLOCKED"
    assert result["expected_owner_absent"] is False


def test_private_lock_is_exclusive_and_released_after_exception(tmp_path):
    lock = tmp_path / "synthetic-resource.lock"
    with pytest.raises(RuntimeError, match="synthetic interruption"):
        with host._lock_at(lock) as lease:
            assert lease["pid"] == os.getpid() and lease["standard_host_lock"] is False
            assert lease["accepted_fitting_gate"] is False
            with pytest.raises(BlockingIOError):
                with host._lock_at(lock):
                    pass
            raise RuntimeError("synthetic interruption")
    with host._lock_at(lock) as lease:
        assert lease["inode"] == lock.stat().st_ino


def test_alias_and_shared_permission_lock_rejected(tmp_path):
    target = tmp_path / "target"
    target.write_text("explicitly synthetic administrative fixture")
    alias = tmp_path / "alias"
    alias.symlink_to(target)
    with pytest.raises(ValueError):
        with host._lock_at(alias):
            pass
    target.chmod(0o644)
    with pytest.raises(ValueError):
        with host._lock_at(target):
            pass
