"""Observe physical-job ownership and hold one Linux geometry resource lease."""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import os
from pathlib import Path, PurePosixPath
import stat
import sys
import uuid

PHYSICAL_PROGRAMS = frozenset((
    "run_sia_geometry_job.jl", "audit_sia_geometry_job.jl",
    "predict_historical_sia_integrals.jl", "predict_historical_sia_points.jl",
    "predict_historical_sidis_points.jl", "run_historical_sidis_point_job.jl",
    "run_historical_sidis_global_job.jl",
    "run_sidis_continuation_point.jl", "run_sidis_continuation_global.jl",
    "audit_sidis_continuation_job.jl", "audit_sidis_continuation_prefix.jl",
    "run_sidis_recovery_point.jl", "run_sidis_recovery_global.jl",
    "audit_sidis_recovery_job.jl", "audit_sidis_recovery_retained.jl",
    "run_sidis_continuation.py", "run_sidis_recovery.py", "run_sia_geometry_jobs.py",
))


def process_stat(text):
    prefix, tail = text.rsplit(")", 1)
    fields = tail.split()
    pid = int(prefix.split("(", 1)[0])
    result = dict(pid=pid, state=fields[0], parent_pid=int(fields[1]), start_ticks=int(fields[19]))
    if (pid <= 0 or result["parent_pid"] < 0 or result["start_ticks"] < 0 or
            len(result["state"]) != 1 or result["state"] not in "RSDZTtWXxKPI"):
        raise ValueError("invalid process identity")
    return result


def _owner(owner):
    if (set(owner) != {"pid", "start_ticks", "boot_id", "command"} or
            any(type(owner[key]) is not int or owner[key] <= 0 for key in ("pid", "start_ticks")) or
            not isinstance(owner["boot_id"], str) or str(uuid.UUID(owner["boot_id"])) != owner["boot_id"] or
            not isinstance(owner["command"], list) or not owner["command"] or
            any(not isinstance(value, str) or not value or "\0" in value for value in owner["command"])):
        raise ValueError("explicit PID/start/boot/command owner tuple required")


def _run_root(value):
    if (not isinstance(value, str) or not value.startswith("/") or value.startswith("//") or str(PurePosixPath(value)) != value or
            value == "/" or ".." in PurePosixPath(value).parts or "\0" in value):
        raise ValueError("explicit canonical run roots required")
    return value


def _mentions_run(argument, roots):
    value = argument.split("=", 1)[1] if argument.startswith("--") and "=" in argument else argument
    return any(value == root or value.startswith(root + "/") for root in roots)


def observe_host(owner, run_roots=(), *, proc=Path("/proc")):
    """A missing owner is not scientific completion; orphan processes still block."""
    _owner(owner)
    roots = tuple(_run_root(value) for value in run_roots)
    proc = Path(proc)
    live = proc == Path("/proc") and proc.resolve() == proc
    result = dict(status="HOST_RESOURCE_OBSERVATION_UNKNOWN_NOT_RELEASED", observed_at=datetime.now(timezone.utc).isoformat(),
        theory_scope="resource_observation_not_physics", live_host_observed=live,
        input_scope="live_Linux_proc" if live else "explicitly_synthetic_proc_fixture_not_live_host",
        blockers=[], observation_errors=[], expected_owner_live=False, expected_owner_absent=False,
        jobs_started_stopped_or_restarted=False, fit_performed=False, production_ready=False, accepted_fitting_gate=False)
    try:
        boot = (proc / "sys/kernel/random/boot_id").read_text(encoding="utf-8").strip()
        if str(uuid.UUID(boot)) != boot:
            raise ValueError("noncanonical boot identity")
        result["boot_id"] = boot
    except (OSError, ValueError) as error:
        result["observation_errors"].append(dict(error_type=type(error).__name__, scope="process_directory_or_boot"))
        return result
    owner_seen = False
    # Retain an owner seen before listing processes. If it exits while the
    # directory is scanned, this observation cannot miss its last fork and clear.
    if boot == owner["boot_id"]:
        try:
            path = proc / str(owner["pid"])
            identity = process_stat((path / "stat").read_text(encoding="utf-8"))
            if identity["pid"] == owner["pid"] and identity["start_ticks"] == owner["start_ticks"] and identity["state"] != "Z":
                owner_seen = True
                result["expected_owner_live"] = (path / "cmdline").read_bytes() == ("\0".join(owner["command"]) + "\0").encode()
                result["blockers"].append({**identity, "reason": "expected_owner_live" if result["expected_owner_live"] else
                    "expected_owner_command_mismatch_not_release"})
        except FileNotFoundError:
            if owner_seen:
                result["blockers"].append({**identity, "reason": "owner_seen_before_directory_scan_not_release"})
        except (OSError, ValueError, UnicodeError, IndexError) as error:
            result["observation_errors"].append(dict(pid=owner["pid"], error_type=type(error).__name__))
    try:
        entries = list(proc.iterdir())
    except OSError as error:
        result["observation_errors"].append(dict(error_type=type(error).__name__, scope="process_directory"))
        return result
    for path in sorted((p for p in entries if p.name.isdigit()), key=lambda p: int(p.name)):
        if int(path.name) == os.getpid():
            continue
        try:
            identity = process_stat((path / "stat").read_text(encoding="utf-8"))
            if identity["pid"] != int(path.name):
                raise ValueError("process directory and stat disagree")
            if identity["state"] == "Z":
                continue
            command = (path / "cmdline").read_bytes().decode("utf-8").rstrip("\0").split("\0")
            command = [] if command == [""] else command
            exact_process = (boot == owner["boot_id"] and identity["pid"] == owner["pid"] and identity["start_ticks"] == owner["start_ticks"])
            reason = None
            if exact_process:
                owner_seen = True
                result["expected_owner_live"] = command == owner["command"]
                reason = "expected_owner_live" if result["expected_owner_live"] else "expected_owner_command_mismatch_not_release"
            elif command == owner["command"]:
                reason = "replacement_owner_command_present_not_old_owner_completion"
            elif any(PurePosixPath(argument).name in PHYSICAL_PROGRAMS for argument in command):
                reason = "physical_producer_or_pair_reader_present"
            elif any(_mentions_run(argument, roots) for argument in command):
                reason = "bound_run_process_present_including_orphans"
            if reason:
                result["blockers"].append({**identity, "reason": reason})
            elif not command and identity["parent_pid"] not in (0, 2):
                result["observation_errors"].append(dict(pid=identity["pid"], error_type="empty_userspace_command_unknown"))
        except FileNotFoundError:
            continue
        except (OSError, ValueError, UnicodeError, IndexError) as error:
            result["observation_errors"].append(dict(pid=int(path.name), error_type=type(error).__name__))
    result["expected_owner_absent"] = not owner_seen
    result["blockers"] = list({(row["pid"], row["start_ticks"]): row for row in result["blockers"]}.values())
    if result["observation_errors"]:
        return result
    result["status"] = "HOST_RESOURCE_GUARD_BLOCKED" if result["blockers"] else "HOST_RESOURCE_GUARD_CLEAR_NOT_SCIENTIFIC_ADMISSION"
    return result


def standard_lock_path():
    if not sys.platform.startswith("linux"):
        raise ValueError("physical geometry resource leases require Linux")
    return Path("/tmp") / f"collinsflow-sia-geometry-{os.getuid()}.lock"


@contextmanager
def _lock_at(path):
    """Private test seam; production uses the one fixed per-user host lock."""
    import fcntl

    path = Path(path)
    if not path.is_absolute() or path.parent.resolve() != path.parent or path.resolve() != path:
        raise ValueError("canonical resource lock path required")
    descriptor = os.open(path, os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW | os.O_CLOEXEC, 0o600)
    try:
        info = os.fstat(descriptor)
        if not stat.S_ISREG(info.st_mode) or info.st_uid != os.getuid() or info.st_nlink != 1 or info.st_mode & 0o077:
            raise ValueError("ordinary private single-link resource lock required")
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        identity = process_stat(Path("/proc/self/stat").read_text(encoding="utf-8"))
        lease = dict(scope="geometry_resource_lease_not_scientific_admission", lock_path=str(path),
            standard_host_lock=path == standard_lock_path(), file_descriptor=descriptor, device=info.st_dev, inode=info.st_ino,
            pid=identity["pid"], start_ticks=identity["start_ticks"],
            boot_id=Path("/proc/sys/kernel/random/boot_id").read_text(encoding="utf-8").strip(),
            fit_performed=False, production_ready=False, accepted_fitting_gate=False)
        yield lease
    finally:
        os.close(descriptor)


def geometry_lock():
    """Hold across both Julia commands; process scans must precede dispatch."""
    return _lock_at(standard_lock_path())
