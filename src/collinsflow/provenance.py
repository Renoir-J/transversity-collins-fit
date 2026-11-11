"""Content-addressed artifacts and receipts; no executable serialization."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as stream:
            json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def computational_sources(root):
    root = Path(root)
    files = []
    for directory in ("src", "numerics/src", "numerics/bin", "fortran", "workflow", "config", "fixtures"):
        files.extend(p for p in (root / directory).rglob("*") if p.is_file()
                     and not any(part == "__pycache__" or part.endswith(".egg-info") for part in p.parts)
                     and p.suffix not in {".so", ".dll", ".dylib", ".o", ".mod"}
                     and p.name not in {"github-state.json", "project-tasks.json", "actions.lock.json"})
    for name in ("pyproject.toml", "requirements.lock", "numerics/Project.toml", "numerics/Manifest.toml",
                 "references/source-inventory.json", "references/code-review.json"):
        if (root / name).is_file():
            files.append(root / name)
    return sorted(set(files))


def source_fingerprint(root):
    root = Path(root)
    return digest({p.relative_to(root).as_posix(): file_hash(p) for p in computational_sources(root)})


def git_commit(root):
    result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True)
    return result.stdout.strip() if result.returncode == 0 else "UNCOMMITTED"


def environment(root, julia):
    version = subprocess.run([julia, "--startup-file=no", "--version"], capture_output=True,
                             text=True, check=True, timeout=60).stdout.strip()
    details = {"python": platform.python_version(), "platform": platform.platform(), "julia": version,
               # Numerical workers are explicitly launched with these values; the coordinator may differ.
               "threads": {key: "1" for key in
                           ("JULIA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")}}
    for name in ("requirements.lock", "numerics/Manifest.toml"):
        details[name] = file_hash(Path(root) / name) if (Path(root) / name).exists() else "MISSING"
    return details


def verify_receipt(path, expected_scope, fingerprint, environment_digest=None, _active=None):
    active = set() if _active is None else set(_active)
    resolved = str(Path(path).resolve())
    if resolved in active:
        raise ValueError("cyclic receipt dependencies")
    active.add(resolved)
    receipt = read_json(path)
    if receipt["status"] != "PASS" or receipt["scope"] != expected_scope:
        raise ValueError(f"non-PASS or wrong-scope receipt: {path}")
    if receipt["source_fingerprint"] != fingerprint:
        raise ValueError(f"stale code/configuration receipt: {path}")
    if environment_digest and receipt["environment_digest"] != environment_digest:
        raise ValueError(f"stale environment receipt: {path}")
    for filename, expected in receipt["input_hashes"].items():
        if not Path(filename).is_file() or file_hash(filename) != expected:
            raise ValueError(f"changed receipt input: {filename}")
        if Path(filename).parent.name == "receipts" and Path(filename).suffix == ".json":
            verify_receipt(filename, expected_scope, fingerprint, environment_digest, active)
    for output in receipt["outputs"]:
        if not Path(output["path"]).is_file() or file_hash(output["path"]) != output["sha256"]:
            raise ValueError(f"changed receipt output: {output['path']}")
    if not receipt["checks"] or not all(c["passed"] for c in receipt["checks"]):
        raise ValueError(f"missing or failed checks: {path}")
    return receipt

