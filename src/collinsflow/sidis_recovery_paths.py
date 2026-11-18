"""Owned canonical job locations; this view is never dispatch or physics admission."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
POLICY_SHA = "b9b704ab85b473021de062ef555a1979d2522a9760de1122ffb3fa350b4bbcf2"
POINT_PHASES = ("native123", "native4", "tailref123", "reference4")
GLOBAL_PHASES = ("native-select3", "native-selection", "source-math", "derivative-select3", "diagnostic")
FIELDS = ("phase", "index", "folder", "producer_root", "producer_manifest", "producer_manifest_sha256", "protocol")
PROTOCOLS = ("frozen_original", "continuation_v1", "recovery_v2")
FLAGS = ("production_ready", "physical_extraction", "historical_fit_ready", "experimental_data_admitted",
         "fit_performed", "accepted_prediction_diagnostic", "accepted_fitting_gate", "G2PERF_pass",
         "interval_coverage_verified", "all132_coverage_verified")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def policy():
    path = ROOT / "config/sidis-recovery-contract.json"
    if digest(path) != POLICY_SHA:
        raise ValueError("prospective recovery policy changed")
    data = json.loads(path.read_text(encoding="utf-8"))
    if (data["retained_jobs"]["point_phases"] != list(POINT_PHASES) or
            any(data.get(key) is not False for key in FLAGS)):
        raise ValueError("wrong or self-admitting recovery policy")
    return data


def canonical(path, *, directory):
    if not isinstance(path, (str, Path)) or (isinstance(path, str) and str(Path(path)) != path):
        raise ValueError("explicit canonical path required")
    value = Path(path)
    if (not value.is_absolute() or value.resolve() != value or value.is_symlink() or
            any(character in str(value) for character in ("\0", "\n", "\r")) or
            not (value.is_dir() if directory else value.is_file())):
        raise ValueError("ordinary canonical file or directory required")
    return value


def point_index(phase, index):
    if phase not in POINT_PHASES or type(index) is not int or not 1 <= index <= 132:
        raise ValueError("declared point phase and exact index1..132 required")


def within_directory(path, root):
    """Also detect filesystem aliases such as case variants on a Windows volume."""
    target = root.stat()
    identity = target.st_dev, target.st_ino
    return any((value.stat().st_dev, value.stat().st_ino) == identity for value in (path, *path.parents))


def _entry(root, row):
    if not isinstance(row, dict) or set(row) != set(FIELDS):
        raise ValueError("complete location schema required; no admission fields allowed")
    point_index(row["phase"], row["index"])
    if row["protocol"] not in PROTOCOLS:
        raise ValueError("unsupported producing protocol")
    producer = canonical(row["producer_root"], directory=True)
    folder = canonical(row["folder"], directory=True)
    manifest = canonical(row["producer_manifest"], directory=False)
    if folder != producer / row["phase"] / f"p{row['index']:04d}":
        raise ValueError("folder does not match its exact producing phase/index")
    protected = producer.parent if row["protocol"] == "frozen_original" else producer
    if within_directory(root, protected) or within_directory(protected, root):
        raise ValueError("current and protected producing roots must be disjoint")
    if row["protocol"] == "frozen_original":
        if row["phase"] != "native123" or row["index"] > 13 or manifest != producer.parent / "invocation.json":
            raise ValueError("only the declared original13 native123 folders belong to the original anchor")
    elif manifest != producer / ("continuation.json" if row["protocol"] == "continuation_v1" else "recovery.json"):
        raise ValueError("manifest does not belong to the declared producing root/protocol")
    sha = row["producer_manifest_sha256"]
    if not isinstance(sha, str) or not re.fullmatch("[0-9a-f]{64}", sha) or digest(manifest) != sha:
        raise ValueError("changed externally bound producing manifest")
    return copy.deepcopy(row)


class RecoveryRunPaths:
    """An owned path view. Complete metadata and independent arithmetic stay mandatory."""

    def __init__(self, root, entries):
        policy()
        self.root = canonical(root, directory=True)
        if not isinstance(entries, list):
            raise ValueError("ordered retained-job entries required")
        self.entries = [_entry(self.root, row) for row in entries]
        keys = [(POINT_PHASES.index(row["phase"]), row["index"]) for row in self.entries]
        if keys != sorted(set(keys)):
            raise ValueError("unique source-ordered retained phase/index entries required")
        producers = {}
        for row in self.entries:
            producer = Path(row["producer_root"])
            identity = row["protocol"], row["producer_manifest"], row["producer_manifest_sha256"]
            protected = producer.parent if row["protocol"] == "frozen_original" else producer
            if producer in producers:
                if producers[producer][0] != identity:
                    raise ValueError("one producing root cannot claim different protocols/manifests")
            else:
                if any(within_directory(protected, old[1]) or within_directory(old[1], protected) for old in producers.values()):
                    raise ValueError("protected producing roots must be disjoint")
                producers[producer] = identity, protected
        self.signature = self._signature()

    def _signature(self):
        payload = [str(self.root), [[row[key] for key in FIELDS] for row in self.entries]]
        return hashlib.sha256(json.dumps(payload, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode()).hexdigest()

    def valid(self):
        policy()
        if self._signature() != self.signature or canonical(self.root, directory=True) != self.root:
            raise ValueError("mutated owned recovery path view")
        for row in self.entries:
            _entry(self.root, row)
        return self

    def retained(self, phase, index):
        point_index(phase, index)
        self.valid()
        return next((copy.deepcopy(row) for row in self.entries if row["phase"] == phase and row["index"] == index), None)

    def current_path(self, *parts):
        self.valid()
        value = self.root.joinpath(*parts)
        if not value.is_relative_to(self.root) or value.resolve() != value:
            raise ValueError("current-run paths cannot escape or traverse aliases")
        return value

    def folder(self, phase, index=0):
        if phase in POINT_PHASES:
            row = self.retained(phase, index)
            return Path(row["folder"]) if row else self.current_path(phase, f"p{index:04d}")
        if phase not in GLOBAL_PHASES or type(index) is not int or index != 0:
            raise ValueError("declared fresh global phase with index0 required")
        return self.current_path(phase)

    def fresh_folder(self, phase, index=0):
        folder = self.folder(phase, index)
        if not folder.is_relative_to(self.root) or folder.exists() or folder.is_symlink():
            raise ValueError("never dispatch through retained, aliased or occupied folders")
        return folder

    def metadata(self):
        self.valid()
        return dict(scope="owned_whole_job_path_view_not_admission", root=str(self.root), retained_jobs=copy.deepcopy(self.entries),
                    policy_sha256=POLICY_SHA, numerical_dispatch_ready=False, complete_job_provenance_verified=False,
                    saved_numerical_arithmetic_verified=False, **dict.fromkeys(FLAGS, False))
