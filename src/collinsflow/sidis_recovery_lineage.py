"""Read immutable SIDIS ancestry; inventory evidence never replaces Julia review."""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]


def sibling(name):
    spec = importlib.util.spec_from_file_location(name + "_recovery_adapter", Path(__file__).with_name(name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


old = sibling("sidis_continuation")
paths = sibling("sidis_recovery_paths")
SCOPE = "same_attempt_reconstructed_full132_recovery_not_fitting"
PREDECESSOR = "b524304c005091a5717de33141fc842fa07fbee1"
DECLARATION = "a88e9ba2d0059796a3dadc6c7b86707f36ded616"
FIRST_MANIFEST = "8743c5052a7b950812f21db35c83ae2008c43f2d32dbd6cfd12c602808fea83a"
FIRST_LOSS = "2afc9d062c0124897b45ed7507d09603fc9cfb0c6b7c781bb9064caa8bf7d1b7"
LAUNCH_SHA = "407642f86461642743a21d5b096decb4b850e25bafb185c7797ade6a558bc964"
PREPARED = "PREPARED_IMMUTABLE_RECOVERY_NOT_A_GATE"
COMPLETED = "COMPLETE_JOB_RECORDED_NOT_PRODUCTION_GATE"
ANCHOR_FIELDS = ("original", "original_inventory", "review", "integrity", "loss", "request", "request_sha256",
                 "native", "parton", "dependency_receipts", "numerical_source_fingerprint", "frozen", "frozen_head", "frozen_files")


def sha_file(path, expected):
    value = paths.canonical(path, directory=False)
    if not isinstance(expected, str) or not re.fullmatch("[0-9a-f]{64}", expected) or paths.digest(value) != expected:
        raise ValueError("changed externally bound immutable file")
    return value


def inventory(root):
    """Inventory links as links; never follow an ancestor's saved aliases."""
    root = paths.canonical(root, directory=True)
    result = []
    for current, directories, names in os.walk(root, followlinks=False):
        for name in sorted(directories + names):
            path = Path(current) / name
            if path.is_symlink():
                result.append(dict(path=str(path), kind="symlink_retained_as_link_not_followed", target=os.readlink(path)))
                if name in directories:
                    directories.remove(name)
            elif path.is_file():
                if path.resolve() != path:
                    raise ValueError("ordinary ancestor files cannot traverse aliases")
                result.append(dict(path=str(path), kind="ordinary_file", bytes=path.stat().st_size, sha256=paths.digest(path)))
            elif not path.is_dir():
                raise ValueError("nonordinary ancestor output")
    return sorted(result, key=lambda row: row["path"])


def verify_inventory(root, rows, workspace):
    root = paths.canonical(root, directory=True)
    expected = []
    for row in rows:
        value = Path(row["path"])
        value = value if value.is_absolute() else Path(workspace) / value
        if str(value) != os.path.normpath(str(value)) or not value.is_relative_to(root) or value == root:
            raise ValueError("ancestor inventory path outside its canonical root")
        kind = row.get("kind", "ordinary_file")
        item = dict(path=str(value), kind=kind)
        if kind == "ordinary_file":
            if type(row["bytes"]) is not int or row["bytes"] < 0:
                raise ValueError("exact nonnegative inventory size required")
            item.update(bytes=row["bytes"], sha256=row["sha256"])
        elif kind == "symlink_retained_as_link_not_followed":
            item["target"] = row["target"]
        else:
            raise ValueError("unknown immutable inventory entry")
        expected.append(item)
    expected.sort(key=lambda row: row["path"])
    if len({row["path"] for row in expected}) != len(expected) or inventory(root) != expected:
        raise ValueError("ancestor namespace inventory changed")
    return expected


def complete_locations(root, protocol, manifest, manifest_sha):
    """Select complete point pairs only; terminal failures remain failures."""
    result = []
    for phase in paths.POINT_PHASES:
        phase_root = Path(root) / phase
        if phase_root.is_symlink() or (phase_root.exists() and phase_root.resolve() != phase_root):
            raise ValueError("aliased ancestor phase")
        for folder in sorted(phase_root.glob("*")):
            if not re.fullmatch(r"p\d{4}", folder.name):
                raise ValueError("undeclared ancestor point folder")
            index = int(folder.name[1:])
            paths.point_index(phase, index)
            paths.canonical(folder, directory=True)
            invocation = folder / "invocation.json"
            if not invocation.is_file():
                continue
            inv = old.read_json(invocation)
            if inv.get("status") != COMPLETED:
                if inv.get("completed_at") is not None or (folder / "job.completed").exists():
                    raise ValueError("terminal failed point cannot become a recovery candidate")
                continue
            if (not inv.get("completed_at") or inv.get("phase") != phase or type(inv.get("point_index")) is not int or
                    inv["point_index"] != index or not (folder / "job.completed").is_file()):
                raise ValueError("false or incomplete point completion")
            if [s["name"] for s in inv.get("steps", [])] != ["producer", "audit"]:
                raise ValueError("complete producer/auditor pair required")
            for step in inv["steps"]:
                if (type(step.get("returncode")) is not int or step["returncode"] != 0 or
                        type(step.get("process_returncode")) is not int or step["process_returncode"] != 0 or
                        step.get("interruption") is not None or step.get("termination_signal") is not None):
                    raise ValueError("timeout, signal or effective/raw failure cannot approve recovery")
            result.append(dict(phase=phase, index=index, folder=str(folder), producer_root=str(root),
                producer_manifest=str(manifest), producer_manifest_sha256=manifest_sha, protocol=protocol))
    for phase in paths.GLOBAL_PHASES:
        invocation = Path(root) / phase / "invocation.json"
        if invocation.is_file():
            inv = old.read_json(invocation)
            if inv.get("status") != COMPLETED and inv.get("completed_at") is not None:
                raise ValueError("failed terminal global gate cannot be silently retried")
    return result


def ordered(entries):
    result = sorted(entries, key=lambda row: (paths.POINT_PHASES.index(row["phase"]), row["index"]))
    if len({(row["phase"], row["index"]) for row in result}) != len(result):
        raise ValueError("conflicting retained phase/index origins")
    for phase in paths.POINT_PHASES:
        indices = [row["index"] for row in result if row["phase"] == phase]
        if indices != list(range(1, len(indices) + 1)):
            raise ValueError("retained phase must be a complete source-ordered prefix")
    return result


def original_locations(first):
    root = Path(first["original"]) / "full"
    return [dict(phase="native123", index=index, folder=str(root / "native123" / f"p{index:04d}"),
        producer_root=str(root), producer_manifest=str(root.parent / "invocation.json"),
        producer_manifest_sha256=paths.policy()["logical_attempt"]["original_outer_invocation_sha256"],
        protocol="frozen_original") for index in range(1, 14)]


def source_bindings(frozen, predecessor, adapter, m=None):
    frozen, predecessor, adapter = (paths.canonical(value, directory=True) for value in (frozen, predecessor, adapter))
    fhead, ffiles = old.verify_git(frozen, old.FROZEN_HEAD)
    phead, pfiles = old.verify_git(predecessor, PREDECESSOR)
    ahead, afiles = old.verify_git(adapter, None if m is None else m["adapter_commit"])
    if (len(ffiles) != 335 or len(pfiles) != 348 or any(pfiles.get(k) != v for k, v in ffiles.items()) or
            any(afiles.get(k) != v for k, v in pfiles.items())):
        raise ValueError("all335 frozen and348 predecessor files must remain unchanged")
    import subprocess
    if subprocess.run(["git", "merge-base", "--is-ancestor", DECLARATION, ahead], cwd=adapter, check=False).returncode:
        raise ValueError("recovery implementation must descend from its prospective declaration")
    if m is not None and (m["frozen"] != str(frozen) or m["predecessor"] != str(predecessor) or
            m["adapter"] != str(adapter) or m["frozen_files"] != ffiles or m["predecessor_files"] != pfiles or
            m["adapter_files"] != afiles or m["frozen_head"] != fhead or m["predecessor_head"] != phead):
        raise ValueError("changed explicit source bindings")
    return dict(frozen=str(frozen), frozen_head=fhead, frozen_files=ffiles, predecessor=str(predecessor),
        predecessor_head=phead, predecessor_files=pfiles, adapter=str(adapter), adapter_commit=ahead, adapter_files=afiles)


def validate_ancestors(nodes, frozen, predecessor, current_root):
    """Externally bound files select data, never the code that executes this reader."""
    if not isinstance(nodes, list) or not 1 <= len(nodes) <= 64:
        raise ValueError("one to64 ordered ancestor nodes required")
    workspace = Path(frozen).parents[2]
    seen, data, entries = [], [], []
    for position, node in enumerate(nodes):
        if set(node) != {"protocol", "manifest", "manifest_sha256", "inventory", "inventory_sha256"}:
            raise ValueError("complete externally bound ancestry schema required")
        manifest = sha_file(node["manifest"], node["manifest_sha256"])
        evidence = sha_file(node["inventory"], node["inventory_sha256"])
        root = paths.canonical(manifest.parent, directory=True)
        if (not root.is_relative_to(Path(frozen) / "runs") or
                paths.within_directory(root, current_root) or paths.within_directory(current_root, root) or
                any(paths.within_directory(root, parent) or paths.within_directory(parent, root) for parent in seen)):
            raise ValueError("repeated, nested, aliased or nonfamily ancestor root")
        parent = old.read_json(manifest)
        loss = old.read_json(evidence)
        if position == 0:
            if (node["protocol"] != "continuation_v1" or node["manifest_sha256"] != FIRST_MANIFEST or
                    node["inventory_sha256"] != FIRST_LOSS or manifest.name != "continuation.json" or
                    parent["adapter"] != str(predecessor) or parent["adapter_commit"] != PREDECESSOR or
                    parent["contract_sha256"] != old.CONTRACT_SHA or parent["scope"] != old.SCOPE or
                    parent["status"] != "PREPARED_IMMUTABLE_CONTINUATION_NOT_A_GATE"):
                raise ValueError("the exact supported stopped continuation is the first ancestor")
            first = parent
            entries = original_locations(first)
            if (parent["original"] != str(Path(frozen) / "runs/historical-sidis-jobs-real-data/attempt-2") or
                    paths.within_directory(root, Path(first["original"])) or
                    paths.within_directory(current_root, Path(first["original"])) or
                    paths.within_directory(Path(first["original"]), current_root)):
                raise ValueError("wrong or overlapping original attempt2 anchor")
            if (loss.get("old_owner_absent_verified") is not True or loss.get("source_run") != str(root.relative_to(workspace)) or
                    loss.get("recorded_complete_prefix") != 37 or len(loss["retained_files"]) != 404):
                raise ValueError("wrong supported interruption inventory")
        else:
            if (node["protocol"] != "recovery_v2" or manifest.name != "recovery.json" or
                    parent.get("status") != PREPARED or parent.get("scope") != SCOPE or
                    parent.get("contract_sha256") != paths.POLICY_SHA or parent.get("ancestors") != nodes[:position] or
                    parent.get("retained_jobs") != ordered(entries) or
                    loss.get("status") != "STOPPED_SIDIS_NAMESPACE_INVENTORIED_NOT_A_GATE" or
                    loss.get("manifest_sha256") != node["manifest_sha256"] or loss.get("run") != str(root) or
                    loss.get("old_owner_absent_verified") is not True):
                raise ValueError("inconsistent, cyclic or unreviewed recovery ancestry")
            # An ancestor adapter is checked as immutable evidence. It is never imported.
            _, ancestor_files = old.verify_git(parent["adapter"], parent["adapter_commit"])
            if ancestor_files != parent["adapter_files"] or any(ancestor_files.get(k) != v for k, v in first["adapter_files"].items()):
                raise ValueError("changed ancestor adapter source")
        if (parent["run"] != str(root) or any(parent.get(k) != first[k] for k in ANCHOR_FIELDS) or
                any(parent.get(k) is not False for k in paths.FLAGS)):
            raise ValueError("ancestor belongs to a different or self-admitting numerical family")
        verify_inventory(root, loss["retained_files"], workspace)
        entries += complete_locations(root, node["protocol"], manifest, node["manifest_sha256"])
        entries = ordered(entries)
        seen.append(root)
        data.append(parent)
    if len([row for row in entries if row["protocol"] == "continuation_v1"]) != 24:
        raise ValueError("supported predecessor must retain exactly its24 complete point pairs")
    view = paths.RecoveryRunPaths(current_root, entries)
    return dict(first=first, nodes=data, entries=entries, view=view)


def verify_original(first, frozen, api):
    policy = paths.policy()
    anchor = policy["logical_attempt"]
    for key, expected in (("review", anchor["independent_original13_review_sha256"]),
                          ("integrity", anchor["post_review_original141_integrity_sha256"]),
                          ("loss", anchor["original_loss_inventory_sha256"])):
        sha_file(first[key], expected)
    loss = old.read_json(first["loss"])
    workspace = Path(frozen).parents[2]
    rows = [dict(row, path=str(workspace / row["path"])) for row in loss["interrupted_attempts"]["sidis132"]["retained_files"]]
    if len(rows) != 141 or rows != first["original_inventory"]:
        raise ValueError("changed original141-file anchor")
    old.verify_inventory(first["original"], rows)
    sha_file(Path(first["original"]) / "invocation.json", anchor["original_outer_invocation_sha256"])
    if (api.ROOT != Path(frozen) or api.source_fingerprint(frozen) != policy["frozen_numerics"]["source_fingerprint"] or
            first["numerical_source_fingerprint"] != policy["frozen_numerics"]["source_fingerprint"] or
            first["request"] != str(Path(first["original"]) / "full/request.json")):
        raise ValueError("wrong frozen library or source request")
    sha_file(first["request"], first["request_sha256"])
    if old.read_json(Path(first["original"]) / "invocation.json")["dependency_receipts"] != first["dependency_receipts"]:
        raise ValueError("changed original build identities")
    if set(first["dependency_receipts"]) != {first["native"], first["parton"]}:
        raise ValueError("exact native and parton build receipts required")
    for name, sha in first["dependency_receipts"].items():
        sha_file(name, sha)
    request = old.read_json(first["request"])
    if (request["source_fingerprint"] != first["numerical_source_fingerprint"] or len(request["points"]) != 132 or
            len({point["id"] for point in request["points"]}) != 132):
        raise ValueError("full ordered132-point request required")
    for name, sha in request["source_inputs"].items():
        sha_file(name, sha)
    for index in range(1, 14):
        old._verify_original_job(first, index, request, api)
        paths.canonical(old.point_folder(first, "native123", index) / "job.completed", directory=False)
    return request
