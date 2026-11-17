"""Source-bound planning for complete SIA jobs; no numerical dispatch or selection."""
from __future__ import annotations

import hashlib
import importlib
import json
import math
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone
import uuid

ROOT = Path(__file__).resolve().parents[2]
FROZEN_HEAD = "6c88d12516074889f49a0dfeae456d6cc63491b4"
FROZEN_FP = "01be55cc1449679f6dc26bd43ea9788604938c1fecbb25d276acaa35bc4bae47"
POLICY_HEAD = "655eee588b7bf8c69e120f3bb62fb707b794dfa2"
POLICY_SHA = "da14c5b667c2866aea2aebe546baedb45d97ad7d513c6bed1b1aa884a23bf337"
PARENT_SHA = "29830e78c4e48d2eb8e5a70711885fc961b2b94e0ec3303315befd59a65fdf97"
FLAGS = ("production_ready", "physical_extraction", "historical_fit_ready", "experimental_data_admitted",
         "fit_performed", "G2PERF_pass", "interval_coverage_verified")
REQUEST_SCOPE = "historical_reconstructed_sia_integral_prediction_diagnostic"
SCOPE = "historical_reconstructed_sia_complete_geometry_jobs"
PLAN_STATUS = "SOURCE_PLAN_PREPARED_NOT_NUMERICAL_DISPATCH_OR_GATE"


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], stderr=subprocess.PIPE)


def verify_checkout(root, *, head=None, ancestor=None):
    """Read raw blobs even when the index hides edits with assume-unchanged."""
    root = Path(root).resolve()
    current = git(root, "rev-parse", "HEAD").decode().strip()
    if head is not None and current != head:
        raise ValueError("wrong frozen source commit")
    if ancestor is not None:
        git(root, "merge-base", "--is-ancestor", ancestor, current)
    inventory = {}
    for entry in git(root, "ls-tree", "-r", "-z", current).split(b"\0"):
        if not entry:
            continue
        metadata, relative = entry.split(b"\t", 1)
        mode, kind, blob = metadata.split()
        if mode not in (b"100644", b"100755") or kind != b"blob":
            raise ValueError("only ordinary committed source files are allowed")
        relative = relative.decode()
        path = root / relative
        if path.resolve() != path or not path.is_file():
            raise ValueError("source aliases or missing files are not allowed")
        contents = path.read_bytes()
        actual = hashlib.sha1(b"blob " + str(len(contents)).encode() + b"\0" + contents).hexdigest()
        if actual != blob.decode():
            raise ValueError("working source differs from its raw committed blob: " + relative)
        inventory[relative] = hashlib.sha256(contents).hexdigest()
    if git(root, "ls-files", "--others", "--exclude-standard", "-z").strip(b"\0"):
        raise ValueError("commit the complete additive adapter before source planning")
    return current, inventory


def load_policy(adapter=ROOT):
    path = Path(adapter) / "config/historical-sia-job-contract.json"
    if file_hash(path) != POLICY_SHA:
        raise ValueError("prospective SIA job policy changed")
    policy = read_json(path)
    if policy["scope"] != SCOPE or any(policy.get(key) is not False for key in (*FLAGS, "accepted_fitting_gate")):
        raise ValueError("job policy cannot admit fitting or production")
    return policy


def geometry_partition(request, policy):
    """Map source identities only; source replay must precede physical use."""
    if (request.get("scope") != REQUEST_SCOPE or request.get("theory_scope") != "reconstructed_theory" or
            request.get("source_fingerprint") != FROZEN_FP or any(request.get(key) is not False for key in FLAGS)):
        raise ValueError("wrong or self-admitting full SIA source request")
    points = request["points"]
    expected = [f"KPSY15_SIA_{view}:r{row}" for view in (1002, 1003) for row in range(2, 38)]
    if len(points) != 72 or [point["id"] for point in points] != expected:
        raise ValueError("full ordered72 source rows are required")
    groups = {}
    derivative_ids = policy["independent_and_derivative_checks"]["representative_source_ids"]
    for index, point in enumerate(points, 1):
        ratio = "UC" if index <= 36 else "UL"
        values = [point["kinematics"][name] for name in ("z1", "z2", "Q2_GeV2")]
        if (point["ratio"] != ratio or any(type(value) is not float or not math.isfinite(value) for value in values) or
                not all(0 < value < 1 for value in values[:2]) or values[2] != 112.36):
            raise ValueError("changed literal source geometry or ratio")
        angular = point["source_angular"]
        if (angular["average_type"] != "mean_of_angular_ratio" or angular["frame"] != "hadron2" or
                angular["measure_status"] != "joint_acceptance_unresolved" or
                type(angular["value"]) is not float or not 0 <= angular["value"] <= 1):
            raise ValueError("unreviewed angular source")
        key = tuple(value.hex() for value in values)
        if key not in groups:
            groups[key] = {"index": len(groups) + 1, "geometry": values, "Float64_hex": list(key), "points": []}
        groups[key]["points"].append({"request_index": index, "id": point["id"], "ratio": ratio,
                                     "source_point_sha256": canonical_hash(point),
                                     "angular_source": json.loads(json.dumps(angular)),
                                     "representative_derivative": point["id"] in derivative_ids})
    if len(groups) != 36 or any([p["ratio"] for p in row["points"]] != ["UC", "UL"] for row in groups.values()):
        raise ValueError("exact36 geometries each require one source UC and one source UL")
    selected = [row["index"] for row in groups.values() if any(p["representative_derivative"] for p in row["points"])]
    if len(selected) != 3:
        raise ValueError("all three prospective derivative geometries are required")
    return list(groups.values()), selected


def _trusted_api(frozen):
    sys.path.insert(0, str(frozen / "src"))
    api = importlib.import_module("collinsflow.historical_sia_integrals")
    provenance = importlib.import_module("collinsflow.provenance")
    if any(Path(module.__file__).resolve().parents[2] != frozen for module in (api, provenance)):
        raise ValueError("source replay must import the explicitly verified frozen library")
    return api, provenance


def prepare_plan(frozen, request_path, output, *, adapter=ROOT):
    """Reconstruct a full source request below a fresh run, without Julia jobs."""
    frozen, adapter, request_path = (Path(path).resolve() for path in (frozen, adapter, request_path))
    folder = Path(output).absolute()
    if folder.resolve() != folder or not folder.is_relative_to(adapter / "runs") or folder.exists():
        raise ValueError("use a fresh canonical run folder below the adapter runs directory")
    policy = load_policy(adapter)
    old_head, old_files = verify_checkout(frozen, head=FROZEN_HEAD)
    new_head, new_files = verify_checkout(adapter, ancestor=POLICY_HEAD)
    if new_head == POLICY_HEAD or len(old_files) != 302 or any(new_files.get(name) != sha for name, sha in old_files.items()):
        raise ValueError("an additive committed implementation over all302 unchanged files is required")
    if file_hash(frozen / "config/historical-sia-integral-contract.json") != PARENT_SHA:
        raise ValueError("parent scientific contract changed")
    api, provenance = _trusted_api(frozen)
    if provenance.source_fingerprint(frozen) != FROZEN_FP:
        raise ValueError("frozen numerical fingerprint changed")
    folder.mkdir(parents=True, exist_ok=False)
    invocation = {"status": "FAIL", "record_role": "source_planning_not_numerical_job",
                  "started_at": datetime.now(timezone.utc).isoformat(), **dict.fromkeys(FLAGS, False),
                  "accepted_fitting_gate": False, "numerical_dispatch_ready": False}
    provenance.atomic_json(folder / "invocation.json", invocation)
    try:
        request_sha = file_hash(request_path)
        original = read_json(request_path)
        source = original["source_preparation"]
        fresh = folder / "request.json"
        # verify_request creates siblings of its input. Use the explicit-output
        # preparation instead; its historical-intake verifier is read only.
        api.prepare_request(source["manifest"], source["cache"], source["candidates"], fresh)
        request = read_json(fresh)
        if canonical_hash(request) != canonical_hash(original) or file_hash(request_path) != request_sha:
            raise ValueError("source-derived full request differs from the supplied request")
        if any(file_hash(path) != sha for path, sha in request["source_inputs"].items()):
            raise ValueError("changed source input")
        geometries, derivative_indices = geometry_partition(request, policy)
        plan = {"schema_version": 1, "status": PLAN_STATUS, "scope": SCOPE,
                "theory_scope": "reconstructed_theory", **dict.fromkeys(FLAGS, False),
                "accepted_fitting_gate": False, "numerical_dispatch_ready": False,
                "physical_jobs_started": False, "policy_commit": POLICY_HEAD, "policy_sha256": POLICY_SHA,
                "frozen": str(frozen), "frozen_head": old_head, "numerical_source_fingerprint": FROZEN_FP,
                "frozen_files": old_files, "adapter": str(adapter), "adapter_head": new_head,
                "adapter_files": new_files, "request": str(fresh), "request_sha256": file_hash(fresh),
                "source_request": str(request_path), "source_request_sha256": request_sha, "source_inputs": request["source_inputs"],
                "point_order": [point["id"] for point in request["points"]], "geometries": geometries,
                "derivative_geometry_indices": derivative_indices, "proposed_phases": policy["job_phases"],
                "resources": policy["resources"], "limitations": request["limitations"],
                "implementation_remaining": ["complete numerical producers and auditors", "global selectors and complete ledger",
                                             "fresh cached physical diagnostic and products"]}
        if (verify_checkout(frozen, head=old_head)[1] != old_files or
                verify_checkout(adapter, head=new_head)[1] != new_files or
                any(file_hash(path) != sha for path, sha in request["source_inputs"].items())):
            raise ValueError("source identities changed during planning")
        provenance.atomic_json(folder / "plan.json", plan)
        invocation.update(status=PLAN_STATUS, plan_sha256=file_hash(folder / "plan.json"), geometries=36, points=72)
        return plan
    except BaseException as error:
        invocation["error"] = type(error).__name__ + ": " + str(error)
        raise
    finally:
        invocation["completed_at"] = datetime.now(timezone.utc).isoformat()
        provenance.atomic_json(folder / "invocation.json", invocation)


def verify_plan(path, *, expected_sha256):
    """Replay source metadata below the plan folder; never authorize numerics."""
    path = Path(path).absolute()
    if path.resolve() != path or file_hash(path) != expected_sha256:
        raise ValueError("changed externally bound source plan")
    plan = read_json(path)
    if (type(plan.get("schema_version")) is not int or plan["schema_version"] != 1 or
            plan.get("status") != PLAN_STATUS or plan.get("scope") != SCOPE or
            plan.get("theory_scope") != "reconstructed_theory" or
            any(plan.get(key) is not False for key in (*FLAGS, "accepted_fitting_gate", "numerical_dispatch_ready", "physical_jobs_started"))):
        raise ValueError("source plan cannot admit numerical jobs or fitting")
    frozen, adapter = Path(plan["frozen"]), Path(plan["adapter"])
    if adapter != ROOT or path.parent.resolve() != path.parent or not path.parent.is_relative_to(adapter / "runs"):
        raise ValueError("wrong canonical source-plan namespace")
    policy = load_policy(adapter)
    old_head, old_files = verify_checkout(frozen, head=FROZEN_HEAD)
    new_head, new_files = verify_checkout(adapter, head=plan["adapter_head"], ancestor=POLICY_HEAD)
    if (new_head == POLICY_HEAD or len(old_files) != 302 or plan["frozen_head"] != old_head or
            plan["frozen_files"] != old_files or plan["adapter_files"] != new_files or
            any(new_files.get(name) != sha for name, sha in old_files.items()) or
            plan["policy_commit"] != POLICY_HEAD or plan["policy_sha256"] != POLICY_SHA or
            plan["numerical_source_fingerprint"] != FROZEN_FP or
            file_hash(frozen / "config/historical-sia-integral-contract.json") != PARENT_SHA):
        raise ValueError("source-plan code or policy identity changed")
    request_path = path.parent / "request.json"
    if str(request_path) != plan["request"] or file_hash(request_path) != plan["request_sha256"]:
        raise ValueError("changed source-plan request bytes or path")
    original = Path(plan["source_request"])
    if file_hash(original) != plan["source_request_sha256"]:
        raise ValueError("changed original source request")
    request = read_json(request_path)
    geometries, selected = geometry_partition(request, policy)
    expected = {"geometries": geometries, "derivative_geometry_indices": selected,
                "point_order": [point["id"] for point in request["points"]],
                "source_inputs": request["source_inputs"], "proposed_phases": policy["job_phases"],
                "resources": policy["resources"], "limitations": request["limitations"]}
    if any(canonical_hash(plan[key]) != canonical_hash(value) for key, value in expected.items()):
        raise ValueError("changed source-derived geometry plan or resource policy")
    api, provenance = _trusted_api(frozen)
    if provenance.source_fingerprint(frozen) != FROZEN_FP:
        raise ValueError("changed frozen numerical fingerprint")
    replay = path.parent / "source-replays" / uuid.uuid4().hex
    replay.mkdir(parents=True, exist_ok=False)
    receipt = {"status": "FAIL", "scope": "source_plan_metadata_verification_not_numerical_admission",
               "externally_bound_plan_sha256": expected_sha256, **dict.fromkeys(FLAGS, False),
               "accepted_fitting_gate": False, "numerical_dispatch_ready": False}
    try:
        source = request["source_preparation"]
        fresh = replay / "request.json"
        api.prepare_request(source["manifest"], source["cache"], source["candidates"], fresh)
        if (canonical_hash(read_json(fresh)) != canonical_hash(request) or
                canonical_hash(read_json(original)) != canonical_hash(request) or
                file_hash(path) != expected_sha256 or
                any(file_hash(p) != sha for p, sha in request["source_inputs"].items()) or
                verify_checkout(frozen, head=old_head)[1] != old_files or
                verify_checkout(adapter, head=new_head)[1] != new_files):
            raise ValueError("source replay or current complete source identities differ")
        receipt.update(status="SOURCE_PLAN_REPLAY_VERIFIED_NOT_NUMERICAL_DISPATCH_OR_GATE",
                       fresh_request_sha256=file_hash(fresh), geometries=36, points=72)
        return receipt
    except BaseException as error:
        receipt["error"] = type(error).__name__ + ": " + str(error)
        raise
    finally:
        receipt["completed_at"] = datetime.now(timezone.utc).isoformat()
        provenance.atomic_json(replay / "replay.json", receipt)
