"""Explicitly synthetic source-planning and raw-source identity controls."""
import copy
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

from collinsflow import sia_geometry_jobs as jobs
from collinsflow.provenance import atomic_json, read_json
from test_historical_sia_integrals import synthetic_candidate
from collinsflow.historical_sia_integrals import request_points


def fixture_request():
    candidate, contract, initial = synthetic_candidate()
    return dict(scope=jobs.REQUEST_SCOPE, theory_scope="reconstructed_theory", source_fingerprint=jobs.FROZEN_FP,
                fixture_scope="explicitly_synthetic_source_representation_not_primary_data",
                points=request_points(candidate, contract, initial), **dict.fromkeys(jobs.FLAGS, False))


def policy():
    return jobs.load_policy()


def test_exact_partition_preserves_both_source_angles_and_owned_metadata():
    request = fixture_request()
    request["points"][36]["source_angular"]["value"] = .51
    before = copy.deepcopy(request)
    groups, selected = jobs.geometry_partition(request, policy())
    assert len(groups) == 36 and selected == [1, 15, 36]
    assert [point["request_index"] for point in groups[0]["points"]] == [1, 37]
    assert [point["angular_source"]["value"] for point in groups[0]["points"]] == [.6, .51]
    assert groups[0]["points"][0]["source_point_sha256"] != groups[0]["points"][1]["source_point_sha256"]
    groups[0]["points"][0]["angular_source"]["value"] = .1
    assert request == before


@pytest.mark.parametrize("flag", jobs.FLAGS)
def test_source_or_missing_flags_cannot_self_admit(flag):
    request = fixture_request()
    request[flag] = True
    with pytest.raises(ValueError):
        jobs.geometry_partition(request, policy())
    request.pop(flag)
    with pytest.raises(ValueError):
        jobs.geometry_partition(request, policy())


@pytest.mark.parametrize("mutation", ["omit", "duplicate", "order", "geometry", "rounded-geometry", "ratio",
                                     "scope", "fingerprint", "nan", "bool", "angular", "angular-bool"])
def test_full_source_partition_rejects_substitution(mutation):
    request = fixture_request()
    point = request["points"][71]
    if mutation == "omit":
        request["points"].pop()
    elif mutation == "duplicate":
        request["points"][71] = copy.deepcopy(request["points"][0])
    elif mutation == "order":
        request["points"].reverse()
    elif mutation in ("geometry", "rounded-geometry"):
        point["kinematics"]["z1"] += .01 if mutation == "geometry" else 1e-15
    elif mutation == "ratio":
        point["ratio"] = "UC"
    elif mutation == "scope":
        request["scope"] = "production"
    elif mutation == "fingerprint":
        request["source_fingerprint"] = "0" * 64
    elif mutation == "nan":
        point["kinematics"]["z1"] = float("nan")
    elif mutation == "bool":
        point["kinematics"]["Q2_GeV2"] = True
    elif mutation == "angular":
        point["source_angular"]["average_type"] = "ratio_of_angular_means"
    else:
        point["source_angular"]["value"] = False
    with pytest.raises(ValueError):
        jobs.geometry_partition(request, policy())


def test_canonical_source_comparison_distinguishes_boolean_and_float():
    assert {"lower": 0.0} == {"lower": False}
    assert jobs.canonical_hash({"lower": 0.0}) != jobs.canonical_hash({"lower": False})


def test_raw_blob_verification_rejects_hidden_index_edits(tmp_path):
    def git(*args):
        return subprocess.check_output(["git", "-C", str(tmp_path), *args], stderr=subprocess.PIPE)
    git("init")
    git("config", "user.name", "Synthetic identity fixture")
    git("config", "user.email", "synthetic@example.invalid")
    source = tmp_path / "kernel.txt"
    source.write_text("explicitly synthetic committed source\n", encoding="utf-8")
    git("add", "kernel.txt")
    git("commit", "-m", "synthetic identity fixture")
    head, inventory = jobs.verify_checkout(tmp_path)
    assert inventory == {"kernel.txt": jobs.file_hash(source)}
    git("update-index", "--assume-unchanged", "kernel.txt")
    source.write_text("changed hidden source\n", encoding="utf-8")
    assert not git("diff", "--name-only").strip()
    with pytest.raises(ValueError, match="raw committed blob"):
        jobs.verify_checkout(tmp_path, head=head)


def test_planning_requires_fresh_canonical_new_output(tmp_path):
    occupied = tmp_path / "runs" / "occupied"
    occupied.mkdir(parents=True)
    sentinel = occupied / "old.json"
    sentinel.write_text('{"status":"FAIL"}\n', encoding="utf-8")
    before = sentinel.read_bytes()
    with pytest.raises(ValueError, match="fresh canonical"):
        jobs.prepare_plan(tmp_path, tmp_path / "request.json", occupied, adapter=tmp_path)
    assert sentinel.read_bytes() == before
    with pytest.raises(ValueError, match="fresh canonical"):
        jobs.prepare_plan(tmp_path, tmp_path / "request.json", tmp_path / "outside", adapter=tmp_path)


def test_policy_is_the_published_prospective_artifact():
    data = policy()
    assert data["frozen_numerics"]["head"] == jobs.FROZEN_HEAD
    assert data["frozen_numerics"]["source_fingerprint"] == jobs.FROZEN_FP
    assert read_json(jobs.ROOT / "config/historical-sia-integral-contract.json")["selection"]["total_points"] == 72


@pytest.mark.parametrize("flag", (*jobs.FLAGS, "accepted_fitting_gate", "numerical_dispatch_ready", "physical_jobs_started"))
def test_rehashed_source_plan_cannot_authorize_numerical_jobs(tmp_path, flag):
    path = tmp_path / "plan.json"
    data = dict(schema_version=1, status=jobs.PLAN_STATUS, scope=jobs.SCOPE, theory_scope="reconstructed_theory",
                **dict.fromkeys((*jobs.FLAGS, "accepted_fitting_gate", "numerical_dispatch_ready", "physical_jobs_started"), False))
    data[flag] = True
    path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match="cannot admit"):
        jobs.verify_plan(path, expected_sha256=jobs.file_hash(path))


def test_plan_verification_requires_external_hash_before_any_replay(tmp_path):
    path = tmp_path / "plan.json"
    path.write_text("{}\n", encoding="utf-8")
    with pytest.raises(ValueError, match="externally bound"):
        jobs.verify_plan(path, expected_sha256="0" * 64)
    assert set(tmp_path.iterdir()) == {path}


@pytest.mark.parametrize("version", [None, False, True, 1.0, 2])
def test_plan_schema_rejects_rehashed_boolean_or_other_versions(tmp_path, version):
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(dict(schema_version=version, status=jobs.PLAN_STATUS, scope=jobs.SCOPE,
        theory_scope="reconstructed_theory", **dict.fromkeys((*jobs.FLAGS, "accepted_fitting_gate",
        "numerical_dispatch_ready", "physical_jobs_started"), False))), encoding="utf-8")
    with pytest.raises(ValueError, match="cannot admit"):
        jobs.verify_plan(path, expected_sha256=jobs.file_hash(path))


@pytest.fixture
def synthetic_planning_context(tmp_path, monkeypatch):
    """Mock only source/Git boundaries; write no real source or numerical file."""
    policy_data = jobs.load_policy()
    parent = (jobs.ROOT / "config/historical-sia-integral-contract.json").read_bytes()
    frozen, adapter = tmp_path / "frozen", tmp_path / "adapter"
    (frozen / "config").mkdir(parents=True)
    (frozen / "config/historical-sia-integral-contract.json").write_bytes(parent)
    files = {f"explicitly-synthetic-source-{i}": str(i) for i in range(302)}

    def checkout(path, **kwargs):
        return (jobs.FROZEN_HEAD if Path(path) == frozen else "synthetic-additive-head", copy.deepcopy(files))

    monkeypatch.setattr(jobs, "ROOT", adapter)
    monkeypatch.setattr(jobs, "load_policy", lambda *args: policy_data)
    monkeypatch.setattr(jobs, "verify_checkout", checkout)
    request = fixture_request()
    request.update(source_preparation=dict(manifest="synthetic", cache="synthetic", candidates="synthetic"),
                   source_inputs={}, limitations=["explicitly synthetic metadata race fixture"])
    original = tmp_path / "source.json"
    atomic_json(original, request)

    def prepare(manifest, cache, candidates, output):
        atomic_json(output, request)

    api = SimpleNamespace(prepare_request=prepare)
    provenance = SimpleNamespace(atomic_json=atomic_json, source_fingerprint=lambda path: jobs.FROZEN_FP)
    monkeypatch.setattr(jobs, "_trusted_api", lambda path: (api, provenance))
    return dict(original=original, folder=adapter / "runs/attempt-1", adapter=adapter, frozen=frozen,
                api=api, checkout=checkout, prepare=prepare)


@pytest.mark.parametrize("target", ["original", "saved"])
def test_request_bytes_changed_during_final_replay_are_rejected(synthetic_planning_context, target):
    context = synthetic_planning_context
    original, folder = context["original"], context["folder"]
    jobs.prepare_plan(context["frozen"], original, folder, adapter=context["adapter"])
    plan = folder / "plan.json"
    bound_sha = jobs.file_hash(plan)

    def mutate(manifest, cache, candidates, output):
        context["prepare"](manifest, cache, candidates, output)
        changed = original if target == "original" else folder / "request.json"
        changed.write_text(changed.read_text(encoding="utf-8") + "\n", encoding="utf-8")

    context["api"].prepare_request = mutate
    with pytest.raises(ValueError, match="source replay|request bytes"):
        jobs.verify_plan(plan, expected_sha256=bound_sha)
    receipts = list((folder / "source-replays").glob("*/replay.json"))
    assert len(receipts) == 1 and read_json(receipts[0])["status"] == "FAIL"
    assert read_json(plan)["numerical_dispatch_ready"] is False


@pytest.mark.parametrize("target", ["original", "saved"])
def test_request_bytes_changed_during_final_source_scan_are_rejected(synthetic_planning_context, monkeypatch, target):
    context = synthetic_planning_context
    calls = 0

    def mutate(path, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 4:
            changed = context["original"] if target == "original" else context["folder"] / "request.json"
            changed.write_text(changed.read_text(encoding="utf-8") + "\n", encoding="utf-8")
        return context["checkout"](path, **kwargs)

    monkeypatch.setattr(jobs, "verify_checkout", mutate)
    with pytest.raises(ValueError, match="source identities changed"):
        jobs.prepare_plan(context["frozen"], context["original"], context["folder"], adapter=context["adapter"])
    assert read_json(context["folder"] / "invocation.json")["status"] == "FAIL"
    assert not (context["folder"] / "plan.json").exists()
