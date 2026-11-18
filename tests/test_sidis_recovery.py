"""Synthetic recovery corruption and scheduling controls; no physical admission."""
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from collinsflow import sidis_recovery as r

L = r.lineage


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def fake_pair(root, index=1, phase="native123"):
    folder = root / phase / f"p{index:04d}"
    inv = dict(status=L.COMPLETED, phase=phase, point_index=index, completed_at="SYNTHETIC_NOT_A_TIME",
               steps=[dict(name=name, returncode=0, process_returncode=0, interruption=None, termination_signal=None)
                      for name in ("producer", "audit")], scope="explicitly_synthetic_incomplete_metadata_not_admission")
    write(folder / "invocation.json", inv)
    (folder / "job.completed").write_text("SYNTHETIC_ONLY")
    return folder, inv


@pytest.mark.parametrize("mutation", ("bytes", "size", "added", "missing", "duplicate", "relative-escape", "alias", "link-target"))
def test_immutable_inventory_includes_links_without_following_them(tmp_path, mutation):
    root = tmp_path / "old"
    root.mkdir()
    path = root / "record.json"
    path.write_text("SYNTHETIC_IMMUTABLE")
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "must-not-be-read").write_text("OUTSIDE")
    link = root / "retained-link"
    link.symlink_to(outside, target_is_directory=True)
    rows = L.inventory(root)
    assert len(rows) == 2
    assert L.verify_inventory(root, rows, tmp_path) == rows
    if mutation == "bytes":
        path.write_text("CHANGED")
    elif mutation == "size":
        rows[0]["bytes"] += 1
    elif mutation == "added":
        (root / "added").write_text("NEW")
    elif mutation == "missing":
        path.unlink()
    elif mutation == "duplicate":
        rows.append(copy.deepcopy(rows[0]))
    elif mutation == "relative-escape":
        rows[0]["path"] = "old/../outside/must-not-be-read"
    elif mutation == "alias":
        path.unlink()
        path.symlink_to(outside / "must-not-be-read")
    else:
        link.unlink()
        link.symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(ValueError):
        L.verify_inventory(root, rows, tmp_path)
    assert (outside / "must-not-be-read").read_text() == "OUTSIDE"


@pytest.mark.parametrize("mutation", ("raw-failure", "effective-failure", "timeout", "signal", "boolean-exit", "missing-auditor", "marker", "wrong-index"))
def test_status_or_marker_cannot_approve_an_incomplete_or_interrupted_pair(tmp_path, mutation):
    folder, inv = fake_pair(tmp_path)
    if mutation == "raw-failure":
        inv["steps"][0]["process_returncode"] = 2
    elif mutation == "effective-failure":
        inv["steps"][0]["returncode"] = 124
    elif mutation == "timeout":
        inv["steps"][0]["interruption"] = "deadline"
    elif mutation == "signal":
        inv["steps"][0]["termination_signal"] = 15
    elif mutation == "boolean-exit":
        inv["steps"][0]["returncode"] = False
    elif mutation == "missing-auditor":
        inv["steps"].pop()
    elif mutation == "marker":
        (folder / "job.completed").unlink()
    else:
        inv["point_index"] = True
    write(folder / "invocation.json", inv)
    with pytest.raises(ValueError):
        L.complete_locations(tmp_path, "recovery_v2", tmp_path / "recovery.json", "0" * 64)


def test_stopped_partial_is_preserved_but_never_selected_and_terminal_failure_blocks(tmp_path):
    folder, inv = fake_pair(tmp_path)
    inv.update(status="FAIL", completed_at=None)
    (folder / "job.completed").unlink()
    write(folder / "invocation.json", inv)
    progress = folder / "report.json.progress.json"
    progress.write_text("SYNTHETIC_PARTIAL_LEVEL1_ONLY")
    before = L.inventory(tmp_path)
    assert L.complete_locations(tmp_path, "recovery_v2", tmp_path / "recovery.json", "0" * 64) == []
    assert L.inventory(tmp_path) == before
    inv["completed_at"] = "SYNTHETIC_TERMINAL_FAILURE"
    write(folder / "invocation.json", inv)
    with pytest.raises(ValueError, match="terminal failed point"):
        L.complete_locations(tmp_path, "recovery_v2", tmp_path / "recovery.json", "0" * 64)


def test_failed_global_gate_is_not_recovered_as_a_fresh_chance(tmp_path):
    write(tmp_path / "native-select3/invocation.json", dict(status="FAIL", completed_at="SYNTHETIC_TERMINAL"))
    with pytest.raises(ValueError, match="global gate"):
        L.complete_locations(tmp_path, "recovery_v2", tmp_path / "recovery.json", "0" * 64)


def test_duplicate_or_nonprefix_jobs_cannot_define_a_retained_source_set():
    one = dict(phase="native123", index=1)
    assert L.ordered([one]) == [one]
    for rows in ([one, one], [dict(phase="native123", index=2)]):
        with pytest.raises(ValueError):
            L.ordered(rows)


@pytest.fixture
def synthetic_ancestry(tmp_path, monkeypatch):
    frozen = tmp_path / "frozen"
    predecessor = tmp_path / "adapter"
    original = frozen / "runs/historical-sidis-jobs-real-data/attempt-2"
    stopped = frozen / "runs/stopped"
    current = frozen / "runs/current"
    for folder in (predecessor, original / "full", stopped, current):
        folder.mkdir(parents=True)
    write(original / "invocation.json", dict(scope="explicitly_synthetic_location_fixture_not_a_job"))
    policy = L.paths.policy()
    policy["logical_attempt"]["original_outer_invocation_sha256"] = r.file_hash(original / "invocation.json")
    monkeypatch.setattr(L.paths, "policy", lambda: policy)
    for index in range(1, 14):
        (original / "full/native123" / f"p{index:04d}").mkdir(parents=True)
    parent = {key: "SYNTHETIC_" + key for key in L.ANCHOR_FIELDS}
    parent.update(original=str(original), frozen=str(frozen), adapter=str(predecessor), adapter_commit=L.PREDECESSOR,
        adapter_files={"synthetic.txt": "SYNTHETIC_SHA"}, contract_sha256=L.old.CONTRACT_SHA, scope=L.old.SCOPE,
        status="PREPARED_IMMUTABLE_CONTINUATION_NOT_A_GATE", run=str(stopped), **dict.fromkeys(L.paths.FLAGS, False))
    manifest = stopped / "continuation.json"
    write(manifest, parent)
    for index in range(14, 38):
        fake_pair(stopped, index)
    for index in range(355):
        (stopped / f"synthetic-extra-{index:04d}").write_text("EXPLICITLY_SYNTHETIC_NOT_NUMERICAL")
    assert len(L.inventory(stopped)) == 404
    loss = tmp_path / "loss.json"
    # The fixture substitutes synthetic pins only inside this test process.
    workspace = frozen.parents[2]
    write(loss, dict(old_owner_absent_verified=True, source_run=str(stopped.relative_to(workspace)),
                     recorded_complete_prefix=37, retained_files=L.inventory(stopped)))
    monkeypatch.setattr(L, "FIRST_MANIFEST", r.file_hash(manifest))
    monkeypatch.setattr(L, "FIRST_LOSS", r.file_hash(loss))
    node = dict(protocol="continuation_v1", manifest=str(manifest), manifest_sha256=r.file_hash(manifest),
                inventory=str(loss), inventory_sha256=r.file_hash(loss))
    return frozen, predecessor, original, stopped, current, parent, node


def test_read_only_ancestry_retains37_actual_locations_and_protects_partial38(synthetic_ancestry):
    frozen, predecessor, original, stopped, current, parent, node = synthetic_ancestry
    before = L.inventory(stopped)
    state = L.validate_ancestors([node], frozen, predecessor, current)
    assert len(state["entries"]) == 37
    for index in range(1, 38):
        base = original / "full" if index <= 13 else stopped
        assert state["view"].folder("native123", index) == base / "native123" / f"p{index:04d}"
    assert state["view"].fresh_folder("native123", 38) == current / "native123/p0038"
    assert state["view"].metadata()["complete_job_provenance_verified"] is False
    assert L.inventory(stopped) == before


@pytest.mark.parametrize("mutation", ("repeat", "wrong-pin", "wrong-family", "stale-inventory", "depth"))
def test_untrusted_ancestry_cannot_relabel_sources_or_reuse_a_namespace(synthetic_ancestry, mutation):
    frozen, predecessor, original, stopped, current, parent, node = synthetic_ancestry
    nodes = [node]
    if mutation == "repeat":
        nodes.append(copy.deepcopy(node))
    elif mutation == "wrong-pin":
        node["manifest_sha256"] = "0" * 64
    elif mutation == "wrong-family":
        parent["frozen"] = "SYNTHETIC_OTHER_FAMILY"
        # Pin remains unchanged, so byte tampering fails before family selection.
        write(Path(node["manifest"]), parent)
    elif mutation == "stale-inventory":
        (stopped / "unrecorded").write_text("EXTRA")
    else:
        nodes *= 65
    with pytest.raises(ValueError):
        L.validate_ancestors(nodes, frozen, predecessor, current)


@pytest.mark.parametrize("mutation", (None, "different-family", "cycle", "conflicting-point"))
def test_a_later_shutdown_can_retain_whole_new_jobs_but_not_partial_levels(synthetic_ancestry, monkeypatch, mutation):
    frozen, predecessor, original, stopped, current, first, node = synthetic_ancestry
    state = L.validate_ancestors([node], frozen, predecessor, current)
    parent = {key: first[key] for key in L.ANCHOR_FIELDS}
    parent.update(status=L.PREPARED, scope=L.SCOPE, contract_sha256=L.paths.POLICY_SHA, run=str(current),
        adapter=str(predecessor), adapter_commit="EXPLICITLY_SYNTHETIC_NOT_A_COMMIT", adapter_files=first["adapter_files"],
        ancestors=[node], retained_jobs=state["entries"], **dict.fromkeys(L.paths.FLAGS, False))
    fake_pair(current, 38)
    write(current / "native123/p0039/report.json.progress.json", dict(scope="SYNTHETIC_PARTIAL39_LEVEL1_NOT_A_PAIR"))
    if mutation == "different-family":
        parent["native"] = "EXPLICITLY_SYNTHETIC_DIFFERENT_NATIVE_BUILD"
    elif mutation == "cycle":
        parent["ancestors"] = [node, node]
    elif mutation == "conflicting-point":
        fake_pair(current, 14)
    manifest = current / "recovery.json"
    write(manifest, parent)
    loss = current.parent / "second-loss.json"
    write(loss, dict(status="STOPPED_SIDIS_NAMESPACE_INVENTORIED_NOT_A_GATE", manifest_sha256=r.file_hash(manifest),
        run=str(current), old_owner_absent_verified=True, retained_files=L.inventory(current)))
    next_node = dict(protocol="recovery_v2", manifest=str(manifest), manifest_sha256=r.file_hash(manifest),
                     inventory=str(loss), inventory_sha256=r.file_hash(loss))
    monkeypatch.setattr(L.old, "verify_git", lambda *args: (parent["adapter_commit"], parent["adapter_files"]))
    latest = current.parent / "latest"
    latest.mkdir()
    before = L.inventory(current)
    if mutation is not None:
        with pytest.raises(ValueError):
            L.validate_ancestors([node, next_node], frozen, predecessor, latest)
    else:
        result = L.validate_ancestors([node, next_node], frozen, predecessor, latest)
        assert len(result["entries"]) == 38
        assert result["view"].folder("native123", 38) == current / "native123/p0038"
        assert result["view"].fresh_folder("native123", 39) == latest / "native123/p0039"
        assert result["view"].folder("native-selection") == latest / "native-selection"
    assert L.inventory(current) == before


@pytest.mark.parametrize("native4,reference4", [(False, False), (True, False), (False, True), (True, True)])
def test_recovery_still_schedules_the_entire132_source_policy(tmp_path, monkeypatch, native4, reference4):
    root = tmp_path / "fresh"
    root.mkdir()
    (root / "retained-review").mkdir()
    m = dict(run=str(root), retained_jobs=[dict(phase="native123", index=i) for i in range(1, 38)])
    calls = []
    class SyntheticView:
        def retained(self, phase, index):
            return dict(scope="synthetic_location_only") if phase == "native123" and index <= 37 else None
        def folder(self, phase, index=0):
            return root / phase / f"p{index:04d}"
    monkeypatch.setattr(r, "view", lambda m: SyntheticView())
    monkeypatch.setattr(r, "point_folder", lambda m, phase, index=0: SyntheticView().folder(phase, index))
    monkeypatch.setattr(r, "verify_retained_review", lambda *args: None)
    monkeypatch.setattr(r, "run_job", lambda m, path, phase, index, api, lease: calls.append((phase, index)))
    monkeypatch.setattr(r, "audit_cached", lambda *args: dict(status="SYNTHETIC_SCHEDULING_TRACE_NOT_NUMERICAL_GATE"))
    monkeypatch.setattr(r, "selection", lambda m, phase, api: dict(selected_level=0 if native4 else 3) if phase == "native-select3"
                        else dict(selected_reference_level=0 if reference4 else 3))
    r.run_all(m, root / "recovery.json", None, dict(scope="synthetic_lease_not_live"))
    assert [i for phase, i in calls if phase == "native123"] == list(range(38, 133))
    assert [i for phase, i in calls if phase == "native4"] == (list(range(1, 133)) if native4 else [])
    assert [i for phase, i in calls if phase == "tailref123"] == list(range(1, 133))
    assert [i for phase, i in calls if phase == "reference4"] == (list(range(1, 133)) if reference4 else [])
    assert [phase for phase, index in calls if index == 0] == list(r.GLOBAL_PHASES)
    assert all(index != 38 or phase != "native123" or calls.count((phase, index)) == 1 for phase, index in calls)


@pytest.mark.parametrize("elapsed", (True, float("nan"), float("inf"), -1))
def test_nonfinite_or_boolean_resource_duration_cannot_approve_a_command(tmp_path, elapsed):
    command = ["/synthetic/julia", "explicitly_synthetic_nonexecuted.jl"]
    row = dict(name="producer", command=r.measured_command(tmp_path, "producer", command), returncode=0,
               process_returncode=0, wrapper_elapsed_seconds=elapsed)
    with pytest.raises(ValueError, match="elapsed time"):
        r.verify_step(row, tmp_path, command, None, SimpleNamespace())


def test_readonly_recheck_rejects_manifest_mutation_before_any_writer_or_reader(tmp_path):
    m = dict(status="EXPLICITLY_SYNTHETIC_NOT_A_REAL_MANIFEST")
    path = tmp_path / "recovery.json"
    write(path, {**m, "fit_performed": True})
    with pytest.raises(ValueError, match="changed during"):
        r.recheck(m, path, None)
