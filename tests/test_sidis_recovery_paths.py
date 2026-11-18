import copy
import json

import pytest

from collinsflow import sidis_recovery_paths as paths


def synthetic_case(tmp_path):
    original = tmp_path / "original" / "full"
    predecessor = tmp_path / "predecessor"
    current = tmp_path / "current"
    for directory in (original, predecessor, current):
        directory.mkdir(parents=True)
    original_manifest = original.parent / "invocation.json"
    predecessor_manifest = predecessor / "continuation.json"
    for path in (original_manifest, predecessor_manifest):
        path.write_text(json.dumps({"scope": "explicitly_synthetic_location_fixture_not_a_job", "name": path.name}))
    rows = []
    for index in range(1, 38):
        producer, manifest, protocol = ((original, original_manifest, "frozen_original") if index <= 13 else
                                          (predecessor, predecessor_manifest, "continuation_v1"))
        folder = producer / "native123" / f"p{index:04d}"
        folder.mkdir(parents=True)
        rows.append(dict(phase="native123", index=index, folder=str(folder), producer_root=str(producer),
                         producer_manifest=str(manifest), producer_manifest_sha256=paths.digest(manifest), protocol=protocol))
    partial = predecessor / "native123/p0038/report.json"
    partial.parent.mkdir()
    partial.write_text('{"scope":"explicitly_synthetic_partial_fixture","status":"FAIL"}\n')
    return current, original, predecessor, rows, partial


def test_exact37_locations_and_missing38_never_rewrite_or_admit(tmp_path):
    current, original, predecessor, rows, partial = synthetic_case(tmp_path)
    before = {p: paths.digest(p) for p in tmp_path.rglob("*") if p.is_file()}
    run = paths.RecoveryRunPaths(current, rows)
    for phase in paths.POINT_PHASES:
        for index in range(1, 133):
            expected = original if phase == "native123" and index <= 13 else predecessor if phase == "native123" and index <= 37 else current
            assert run.folder(phase, index) == expected / phase / f"p{index:04d}"
    for phase in paths.GLOBAL_PHASES:
        assert run.folder(phase) == current / phase
    assert run.fresh_folder("native123", 38) == current / "native123/p0038"
    assert run.retained("native123", 38) is None and partial.is_file()
    metadata = run.metadata()
    assert metadata["scope"] == "owned_whole_job_path_view_not_admission"
    assert metadata["numerical_dispatch_ready"] is False
    assert metadata["complete_job_provenance_verified"] is False
    assert metadata["saved_numerical_arithmetic_verified"] is False
    assert all(metadata[key] is False for key in paths.FLAGS)
    assert all(paths.digest(path) == sha for path, sha in before.items())


@pytest.mark.parametrize("mutation", ["boolean", "float", "phase", "missing-field", "admission-field", "protocol", "point-folder", "producer-root", "sha", "sha-type", "order", "duplicate"])
def test_wrong_identity_schema_or_map_rejected(tmp_path, mutation):
    current, _, _, rows, _ = synthetic_case(tmp_path)
    bad = copy.deepcopy(rows)
    if mutation == "boolean":
        bad[0]["index"] = True
    elif mutation == "float":
        bad[0]["index"] = 1.0
    elif mutation == "phase":
        bad[0]["phase"] = "../native123"
    elif mutation == "missing-field":
        del bad[0]["protocol"]
    elif mutation == "admission-field":
        bad[0]["accepted_fitting_gate"] = True
    elif mutation == "protocol":
        bad[0]["protocol"] = "unreviewed"
    elif mutation == "point-folder":
        bad[-1]["folder"] = bad[-2]["folder"]
    elif mutation == "producer-root":
        bad[-1]["producer_root"] = str(current)
    elif mutation == "sha":
        bad[-1]["producer_manifest_sha256"] = "0" * 64
    elif mutation == "sha-type":
        bad[-1]["producer_manifest_sha256"] = True
    elif mutation == "order":
        bad.reverse()
    else:
        bad.append(copy.deepcopy(bad[-1]))
    with pytest.raises(ValueError):
        paths.RecoveryRunPaths(current, bad)


def test_owned_snapshot_and_producing_manifest_rechecked(tmp_path):
    current, _, predecessor, rows, _ = synthetic_case(tmp_path)
    run = paths.RecoveryRunPaths(current, rows)
    rows[0]["index"] = 999
    assert run.folder("native123", 1).name == "p0001"
    returned = run.retained("native123", 1)
    returned["folder"] = "changed"
    assert run.folder("native123", 1).name == "p0001"
    manifest = predecessor / "continuation.json"
    original_bytes = manifest.read_bytes()
    manifest.write_bytes(original_bytes + b"\n")
    with pytest.raises(ValueError):
        run.folder("native123", 38)
    manifest.write_bytes(original_bytes)
    run.entries[-1]["index"] = 36
    with pytest.raises(ValueError):
        run.valid()


def test_aliases_overlap_escape_and_occupied_failure_rejected(tmp_path):
    current, original, predecessor, rows, _ = synthetic_case(tmp_path)
    run = paths.RecoveryRunPaths(current, rows)
    with pytest.raises(ValueError):
        run.fresh_folder("native123", 14)
    nested = original.parent / "other-output"
    nested.mkdir()
    with pytest.raises(ValueError):
        paths.RecoveryRunPaths(nested, rows)
    alias = tmp_path / "alias"
    alias.symlink_to(predecessor, target_is_directory=True)
    bad = copy.deepcopy(rows)
    bad[-1]["folder"] = str(alias / "native123/p0037")
    with pytest.raises(ValueError):
        paths.RecoveryRunPaths(current, bad)
    with pytest.raises(ValueError):
        run.current_path("..", "escaped")
    (current / "alias").symlink_to(predecessor, target_is_directory=True)
    with pytest.raises(ValueError):
        run.current_path("alias", "new-output")
    failure = current / "native123/p0038/report.json"
    failure.parent.mkdir(parents=True)
    failure.write_text("explicitly synthetic occupied failure")
    sha = paths.digest(failure)
    with pytest.raises(ValueError):
        run.fresh_folder("native123", 38)
    assert paths.digest(failure) == sha


def test_later_phase_locations_keep_their_exact_producing_protocol(tmp_path):
    current, _, _, rows, _ = synthetic_case(tmp_path)
    origin = tmp_path / "older-recovery"
    origin.mkdir()
    manifest = origin / "recovery.json"
    manifest.write_text('{"scope":"explicitly_synthetic_recovery_location_fixture"}\n')
    for phase in paths.POINT_PHASES[1:]:
        folder = origin / phase / "p0002"
        folder.mkdir(parents=True)
        rows.append(dict(phase=phase, index=2, folder=str(folder), producer_root=str(origin), producer_manifest=str(manifest),
                         producer_manifest_sha256=paths.digest(manifest), protocol="recovery_v2"))
    run = paths.RecoveryRunPaths(current, rows)
    for phase in paths.POINT_PHASES[1:]:
        assert run.folder(phase, 2) == origin / phase / "p0002"
        with pytest.raises(ValueError):
            run.fresh_folder(phase, 2)
    assert run.metadata()["saved_numerical_arithmetic_verified"] is False


@pytest.mark.parametrize("phase,index", [("native123", False), ("native123", 0), ("native123", 133), ("diagnostic", True), ("diagnostic", 1), ("unknown", 0)])
def test_invalid_requested_phase_index_is_not_a_path(tmp_path, phase, index):
    current, _, _, rows, _ = synthetic_case(tmp_path)
    run = paths.RecoveryRunPaths(current, rows)
    with pytest.raises(ValueError):
        run.folder(phase, index)
