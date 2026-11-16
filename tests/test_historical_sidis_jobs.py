"""Synthetic metadata/controller tests; no physical loader fallback or receipt."""
import csv
import io
import json
import pytest

from collinsflow import historical_sidis_jobs as jobs


def test_frozen_full132_policy_and_source_indices():
    policy = jobs.load_policy()
    assert policy["execution"]["source_points"] == 132
    assert jobs.point_key(1) == "p0001" and jobs.point_key(132) == "p0132"
    for value in (0, 133, True, 1.0, "1"):
        with pytest.raises(ValueError):
            jobs.point_key(value)


def test_canonical_scope_and_full_product_inventory():
    with pytest.raises(ValueError):
        jobs.folder_for(jobs.ROOT/"data", "native123", 1)
    with pytest.raises(ValueError):
        jobs.folder_for(jobs.ROOT/"runs/synthetic", "diagnostic", False)
    folder = jobs.folder_for(jobs.ROOT/"runs/synthetic", "diagnostic")
    assert len(jobs._products(folder, False)) == 6
    assert len(jobs._products(folder, True)) == 10


@pytest.mark.parametrize("selected", [0, -1, True, float("nan"), float("inf"), 101])
def test_caller_cannot_loosen_or_invent_engineering_budget(selected):
    with pytest.raises(ValueError):
        jobs._deadline(100, selected)


def synthetic_metadata(tmp_path):
    request = {"parameters": {"fixture_scope": "explicitly_synthetic_metadata_not_physical_parameters"},
               "points": [{"id": f"SYNTHETIC_METADATA:{i}"} for i in range(132)]}
    flags = dict.fromkeys(jobs.FLAGS, False)
    flags.update(accepted_prediction_diagnostic=True, all132_coverage_verified=True)
    report = {"scope": jobs.SCIENTIFIC_SCOPE, "phase": "diagnostic", "source_fingerprint": "synthetic_fp",
              "fixture_scope": "explicitly_synthetic_metadata_no_pipeline_receipt", "source_points": 132,
              "parameters": request["parameters"], "request_sha256": "synthetic_request_sha", "policy_sha256": jobs.POLICY_SHA,
              "status": "CONDITIONAL_FULL132_POINT_DIAGNOSTIC_RECORDED_NOT_FITTING", "coordinates": 13,
              "full132_native_tail_reference_derivative_policy_completed": True, "native_level": 3, "reference_level": 3,
              "source_ids": [p["id"] for p in request["points"]],
              "points": [{**p, "native_level": 3, "reference_level": 3} for p in request["points"]], **flags}
    (tmp_path/"report.json").write_text(json.dumps(report), encoding="utf-8")
    audit = {"scope": jobs.SCIENTIFIC_SCOPE, "status": jobs.AUDITED, "source_fingerprint": "synthetic_fp",
             "request_sha256": "synthetic_request_sha", "policy_sha256": jobs.POLICY_SHA,
             "report_sha256": jobs.file_hash(tmp_path/"report.json"), **flags}
    (tmp_path/"audit.json").write_text(json.dumps(audit), encoding="utf-8")
    return request, report, audit


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "mixed_native", "mixed_reference", "missing_coordinate", "production"])
def test_synthetic_metadata_cannot_admit_incomplete_or_mixed_full132(tmp_path, mutation):
    request, report, audit = synthetic_metadata(tmp_path)
    if mutation == "missing":
        report["points"].pop()
    elif mutation == "duplicate":
        report["points"][-1]["id"] = report["points"][0]["id"]
    elif mutation == "mixed_native":
        report["points"][-1]["native_level"] = 4
    elif mutation == "mixed_reference":
        report["points"][-1]["reference_level"] = 4
    elif mutation == "missing_coordinate":
        report["coordinates"] = 12
    else:
        report["production_ready"] = True
    (tmp_path/"report.json").write_text(json.dumps(report), encoding="utf-8")
    audit["report_sha256"] = jobs.file_hash(tmp_path/"report.json")
    (tmp_path/"audit.json").write_text(json.dumps(audit), encoding="utf-8")
    with pytest.raises(ValueError):
        jobs._report_audit(tmp_path, "diagnostic", 0, request, "synthetic_fp", "synthetic_request_sha")


def test_failed_predecessor_rejects_forged_marker_before_dispatch(tmp_path, monkeypatch):
    policy = jobs.load_policy()
    request = {"source_preparation": {"manifest": str(tmp_path/"config/historical-sources.json")},
               "source_fingerprint": "explicitly_synthetic_test_layout"}
    monkeypatch.setattr(jobs, "ROOT", tmp_path)
    monkeypatch.setattr(jobs, "load_policy", lambda: policy)
    monkeypatch.setattr(jobs, "source_fingerprint", lambda root: request["source_fingerprint"])
    monkeypatch.setattr(jobs, "verify_request", lambda path: request)
    def rejected(*args, **kwargs):
        raise ValueError("explicitly synthetic failed predecessor")
    monkeypatch.setattr(jobs, "verify_metadata", rejected)
    def forbidden_dispatch(*args, **kwargs):
        raise AssertionError("numerical dispatch occurred after a forged predecessor marker")
    monkeypatch.setattr(jobs, "_timed_job", forbidden_dispatch)
    with pytest.raises(ValueError, match="failed predecessor"):
        jobs.run_complete_job(tmp_path/"request.json", tmp_path/"native.json", tmp_path/"parton.json",
                              tmp_path/"runs/control", "native123", 2)
    folder = tmp_path/"runs/control/native123/p0002"
    invocation = jobs.read_json(folder/"invocation.json")
    assert invocation["status"] == "FAIL" and invocation["steps"] == []
    assert "completed_at" in invocation
    assert jobs.read_json(folder/"initial-invocation.json")["record_role"] == "initial_sentinel_not_an_executed_failure"


def test_csv_preserves_unknown_sigma_composition_without_filling_statistical_error():
    point = {"id": "SYNTHETIC_CSV:1", "dataset_id": "SYNTHETIC", "publication": "1408.4405", "target": "proton",
             "hadron": "pi+", "projection": "x", "kinematics": {"x": .1, "z": .3, "Q2_GeV2": 2., "pT_GeV": .4, "y": .5},
             "measured": {"central_value": .1, "statistical_error": None, "systematic_components": {}},
             "source_row": {"reported_error_u": .02, "uncertainty_status": "reported_error_u_composition_unresolved"},
             "model_convention": {"formula": "-R"}, "UU": 1., "FUT": -.03, "amplitude": .03,
             "unweighted_observed_minus_predicted": .07, "native_level": 3, "reference_level": 3}
    row = next(csv.DictReader(io.StringIO(jobs.comparison_csv({"points": [point]}))))
    assert row["statistical_error"] == "" and row["reported_error_u"] == "0.02"
    assert row["uncertainty_status"] == "reported_error_u_composition_unresolved"
    assert row["systematic_components"] == "{}" and row["source_formula"] == "-R"
