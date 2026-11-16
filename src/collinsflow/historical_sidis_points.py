"""Source preparation for the full132 fixed-mean SIDIS diagnostic, without fitting."""
from __future__ import annotations

import copy
from pathlib import Path
import tempfile

from .historical import SCOPE as CANDIDATE_SCOPE, verify_historical_intake
from .historical_points import finite_number
from .provenance import atomic_json, file_hash, read_json, source_fingerprint

SCOPE = "historical_reconstructed_sidis_point_prediction_diagnostic"
CONTRACT = "config/historical-sidis-point-contract.json"
FALSE_FLAGS = ("production_ready", "physical_extraction", "historical_fit_ready",
               "experimental_data_admitted", "fit_performed", "G2PERF_pass", "interval_coverage_verified")
MEASURED_FIELDS = ("central_value", "statistical_error", "systematic_components")


def geometry_key(point):
    """Exact target/means identity; no rounding or charge/projection deduplication."""
    return (point["target"], *(point["kinematics"][k] for k in ("x", "z", "Q2_GeV2", "pT_GeV")))


def request_points(candidate, contract, initial):
    """Validate representation only; prepare_request separately replays raw sources."""
    if (contract.get("id") != "KPSY-HISTORICAL-SIDIS-POINT-001" or contract.get("scope") != SCOPE or
            any(contract.get(k) is not False for k in FALSE_FLAGS) or
            candidate.get("scope") != CANDIDATE_SCOPE or
            any(candidate.get(k) is not False for k in FALSE_FLAGS[:3])):
        raise ValueError("SIDIS diagnostic cannot approve another scope or fitting")
    if (initial["parameter_order"] != contract["prediction"]["parameter_order"] or
            initial["published_parameters"] != contract["prediction"]["parameter_values"] or
            len(initial["parameter_order"]) != 13 or len(set(initial["parameter_order"])) != 13):
        raise ValueError("unchanged thirteen-coordinate parameter source required")
    for value in initial["published_parameters"]:
        finite_number(value, "published parameter")
    views = {v["id"]: v for v in candidate["views"]}
    if len(views) != len(candidate["views"]):
        raise ValueError("duplicate candidate view")
    points, ids, excluded = [], set(), 0
    for selected in contract["selection"]["source_views"]:
        view = views.get(selected["id"])
        convention = contract["model_conventions"][selected["publication"]]
        if (view is None or view["source_id"] != selected["id"] or view["process"] != "SIDIS" or
                any(view[k] != selected[k] for k in ("publication", "target", "hadron", "projection", "source_sha256")) or
                view["depolarization_removed"] is not selected["original_depolarization_removed"] or
                convention["id"] != selected["source_convention_id"] or
                convention["original_depolarization_removed"] is not selected["original_depolarization_removed"] or
                view["source_observable"] != "AUTcollins" or view["output_units"] != "fraction" or
                view["production_ready"] is not False or view["historical_fit_ready"] is not False or
                not view["blockers"] or len(view["points"]) != selected["points"]):
            raise ValueError("mismatched SIDIS source selection or unresolved source label")
        excluded += len(view["excluded"])
        for row in view["points"]:
            row_number = row["locator"]["excel_row"]
            if (isinstance(row_number, bool) or not isinstance(row_number, int) or row_number < 2 or
                    row["id"] != f"{selected['id']}:r{row_number}" or row["id"] in ids):
                raise ValueError("duplicate or changed source row identity")
            ids.add(row["id"])
            means = row["averages"]
            if set(means) != {"x", "y", "z", "Q2_GeV2", "pT_GeV"}:
                raise ValueError("complete unchanged source means required")
            for key in ("x", "y", "z"):
                if not 0 < finite_number(means[key], key) < 1:
                    raise ValueError("invalid source fraction")
            finite_number(means["Q2_GeV2"], "Q2", positive=True)
            finite_number(means["pT_GeV"], "P", nonnegative=True)
            if view["projection"] == "pt" and means["pT_GeV"] >= .8:
                raise ValueError("changed strict reported-P selection")
            finite_number(row["central_value"], "source central value")
            if selected["source_convention_id"] == "COMPASS2014":
                if (row["statistical_error"] is not None or row["systematic_components"] or
                        row.get("uncertainty_status") != "reported_error_u_composition_unresolved"):
                    raise ValueError("unresolved proton Sigma must not become a statistical error")
                finite_number(row["reported_error_u"], "unresolved reported error_u", positive=True)
            else:
                finite_number(row["statistical_error"], "statistical error", positive=True)
                expected = {"experimental"} if selected["source_convention_id"] == "HERMES2010" else set()
                if set(row["systematic_components"]) != expected:
                    raise ValueError("source systematic components changed")
                for value in row["systematic_components"].values():
                    finite_number(value, "source systematic", nonnegative=True)
            points.append({"id": row["id"], "dataset_id": selected["id"],
                **{k: view[k] for k in ("publication", "target", "hadron", "projection", "source_sha256")},
                "charge": {"pi+": 1, "pi-": -1}[view["hadron"]],
                "kinematics": copy.deepcopy(means), "source_row": copy.deepcopy(row),
                "measured": {k: copy.deepcopy(row[k]) for k in MEASURED_FIELDS},
                "original_depolarization_removed": view["depolarization_removed"],
                "model_convention": copy.deepcopy(convention), "source_blockers": list(view["blockers"])})
    selection = contract["selection"]
    if (len(points) != selection["points"] or selection["points"] != 132 or
            len(points)+excluded != selection["raw_points"] or selection["raw_points"] != 146 or
            excluded != selection["excluded_points"] or excluded != 14 or
            len(selection["source_views"]) != selection["views"] or selection["views"] != 18 or
            len(set(map(geometry_key, points))) != contract["prediction"]["distinct_source_tuples"]):
        raise ValueError("incomplete SIDIS selection, exclusion ledger or exact geometry inventory")
    return points


def prepare_request(manifest_path, cache, candidates, primary_cache, output):
    """Complete source job; this function executes no theory or numerical callback."""
    manifest_path, cache, candidates, primary_cache = map(
        lambda p: Path(p).resolve(), (manifest_path, cache, candidates, primary_cache))
    root = manifest_path.parent.parent
    fingerprint = source_fingerprint(root)
    verified = verify_historical_intake(manifest_path, cache, candidates)
    contract_path, initial_path = root/CONTRACT, root/"config/initial-contract.json"
    contract, initial = read_json(contract_path), read_json(initial_path)
    primary_manifest = root/contract["primary_manifest"]
    if file_hash(primary_manifest) != contract["primary_manifest_sha256"]:
        raise ValueError("changed primary convention manifest")
    for path, digest in contract["parent_sha256"].items():
        if file_hash(root/path) != digest:
            raise ValueError("changed SIDIS parent contract")
    primary = read_json(primary_manifest)
    if (primary["id"] != "KPSY-HISTORICAL-SIDIS-PRIMARY-001" or len(primary["sources"]) != 9 or
            primary["production_ready"] is not False or primary["physical_extraction"] is not False):
        raise ValueError("wrong primary evidence scope")
    primary_inputs = {}
    for source in primary["sources"]:
        filename = source["filename"]
        if Path(filename).name != filename or not filename or filename in {".", ".."}:
            raise ValueError("unsafe primary filename")
        path = primary_cache/filename
        if not path.is_file() or file_hash(path) != source["sha256"]:
            raise ValueError("missing or changed primary convention bytes")
        primary_inputs[str(path)] = source["sha256"]
    if len(primary_inputs) != 9:
        raise ValueError("duplicate primary input")
    candidate_path, report_path = candidates/"candidates.json", candidates/"report.json"
    candidate, report = read_json(candidate_path), read_json(report_path)
    points = request_points(candidate, contract, initial)
    selected_ids = {v["id"] for v in contract["selection"]["source_views"]}
    views = {v["id"]: v for v in candidate["views"]}
    bound_paths = (manifest_path, candidate_path, report_path, contract_path, initial_path, primary_manifest)
    inputs = {**report["source_hashes"], **primary_inputs, **{str(p): file_hash(p) for p in bound_paths},
              **{str(root/p): digest for p, digest in contract["parent_sha256"].items()}}
    request = {"schema_version": 1, "scope": SCOPE, "theory_scope": "reconstructed_theory",
        "observation_scope": CANDIDATE_SCOPE, **{k: False for k in FALSE_FLAGS},
        "accepted_fitting_gate": False, "likelihood_constructed": False,
        "contract_id": contract["id"], "contract_sha256": file_hash(contract_path),
        "initial_contract_sha256": file_hash(initial_path), "primary_manifest_sha256": file_hash(primary_manifest),
        "source_fingerprint": fingerprint,
        "parameters": {"order": initial["parameter_order"], "values": initial["published_parameters"],
                       "role": initial["published_parameter_role"]},
        "model_conventions": copy.deepcopy(contract["model_conventions"]),
        "target_weights": copy.deepcopy(contract["prediction"]["target_weights"]),
        "view_evidence": {key: copy.deepcopy({k: v for k, v in views[key].items() if k != "points"})
                          for key in sorted(selected_ids)},
        "points": points, "source_inputs": inputs,
        "source_preparation": {"status": "CURRENT_HISTORICAL_SOURCE_REPLAY_VERIFIED",
            "verification_kind": "Python_source_preparation_attestation_not_Julia_raw_table_replay",
            "scope": verified["scope"], "historical_counts": verified["counts"],
            "manifest": str(manifest_path), "cache": str(cache), "candidates": str(candidates),
            "candidate_sha256": file_hash(candidate_path), "report_sha256": file_hash(report_path),
            "primary_cache": str(primary_cache), "primary_inputs": primary_inputs,
            "primary_verification": "pinned_convention_bytes_not_numerical_transcription_or_dataset_admission"},
        "limitations": contract["limitations"]}
    if source_fingerprint(root) != fingerprint or any(file_hash(p) != h for p, h in inputs.items()):
        raise ValueError("SIDIS sources changed during preparation")
    atomic_json(Path(output), request)
    return {"status": "FULL_SIDIS_REQUEST_PREPARED_NOT_PREDICTED", "scope": SCOPE, "points": len(points),
            "distinct_geometries": len(set(map(geometry_key, points))),
            "primary_files": len(primary_inputs), "request_sha256": file_hash(output),
            **{k: False for k in FALSE_FLAGS}}


def verify_request(path):
    """Fresh full-source replay, preserving the original request bytes."""
    path = Path(path)
    request = read_json(path)
    source = request["source_preparation"]
    with tempfile.TemporaryDirectory(prefix="historical-sidis-request-", dir=path.parent) as folder:
        fresh = Path(folder)/"request.json"
        prepare_request(source["manifest"], source["cache"], source["candidates"], source["primary_cache"], fresh)
        if read_json(fresh) != request:
            raise ValueError("changed source-derived full SIDIS request")
    return request
