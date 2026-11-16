"""Source-bound BaBar fixed-point jobs; no real-data likelihood admission."""
from __future__ import annotations

import copy
import csv
import os
from pathlib import Path
import subprocess
import tempfile

from .historical import SCOPE as CANDIDATE_SCOPE, verify_historical_intake
from .historical_points import finite_number
from .provenance import atomic_json, file_hash, git_commit, read_json, source_fingerprint, utc_now

SCOPE = "historical_reconstructed_sia_point_prediction_diagnostic"
CONTRACT = "config/historical-sia-point-contract.json"
FALSE_FLAGS = ("production_ready", "physical_extraction", "historical_fit_ready", "experimental_data_admitted",
               "fit_performed", "G2PERF_pass", "interval_coverage_verified")
MODEL_CONVENTION = {
    "frame": "hadron2", "harmonic": "cos_2phi0", "double_ratio_approximation": "linearized",
    "angular_definition": "mean_of_angular_ratio", "kinematic_prescription": "fixed_reported_means",
    "legacy_Q2_status": "source_literal_not_checked_by_selected_primary_table",
    "source_mapping_status": "conditional_fixed_mean_factorized_angular_model_not_acceptance_admission",
}


def request_points(candidate, contract, initial):
    """Check representation and selected source identities; preparation also replays all raw sources."""
    if (contract.get("id") != "KPSY-HISTORICAL-SIA-POINT-001" or contract.get("scope") != SCOPE or
            any(contract.get(k) is not False for k in FALSE_FLAGS) or candidate.get("scope") != CANDIDATE_SCOPE or
            any(candidate.get(k) is not False for k in ("production_ready", "physical_extraction", "historical_fit_ready"))):
        raise ValueError("SIA point diagnostic cannot approve fitting or substitute scope")
    if (len(initial["parameter_order"]) != 13 or len(set(initial["parameter_order"])) != 13 or
            len(initial["published_parameters"]) != 13):
        raise ValueError("complete fixed thirteen-coordinate parameter source required")
    for value in initial["published_parameters"]:
        finite_number(value, "published parameter")
    views = {view["id"]: view for view in candidate["views"]}
    if len(views) != len(candidate["views"]):
        raise ValueError("duplicate historical candidate view")
    selection, primary = contract["selection"], contract["primary_definition"]
    points = []
    for selected in selection["source_views"]:
        view = views.get(selected["id"])
        if (selected["ratio"] not in ("UL", "UC") or view is None or view["process"] != "SIA" or
                view["projection"] != "pT" or view["ratio"] != selected["ratio"] or view["frame"] != "hadron2" or
                view["publication"] != "1309.5278" or view["overlap_group"] != "babar_1309.5278" or
                view["source_id"] != selected["id"] or view["source_sha256"] != selected["sha256"] or
                view["source_observable"] != f"A{selected['ratio']}-0-PT" or view["input_units"] != "percent" or
                view["output_units"] != "fraction" or view["production_ready"] is not False or
                view["historical_fit_ready"] is not False or not view["blockers"] or not view["overlap_group"] or
                len(view["points"]) != selected["points"] or view["excluded"]):
            raise ValueError("missing or mismatched primary BaBar RF0 P-point selection")
        crosscheck = view.get("primary_numeric_crosscheck", {})
        if (crosscheck.get("publication") != selection["publication"] or crosscheck.get("frame") != "hadron2" or
                crosscheck.get("sha256") != primary["archive_sha256"] or crosscheck.get("member") != primary["member"] or
                crosscheck.get("member_sha256") != primary["member_sha256"] or crosscheck.get("table") != "tab:totPt" or
                crosscheck.get("raw_points") != selected["points"] or crosscheck.get("numeric_cells") != 9*selected["points"] or
                set(crosscheck.get("compared_columns", [])) != set(primary["checked_columns"])):
            raise ValueError("complete primary P-table crosscheck identity required")
        for index, row in enumerate(view["points"], 2):
            if row["id"] != f"{selected['id']}:r{index}" or row.get("locator", {}).get("excel_row") != index:
                raise ValueError("changed SIA source point identity or order")
            k, m, angular = row["averages"], row["systematic_components"], row["angular_source"]
            if set(k) != {"z1", "z2", "Q2_GeV2", "pT_GeV"}:
                raise ValueError("complete SIA reported mean kinematics required")
            for key in ("z1", "z2"):
                if not 0 < finite_number(k[key], key) < 1:
                    raise ValueError("invalid SIA momentum fraction")
            if finite_number(k["Q2_GeV2"], "Q2", positive=True) != 112.36:
                raise ValueError("preserve the literal conditional legacy Q2")
            finite_number(k["pT_GeV"], "pT", nonnegative=True)
            finite_number(row["central_value"], "central value")
            finite_number(row["statistical_error"], "statistical error", positive=True)
            if (set(m) != {"reported_systematic"} or
                    row.get("uncertainty_source_columns") != {"statistical": "stat_u", "reported_systematic": "sys_u"}):
                raise ValueError("preserve the separate reported systematic error")
            finite_number(m["reported_systematic"], "reported systematic", nonnegative=True)
            if (angular.get("column") != "S2/1+C2" or angular.get("average_type") != "mean_of_angular_ratio" or
                    angular.get("frame") != "hadron2" or angular.get("measure_status") != "joint_acceptance_unresolved" or
                    not angular.get("averaging_measure") or not 0 <= finite_number(angular["value"], "angular scalar") <= 1):
                raise ValueError("reviewed mean-of-angular-ratio and unresolved joint measure required")
            locator, bins = row.get("primary_locator", {}), row.get("primary_bins", {})
            if (locator.get("publication") != selection["publication"] or locator.get("member") != primary["member"] or
                    locator.get("table") != "tab:totPt" or locator.get("table_row") != index-1 or
                    not isinstance(locator.get("source_line"), int) or isinstance(locator["source_line"], bool) or
                    locator["source_line"] <= 0 or set(bins) != {"pT"} or
                    row.get("bins", {}).get("pT_GeV") != bins["pT"] or len(bins["pT"]) != 2):
                raise ValueError("primary RF0 point/bin locator required")
            low, high = [finite_number(x, "primary P bound", nonnegative=True) for x in bins["pT"]]
            if not low <= k["pT_GeV"] <= high or low >= high:
                raise ValueError("inconsistent source P mean and bounds")
            points.append({"id": row["id"], "dataset_id": view["id"], "publication": view["publication"],
                "primary_publication": crosscheck["publication"], "source_sha256": view["source_sha256"],
                "ratio": view["ratio"], "frame": view["frame"], "harmonic": "cos_2phi0",
                "kinematics": copy.deepcopy(k), "measured": {key: copy.deepcopy(row[key]) for key in
                    ("central_value", "statistical_error", "systematic_components")},
                "source_locator": copy.deepcopy(row["locator"]), "primary_locator": copy.deepcopy(locator),
                "source_bins": copy.deepcopy(row["bins"]), "primary_bins": copy.deepcopy(bins),
                "source_angular": copy.deepcopy(angular), "source_primary_crosscheck": copy.deepcopy(crosscheck),
                "error_source_columns": copy.deepcopy(row["uncertainty_source_columns"]),
                "overlap_group": view["overlap_group"], "source_blockers": list(view["blockers"]),
                "legacy_Q2_status": MODEL_CONVENTION["legacy_Q2_status"]})
    geometries = {tuple(p["kinematics"][k] for k in ("z1", "z2", "Q2_GeV2", "pT_GeV")) for p in points}
    if (len(points) != selection["total_points"] or len({p["id"] for p in points}) != len(points) or
            len(geometries) != selection["exact_shared_geometry_count"]):
        raise ValueError("incomplete or mismatched exact SIA point geometry selection")
    return points


def prepare_request(manifest_path, cache, candidates, output):
    """Rebuild all source preparation before handing fixed source rows to Julia."""
    output = Path(output)
    if output.exists():
        raise ValueError("use a fresh SIA request output; preserve earlier source evidence")
    manifest_path, cache, candidates = [Path(p).resolve() for p in (manifest_path, cache, candidates)]
    root = manifest_path.parent.parent
    fingerprint = source_fingerprint(root)
    verified = verify_historical_intake(manifest_path, cache, candidates)
    candidate_path, report_path = candidates/"candidates.json", candidates/"report.json"
    candidate, report = read_json(candidate_path), read_json(report_path)
    contract_path, initial_path = root/CONTRACT, root/"config/initial-contract.json"
    contract, initial = read_json(contract_path), read_json(initial_path)
    points = request_points(candidate, contract, initial)
    if report["primary_BaBar_raw_points_crosschecked"] != 90 or report["primary_BaBar_numeric_cells_crosschecked"] != 954:
        raise ValueError("incomplete existing BaBar primary source replay")
    inputs = {**report["source_hashes"], **{str(p): file_hash(p) for p in
        (candidate_path, report_path, contract_path, initial_path)}}
    request = {"schema_version": 1, "scope": SCOPE, "theory_scope": "reconstructed_theory",
        "observation_scope": CANDIDATE_SCOPE, **{key: False for key in FALSE_FLAGS},
        "contract_id": contract["id"], "contract_sha256": file_hash(contract_path),
        "initial_contract_sha256": file_hash(initial_path), "source_fingerprint": fingerprint,
        "parameters": {"order": initial["parameter_order"], "values": initial["published_parameters"],
                       "role": initial["published_parameter_role"]},
        "model_convention": dict(MODEL_CONVENTION), "points": points, "source_inputs": inputs,
        "source_preparation": {"status": "CURRENT_HISTORICAL_SOURCE_REPLAY_VERIFIED",
            "verification_kind": "Python_source_preparation_attestation_not_Julia_raw_table_replay",
            "scope": verified["scope"], "historical_counts": verified["counts"],
            "manifest": str(manifest_path), "cache": str(cache), "candidates": str(candidates),
            "candidate_sha256": file_hash(candidate_path), "report_sha256": file_hash(report_path)},
        "limitations": contract["limitations"]}
    if source_fingerprint(root) != fingerprint or any(file_hash(p) != h for p, h in inputs.items()):
        raise ValueError("SIA sources changed during request preparation")
    atomic_json(output, request)
    return {"status": "PREDICTION_REQUEST_PREPARED", "scope": SCOPE, "points": len(points),
            **{key: False for key in FALSE_FLAGS}, "request_sha256": file_hash(output)}


def verify_request(path):
    """Reconstruct source-derived values rather than trusting a rehashed request."""
    path = Path(path)
    request = read_json(path)
    source = request["source_preparation"]
    with tempfile.TemporaryDirectory(prefix="historical-sia-request-", dir=path.parent) as folder:
        fresh = Path(folder)/"request.json"
        prepare_request(source["manifest"], source["cache"], source["candidates"], fresh)
        if read_json(fresh) != request:
            raise ValueError("changed source-derived SIA request")
    return request


def render_comparison(report_path, output):
    """Render completed Julia fields; no Python prediction, fitting or error combination."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    report = read_json(report_path)
    if report["scope"] != SCOPE or report["status"] != "PASS" or any(report.get(k) is not False for k in FALSE_FLAGS):
        raise ValueError("only a completed scoped SIA diagnostic can render comparison artifacts")
    output = Path(output)
    table = output/"point-comparison.csv"
    with table.open("w", encoding="utf-8", newline="") as f:
        f.write("# scope="+SCOPE+"; fixed rounded parameters; no fit or likelihood\n")
        writer = csv.writer(f)
        writer.writerow(["point_id", "ratio", "z1", "z2", "Q2_GeV2", "pT_GeV", "angular_mean_of_ratio",
                         "measured", "statistical_error", "reported_systematic", "fixed_theory",
                         "unweighted_observed_minus_predicted"])
        for p in report["points"]:
            k, m = p["kinematics"], p["measured"]
            writer.writerow([p["id"], p["ratio"], *[k[x] for x in ("z1", "z2", "Q2_GeV2", "pT_GeV")],
                p["source_angular"]["value"], m["central_value"], m["statistical_error"],
                m["systematic_components"]["reported_systematic"], p["amplitude"], p["unweighted_observed_minus_predicted"]])
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), sharex=True)
    for ax, ratio in zip(axes, ("UL", "UC")):
        rows = [p for p in report["points"] if p["ratio"] == ratio]
        xs = [p["kinematics"]["pT_GeV"] for p in rows]
        ax.errorbar(xs, [p["measured"]["central_value"] for p in rows],
            yerr=[p["measured"]["statistical_error"] for p in rows], fmt="o", capsize=3,
            color="#1b5a85", label="BaBar RF0: reported statistical error")
        ax.scatter(xs, [p["amplitude"] for p in rows], marker="D", facecolors="none", edgecolors="#c25b12",
                   s=45, zorder=4, label="Fixed rounded KPSY parameters")
        ax.axhline(0, color="0.55", linewidth=.7)
        ax.set(title=ratio, xlabel=r"Reported $\langle P_{h\perp}\rangle$ [GeV]", ylabel="Cosine amplitude (fraction)")
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=8, loc="best")
    fig.suptitle("BaBar RF0 fixed-mean comparison\nReconstructed historical diagnostic", fontsize=12)
    fig.text(.5, .025, "Conditional factorized angular/linearized model; joint acceptance and legacy Q precision remain unresolved.\n"
        "Statistical bars shown; systematic errors remain separate in the table. No fit or likelihood.",
        ha="center", va="bottom", fontsize=8)
    fig.tight_layout(rect=(0, .15, 1, .90))
    figure = output/"point-comparison.png"
    fig.savefig(figure, dpi=180)
    plt.close(fig)
    return [table, figure]


def _julia_environment():
    return {**os.environ, **{key: "1" for key in ("JULIA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")},
            "MPLBACKEND": "Agg"}


def _julia_command(julia, root, program, *arguments):
    return [julia, "--startup-file=no", "--project="+str(root/"numerics"),
            str(root/"numerics/bin"/program), *map(str, arguments)]


def run_point_job(request_path, native_build, parton_build, output, *, julia="julia"):
    """Schedule complete Julia prediction/audit jobs, then render their outputs."""
    request_path, native_build, parton_build = [Path(p).resolve() for p in (request_path, native_build, parton_build)]
    request = verify_request(request_path)
    root = Path(request["source_preparation"]["manifest"]).parent.parent
    output = Path(output).resolve()
    if output.is_dir() and any(output.iterdir()):
        raise ValueError("use a fresh SIA diagnostic output directory; preserve earlier artifacts")
    output.mkdir(parents=True, exist_ok=True)
    report_path, audit_path, invocation_path = [output/name for name in ("predictions.json", "audit.json", "invocation.json")]
    if any(path.exists() for path in (report_path, audit_path, invocation_path)):
        raise ValueError("use a fresh historical prediction output directory")
    fingerprint = source_fingerprint(root)
    hashes = {str(p): file_hash(p) for p in (request_path, native_build, parton_build)}
    invocation = {"status": "FAIL", "scope": SCOPE, **{key: False for key in FALSE_FLAGS},
        "G2PERF_pass": False, "source_fingerprint": fingerprint, "code_commit": git_commit(root),
        "input_hashes": hashes, "steps": [], "outputs": {}}
    commands = [("predictions", _julia_command(julia, root, "predict_historical_sia_points.jl",
                   request_path, native_build, parton_build, report_path)),
                ("audit", _julia_command(julia, root, "audit_historical_sia_points.jl", request_path, report_path, audit_path))]
    try:
        for name, command in commands:
            log = output/(name+".log")
            with log.open("w", encoding="utf-8", newline="\n") as stream:
                completed = subprocess.run(command, cwd=root, env=_julia_environment(),
                                           stdout=stream, stderr=subprocess.STDOUT)
            invocation["steps"].append({"name": name, "command": command, "returncode": completed.returncode,
                                        "log_sha256": file_hash(log)})
            atomic_json(invocation_path, invocation)
            if completed.returncode != 0:
                raise ValueError("historical point "+name+" failed; preserve report and log: "+str(log))
        audit = read_json(audit_path)
        if (audit["status"] != "PASS" or audit["accepted_prediction_diagnostic"] is not True or
                audit["accepted_fitting_gate"] is not False):
            raise ValueError("unaccepted historical diagnostic audit")
        products = render_comparison(report_path, output)
        verify_request(request_path)
        if source_fingerprint(root) != fingerprint or any(file_hash(p) != sha for p, sha in hashes.items()):
            raise ValueError("historical diagnostic inputs/code changed during complete job")
        invocation["outputs"] = {str(p): file_hash(p) for p in (report_path, audit_path, *products)}
        invocation["status"] = "PASS"
    except Exception as exc:
        invocation["error"] = type(exc).__name__+": "+str(exc)
        raise
    finally:
        invocation["completed_at"] = utc_now()
        atomic_json(invocation_path, invocation)
    return {"status": invocation["status"], "scope": SCOPE, "points": len(request["points"]),
            **{key: False for key in FALSE_FLAGS}, "G2PERF_pass": False, "output": str(output)}


def verify_point_job(request_path, output, *, julia="julia"):
    """Check cached bindings, replay source preparation and audit stored Julia operands."""
    request_path, output = Path(request_path).resolve(), Path(output).resolve()
    request = verify_request(request_path)
    root = Path(request["source_preparation"]["manifest"]).parent.parent
    invocation = read_json(output/"invocation.json")
    if (invocation.get("status") != "PASS" or invocation.get("scope") != SCOPE or
            any(invocation.get(k) is not False for k in FALSE_FLAGS) or
            invocation.get("G2PERF_pass") is not False or invocation["source_fingerprint"] != source_fingerprint(root)):
        raise ValueError("stale or wrong-scope cached historical diagnostic")
    for path, sha in invocation["input_hashes"].items():
        if file_hash(path) != sha:
            raise ValueError("changed historical diagnostic input")
    if invocation["input_hashes"].get(str(request_path)) != file_hash(request_path):
        raise ValueError("cached historical request identity changed")
    expected = {str(output/name) for name in ("predictions.json", "audit.json", "point-comparison.csv", "point-comparison.png")}
    if set(invocation["outputs"]) != expected:
        raise ValueError("incomplete cached diagnostic output bindings")
    for path, sha in invocation["outputs"].items():
        if file_hash(path) != sha:
            raise ValueError("changed historical diagnostic output")
    original_audit = read_json(output/"audit.json")
    if (original_audit.get("status") != "PASS" or original_audit.get("scope") != SCOPE or
            original_audit.get("accepted_prediction_diagnostic") is not True or
            original_audit.get("accepted_fitting_gate") is not False or
            any(original_audit.get(k) is not False for k in FALSE_FLAGS) or
            original_audit.get("report_sha256") != file_hash(output/"predictions.json") or
            original_audit.get("request_sha256") != file_hash(request_path)):
        raise ValueError("invalid original diagnostic audit metadata")
    for step in invocation["steps"]:
        if step["returncode"] != 0 or file_hash(output/(step["name"]+".log")) != step["log_sha256"]:
            raise ValueError("failed or changed historical diagnostic step")
    if [s["name"] for s in invocation["steps"]] != ["predictions", "audit"]:
        raise ValueError("incomplete cached diagnostic step ledger")
    with tempfile.TemporaryDirectory(prefix="historical-numerical-audit-", dir=output.parent) as folder:
        audit_path = Path(folder)/"audit.json"
        command = _julia_command(julia, root, "audit_historical_sia_points.jl", request_path, output/"predictions.json", audit_path)
        completed = subprocess.run(command, cwd=root, env=_julia_environment(), capture_output=True)
        if completed.returncode != 0:
            failure = output/("cache-audit-failure-"+file_hash(output/"predictions.json")[:16]+".json")
            if not failure.exists():
                atomic_json(failure, {"status": "FAIL", "scope": SCOPE, "production_ready": False,
                    "report_sha256": file_hash(output/"predictions.json"),
                    "audit": read_json(audit_path) if audit_path.exists() else None,
                    "log": completed.stdout.decode("utf-8", errors="replace")+completed.stderr.decode("utf-8", errors="replace")})
            raise ValueError("cached diagnostic numerical audit rejected report")
        audit = read_json(audit_path)
        if audit["status"] != "PASS" or audit["accepted_fitting_gate"] is not False:
            raise ValueError("cached diagnostic cannot authorize fitting")
    return {"status": "CACHED_DIAGNOSTIC_REVALIDATED", "scope": SCOPE, **{key: False for key in FALSE_FLAGS},
            "G2PERF_pass": False, "points": len(request["points"])}
