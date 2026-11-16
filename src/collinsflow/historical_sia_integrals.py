"""Source-bound conditional BaBar momentum integrals; no likelihood admission."""
from __future__ import annotations

import copy
import csv
import hashlib
import os
from pathlib import Path
import subprocess
import tempfile
import tarfile
import time

from .historical import SCOPE as CANDIDATE_SCOPE, verify_historical_intake
from .historical_points import finite_number
from .job_control import run_complete_process
from .provenance import atomic_json, canonical, file_hash, git_commit, read_json, source_fingerprint, utc_now

SCOPE = "historical_reconstructed_sia_integral_prediction_diagnostic"
CONTRACT = "config/historical-sia-integral-contract.json"
FALSE_FLAGS = ("production_ready", "physical_extraction", "historical_fit_ready", "experimental_data_admitted",
               "fit_performed", "G2PERF_pass", "interval_coverage_verified")
MODEL_CONVENTION = {
    "frame": "hadron2", "harmonic": "cos_2phi0", "double_ratio_approximation": "linearized",
    "angular_definition": "mean_of_angular_ratio", "kinematic_prescription": "fixed_reported_z_means_finite_qT_integral",
    "legacy_Q2_status": "source_literal_not_checked_by_selected_primary_table",
    "source_mapping_status": "conditional_fixed_z_uniform_momentum_factorized_angular_model_not_acceptance_admission",
}


def momentum_window(contract):
    """Explicit model settings, not a recovered experimental acceptance."""
    prediction = contract["prediction"]
    for bound in prediction["qT_window_GeV"]:
        finite_number(bound, "fixed qT bound", nonnegative=True)
    if prediction["qT_window_GeV"] != [0., 3.5] or not prediction["lower_limit_source"]:
        raise ValueError("preserve the declared conditional qT window")
    return {"coordinate": "qT=P_hperp/z1 at fixed theory z1", "qT_bounds_GeV": [0., 3.5],
            "measure": "P*dP=z1^2*qT*dqT", "source": contract["id"]+":conditional_fixed_z_qT_window",
            "lower_limit_source": prediction["lower_limit_source"],
            "acceptance": "uniform momentum acceptance; explicit reconstruction assumption only"}


def verify_primary_upper_cut(report, contract):
    primary = contract["primary_definition"]
    files = [Path(path) for path, digest in report["source_hashes"].items() if digest == primary["archive_sha256"]]
    if len(files) != 1 or file_hash(files[0]) != primary["archive_sha256"]:
        raise ValueError("missing pinned primary upper-cut source")
    with tarfile.open(files[0]) as archive:
        stream = archive.extractfile(primary["member"])
        if stream is None:
            raise ValueError("missing primary TeX member")
        raw = stream.read()
    if hashlib.sha256(raw).hexdigest() != primary["member_sha256"]:
        raise ValueError("changed primary upper-cut member")
    line = primary["upper_cut"]["source_line"]
    if raw.decode("utf-8").splitlines()[line-1].strip() != r"and a safe compromise is to require $Q_t<3.5$ \gevc.":
        raise ValueError("changed source upper-cut literal")
    return {"archive_sha256": primary["archive_sha256"], "member_sha256": primary["member_sha256"],
            "source_line": line, "qT_upper_GeV": 3.5, "definition": "event_cut_not_a_mean_ratio_cut",
            "experimental_lower_acceptance_verified": False}


def request_points(candidate, contract, initial):
    """Check representation and selected source identities; preparation also replays all raw sources."""
    if (contract.get("id") != "KPSY-HISTORICAL-SIA-INTEGRAL-001" or contract.get("scope") != SCOPE or
            any(contract.get(k) is not False for k in FALSE_FLAGS) or candidate.get("scope") != CANDIDATE_SCOPE or
            any(candidate.get(k) is not False for k in ("production_ready", "physical_extraction", "historical_fit_ready"))):
        raise ValueError("SIA point diagnostic cannot approve fitting or substitute scope")
    if (len(initial["parameter_order"]) != 13 or len(set(initial["parameter_order"])) != 13 or
            len(initial["published_parameters"]) != 13):
        raise ValueError("complete fixed thirteen-coordinate parameter source required")
    for value in initial["published_parameters"]:
        finite_number(value, "published parameter")
    momentum_window(contract)
    views = {view["id"]: view for view in candidate["views"]}
    if len(views) != len(candidate["views"]):
        raise ValueError("duplicate historical candidate view")
    selection, primary = contract["selection"], contract["primary_definition"]
    points = []
    for selected in selection["source_views"]:
        view = views.get(selected["id"])
        if (selected["ratio"] not in ("UL", "UC") or view is None or view["process"] != "SIA" or
                view["projection"] != "z1_z2" or view["ratio"] != selected["ratio"] or view["frame"] != "hadron2" or
                view["publication"] != "1309.5278" or view["overlap_group"] != "babar_1309.5278" or
                view["source_id"] != selected["id"] or view["source_sha256"] != selected["sha256"] or
                view["source_observable"] != f"A{selected['ratio']}-0-PT-INT" or
                view["prediction_prescription"] != "P-integrated source view; theory adapter pending" or view["input_units"] != "percent" or
                view["output_units"] != "fraction" or view["production_ready"] is not False or
                view["historical_fit_ready"] is not False or not view["blockers"] or not view["overlap_group"] or
                len(view["points"]) != selected["points"] or view["excluded"]):
            raise ValueError("missing or mismatched primary BaBar RF0 integral selection")
        crosscheck = view.get("primary_numeric_crosscheck", {})
        systematic_column = "syt_u" if selected["id"] == "KPSY15_SIA_1003" else "sys_u"
        columns = {"S2/1+C2", "pT", "stat_u", systematic_column, "value", "z1", "z1min", "z1max", "z2", "z2min", "z2max"}
        if (crosscheck.get("publication") != selection["publication"] or crosscheck.get("frame") != "hadron2" or
                crosscheck.get("sha256") != primary["archive_sha256"] or crosscheck.get("member") != primary["member"] or
                crosscheck.get("member_sha256") != primary["member_sha256"] or crosscheck.get("table") != selected["primary_table"] or
                crosscheck.get("raw_points") != selected["points"] or crosscheck.get("numeric_cells") != 11*selected["points"] or
                set(crosscheck.get("compared_columns", [])) != columns):
            raise ValueError("complete primary z-table crosscheck identity required")
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
                    row.get("uncertainty_source_columns") != {"statistical": "stat_u", "reported_systematic": systematic_column}):
                raise ValueError("preserve the separate reported systematic error")
            finite_number(m["reported_systematic"], "reported systematic", nonnegative=True)
            if (angular.get("column") != "S2/1+C2" or angular.get("average_type") != "mean_of_angular_ratio" or
                    angular.get("frame") != "hadron2" or angular.get("measure_status") != "joint_acceptance_unresolved" or
                    not angular.get("averaging_measure") or not 0 <= finite_number(angular["value"], "angular scalar") <= 1):
                raise ValueError("reviewed mean-of-angular-ratio and unresolved joint measure required")
            locator, bins = row.get("primary_locator", {}), row.get("primary_bins", {})
            if (locator.get("publication") != selection["publication"] or locator.get("member") != primary["member"] or
                    locator.get("table") != selected["primary_table"] or locator.get("table_row") != index-1 or
                    not isinstance(locator.get("source_line"), int) or isinstance(locator["source_line"], bool) or
                    locator["source_line"] <= 0 or set(bins) != {"z1", "z2"} or row.get("bins") != bins):
                raise ValueError("primary RF0 z-bin locator required")
            for key in ("z1", "z2"):
                low, high = [finite_number(x, "primary z bound", positive=True) for x in bins[key]]
                if not 0 < low <= k[key] <= high <= 1 or low >= high:
                    raise ValueError("inconsistent source z mean and bounds")
            points.append({"id": row["id"], "dataset_id": view["id"], "publication": view["publication"],
                "primary_publication": crosscheck["publication"], "source_sha256": view["source_sha256"],
                "ratio": view["ratio"], "frame": view["frame"], "harmonic": "cos_2phi0",
                "kinematics": {key: k[key] for key in ("z1", "z2", "Q2_GeV2")},
                "reported_pT_mean_GeV": k["pT_GeV"], "reported_pT_role": "source_metadata_not_model_momentum_or_bound",
                "mapped_P_bounds_GeV": [0., k["z1"]*3.5], "measured": {key: copy.deepcopy(row[key]) for key in
                    ("central_value", "statistical_error", "systematic_components")},
                "source_locator": copy.deepcopy(row["locator"]), "primary_locator": copy.deepcopy(locator),
                "source_bins": copy.deepcopy(row["bins"]), "primary_bins": copy.deepcopy(bins),
                "source_angular": copy.deepcopy(angular), "source_primary_crosscheck": copy.deepcopy(crosscheck),
                "error_source_columns": copy.deepcopy(row["uncertainty_source_columns"]),
                "overlap_group": view["overlap_group"], "source_blockers": list(view["blockers"]),
                "legacy_Q2_status": MODEL_CONVENTION["legacy_Q2_status"]})
    geometries = {tuple(p["kinematics"][k] for k in ("z1", "z2", "Q2_GeV2")) for p in points}
    if (len(points) != selection["total_points"] or len({p["id"] for p in points}) != len(points) or
            len(geometries) != selection["exact_geometry_count"]):
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
    window_path = root/"config/sia-qt-window-contract.json"
    contract, initial = read_json(contract_path), read_json(initial_path)
    points = request_points(candidate, contract, initial)
    upper_cut = verify_primary_upper_cut(report, contract)
    if report["primary_BaBar_raw_points_crosschecked"] != 90 or report["primary_BaBar_numeric_cells_crosschecked"] != 954:
        raise ValueError("incomplete existing BaBar primary source replay")
    inputs = {**report["source_hashes"], **{str(p): file_hash(p) for p in
        (candidate_path, report_path, contract_path, initial_path, window_path)}}
    request = {"schema_version": 1, "scope": SCOPE, "theory_scope": "reconstructed_theory",
        "observation_scope": CANDIDATE_SCOPE, **{key: False for key in FALSE_FLAGS},
        "contract_id": contract["id"], "contract_sha256": file_hash(contract_path),
        "initial_contract_sha256": file_hash(initial_path), "source_fingerprint": fingerprint,
        "window_contract_sha256": file_hash(window_path), "momentum_window": momentum_window(contract),
        "parameters": {"order": initial["parameter_order"], "values": initial["published_parameters"],
                       "role": initial["published_parameter_role"]},
        "model_convention": dict(MODEL_CONVENTION), "points": points, "source_inputs": inputs,
        "source_preparation": {"status": "CURRENT_HISTORICAL_SOURCE_REPLAY_VERIFIED",
            "verification_kind": "Python_source_preparation_attestation_not_Julia_raw_table_replay",
            "scope": verified["scope"], "historical_counts": verified["counts"],
            "manifest": str(manifest_path), "cache": str(cache), "candidates": str(candidates),
            "candidate_sha256": file_hash(candidate_path), "report_sha256": file_hash(report_path),
            "primary_upper_cut_attestation": upper_cut},
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
        if canonical(read_json(fresh)) != canonical(request):
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
        writer.writerow(["point_id", "ratio", "z1", "z2", "Q2_GeV2", "reported_pT_mean_GeV",
                         "qT_lower_GeV", "qT_upper_GeV", "P_lower_GeV", "P_upper_GeV", "angular_mean_of_ratio",
                         "measured", "statistical_error", "reported_systematic", "fixed_theory",
                         "unweighted_observed_minus_predicted"])
        for p in report["points"]:
            k, m = p["kinematics"], p["measured"]
            writer.writerow([p["id"], p["ratio"], *[k[x] for x in ("z1", "z2", "Q2_GeV2")],
                p["reported_pT_mean_GeV"], *report["momentum_window"]["qT_bounds_GeV"], *p["mapped_P_bounds_GeV"],
                p["source_angular"]["value"], m["central_value"], m["statistical_error"],
                m["systematic_components"]["reported_systematic"], p["amplitude"], p["unweighted_observed_minus_predicted"]])
    bins = sorted({tuple(p["primary_bins"]["z1"]) for p in report["points"]})
    fig, axes = plt.subplots(len(bins), 2, figsize=(11, 2.6*len(bins)), sharex=True, squeeze=False)
    for row_index, interval in enumerate(bins):
        for col, ratio in enumerate(("UL", "UC")):
            ax = axes[row_index, col]
            rows = [p for p in report["points"] if p["ratio"] == ratio and tuple(p["primary_bins"]["z1"]) == interval]
            xs = [p["kinematics"]["z2"] for p in rows]
            ax.errorbar(xs, [p["measured"]["central_value"] for p in rows],
                yerr=[p["measured"]["statistical_error"] for p in rows], fmt="o", capsize=3,
                color="#1b5a85", label="BaBar RF0: statistical error")
            ax.scatter(xs, [p["amplitude"] for p in rows], marker="D", facecolors="none", edgecolors="#c25b12",
                       s=40, zorder=4, label="Fixed conditional integral")
            ax.axhline(0, color="0.55", linewidth=.7)
            ax.set(title=f"{ratio}; z1 source bin [{interval[0]}, {interval[1]}]", ylabel="Cosine amplitude (fraction)")
            ax.grid(alpha=.2)
    for ax in axes[-1]:
        ax.set_xlabel("Reported mean z2")
    axes[0, 0].legend(fontsize=8, loc="best")
    fig.suptitle("BaBar RF0 conditional momentum integrals\nFixed rounded KPSY parameters; reconstructed diagnostic", fontsize=12)
    fig.text(.5, .01, "Declared qT window [0,3.5] GeV; uniform momentum acceptance and fixed z/angular factorization.\n"
        "Joint acceptance and legacy Q precision remain unresolved. Reported P means are metadata.\n"
        "Statistical bars shown; systematic errors remain separate in the table. No fit or likelihood.",
        ha="center", va="bottom", fontsize=8)
    fig.tight_layout(rect=(0, .055, 1, .96))
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


def run_integral_job(request_path, native_build, parton_build, output, *, julia="julia", timeout_seconds=None):
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
    if timeout_seconds is not None:
        finite_number(timeout_seconds, "complete Julia-job time budget", positive=True)
    invocation["Julia_jobs_time_budget_seconds"] = timeout_seconds
    deadline = None if timeout_seconds is None else time.monotonic()+timeout_seconds
    atomic_json(invocation_path, invocation)
    commands = [("predictions", _julia_command(julia, root, "predict_historical_sia_integrals.jl",
                   request_path, native_build, parton_build, report_path)),
                ("audit", _julia_command(julia, root, "audit_historical_sia_integrals.jl", request_path, report_path, audit_path))]
    try:
        for name, command in commands:
            remaining = None if deadline is None else deadline-time.monotonic()
            if remaining is not None and remaining <= 0:
                raise ValueError("complete Julia-job resource budget exhausted before "+name)
            log = output/(name+".log")
            with log.open("w", encoding="utf-8", newline="\n") as stream:
                completed = run_complete_process(command, cwd=root, env=_julia_environment(),
                                                stdout=stream, timeout_seconds=remaining)
            invocation["steps"].append({"name": name, "command": command, "returncode": completed.returncode,
                "process_returncode": completed.process_returncode, "interruption": completed.interruption,
                "termination_signal": completed.termination_signal, "log_sha256": file_hash(log)})
            atomic_json(invocation_path, invocation)
            if completed.returncode != 0:
                raise ValueError("historical integral "+name+" failed; preserve report and log: "+str(log))
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
    except BaseException as exc:
        invocation["error"] = type(exc).__name__+": "+str(exc)
        raise
    finally:
        invocation["completed_at"] = utc_now()
        atomic_json(invocation_path, invocation)
    return {"status": invocation["status"], "scope": SCOPE, "points": len(request["points"]),
            **{key: False for key in FALSE_FLAGS}, "G2PERF_pass": False, "output": str(output)}


def verify_integral_job(request_path, output, *, julia="julia"):
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
        if (step["returncode"] != 0 or step.get("process_returncode") != 0 or
                step.get("interruption") is not None or step.get("termination_signal") is not None or
                file_hash(output/(step["name"]+".log")) != step["log_sha256"]):
            raise ValueError("failed or changed historical diagnostic step")
    if [s["name"] for s in invocation["steps"]] != ["predictions", "audit"]:
        raise ValueError("incomplete cached diagnostic step ledger")
    with tempfile.TemporaryDirectory(prefix="historical-numerical-audit-", dir=output.parent) as folder:
        audit_path = Path(folder)/"audit.json"
        command = _julia_command(julia, root, "audit_historical_sia_integrals.jl", request_path, output/"predictions.json", audit_path)
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
