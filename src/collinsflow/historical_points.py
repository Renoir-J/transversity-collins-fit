"""Prepare source-bound real-row prediction requests; never approve a likelihood."""
from __future__ import annotations

import copy
import csv
import math
import os
from pathlib import Path
import subprocess
import tempfile

from .historical import SCOPE as CANDIDATE_SCOPE, verify_historical_intake
from .provenance import atomic_json, file_hash, git_commit, read_json, source_fingerprint, utc_now

SCOPE = "historical_reconstructed_point_prediction_diagnostic"
CONTRACT = "config/historical-point-contract.json"
FALSE_FLAGS = ("production_ready", "physical_extraction", "historical_fit_ready", "experimental_data_admitted", "fit_performed")


def finite_number(value, name, *, positive=False, nonnegative=False):
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or
            (positive and value <= 0) or (nonnegative and value < 0)):
        raise ValueError("invalid historical point "+name)
    return value


def request_points(candidate, contract, initial):
    """Validate representation/selection; source replay belongs to prepare_request."""
    if (contract.get("id") != "KPSY-HISTORICAL-POINT-001" or contract.get("scope") != SCOPE or
            any(contract.get(k) is not False for k in FALSE_FLAGS) or
            candidate.get("scope") != CANDIDATE_SCOPE or
            any(candidate.get(k) is not False for k in ("production_ready", "physical_extraction", "historical_fit_ready"))):
        raise ValueError("historical prediction request cannot approve fitting or substitute scope")
    convention = contract["prediction"]
    if (convention["phase"] != "sin_phi_sum" or convention["normalization"] != "fixed_y_DNN_included" or
            convention["source_mapping_status"] != "conditional_historical_assumption_not_experimental_average_admission"):
        raise ValueError("changed conditional historical model convention")
    if len(initial["parameter_order"]) != len(initial["published_parameters"]) or len(initial["parameter_order"]) != 13:
        raise ValueError("complete fixed thirteen-coordinate parameter source required")
    if len(set(initial["parameter_order"])) != 13:
        raise ValueError("duplicate historical parameter coordinate")
    for value in initial["published_parameters"]:
        finite_number(value, "published parameter")
    selection = contract["selection"]
    views = {view["id"]: view for view in candidate["views"]}
    if len(views) != len(candidate["views"]):
        raise ValueError("duplicate historical candidate view")
    points = []
    ids = set()
    for selected in selection["source_views"]:
        view = views.get(selected["id"])
        if (view is None or view["process"] != "SIDIS" or view["projection"] != "x" or
                view["publication"] != selection["publication"] or view["source_id"] != selected["id"] or
                view["target"] != selection["target"] or view.get("target_material") != "3He" or
                view.get("target_result") != selection["target_result"] or view["hadron"] != selected["hadron"] or
                selected["charge"] not in (-1, 1) or isinstance(selected["charge"], bool) or
                {1: "pi+", -1: "pi-"}[selected["charge"]] != selected["hadron"] or
                view["source_sha256"] != selected["sha256"] or view["source_observable"] != "AUTcollins" or
                view["output_units"] != "fraction" or view.get("depolarization_removed") is not None or
                view["production_ready"] is not False or view["historical_fit_ready"] is not False or
                not view["blockers"] or len(view["points"]) != selected["points"] or view["excluded"]):
            raise ValueError("missing or mismatched primary HallA point selection")
        for row_index, row in enumerate(view["points"], 2):
            if (row["id"] in ids or row["id"] != f"{selected['id']}:r{row_index}" or
                    row.get("locator", {}).get("excel_row") != row_index):
                raise ValueError("duplicate or mismatched historical point identifier")
            ids.add(row["id"])
            averages = row["averages"]
            if set(averages) != {"x", "y", "z", "Q2_GeV2", "pT_GeV"}:
                raise ValueError("complete reported mean kinematics required")
            for key in ("x", "y", "z"):
                value = finite_number(averages[key], key)
                if not 0 < value < 1:
                    raise ValueError("invalid historical point fraction")
            finite_number(averages["Q2_GeV2"], "Q2", positive=True)
            finite_number(averages["pT_GeV"], "pT", nonnegative=True)
            finite_number(row["central_value"], "central value")
            finite_number(row["statistical_error"], "statistical error", positive=True)
            if set(row["systematic_components"]) != {"experimental", "model"}:
                raise ValueError("separate experimental and model source errors required")
            for key, value in row["systematic_components"].items():
                finite_number(value, key+" systematic", nonnegative=True)
            remap = row.get("source_remap", {})
            if (remap.get("kind") != "reviewed_legacy_column_remap" or
                    remap.get("primary_tables") != ["table:kine", "table:results"] or
                    len(remap.get("primary_lines", [])) != 2):
                raise ValueError("primary HallA remap provenance required")
            points.append({"id": row["id"], "dataset_id": selected["id"],
                "publication": view["publication"], "source_sha256": view["source_sha256"],
                "target": view["target"], "target_material": view["target_material"],
                "target_result": view["target_result"], "hadron": view["hadron"], "charge": selected["charge"],
                "kinematics": copy.deepcopy(averages),
                "measured": {key: copy.deepcopy(row[key]) for key in
                    ("central_value", "statistical_error", "systematic_components")},
                "source_locator": copy.deepcopy(row["locator"]), "source_remap": copy.deepcopy(remap),
                "original_depolarization_removed": view["depolarization_removed"],
                "source_blockers": list(view["blockers"])})
    if len(points) != selection["total_points"]:
        raise ValueError("incomplete requested historical points")
    return points


def prepare_request(manifest_path, cache, candidates, output):
    """Whole source-preparation job; all numerical prediction remains in Julia."""
    manifest_path, cache, candidates = Path(manifest_path).resolve(), Path(cache).resolve(), Path(candidates).resolve()
    root = manifest_path.parent.parent
    fingerprint = source_fingerprint(root)
    verified = verify_historical_intake(manifest_path, cache, candidates)
    candidate_path, report_path = candidates/"candidates.json", candidates/"report.json"
    candidate, report = read_json(candidate_path), read_json(report_path)
    contract_path, initial_path = root/CONTRACT, root/"config/initial-contract.json"
    contract, initial = read_json(contract_path), read_json(initial_path)
    points = request_points(candidate, contract, initial)
    if report["primary_HallA_points_crosschecked"] != len(points):
        raise ValueError("incomplete independently checked HallA source rows")
    inputs = {**report["source_hashes"], **{str(path): file_hash(path) for path in
        (candidate_path, report_path, contract_path, initial_path)}}
    request = {"schema_version": 1, "scope": SCOPE, "theory_scope": "reconstructed_theory",
        "observation_scope": CANDIDATE_SCOPE, **{key: False for key in FALSE_FLAGS},
        "contract_id": contract["id"], "contract_sha256": file_hash(contract_path),
        "initial_contract_sha256": file_hash(initial_path), "source_fingerprint": fingerprint,
        "parameters": {"order": initial["parameter_order"], "values": initial["published_parameters"],
                       "role": initial["published_parameter_role"]},
        "model_convention": {key: contract["prediction"][key] for key in
            ("phase", "normalization", "model_y", "source_mapping_status")},
        "points": points, "source_inputs": inputs,
        "source_preparation": {"status": "CURRENT_HISTORICAL_SOURCE_REPLAY_VERIFIED",
            "verification_kind": "Python_source_preparation_attestation_not_Julia_raw_table_replay",
            "scope": verified["scope"], "historical_counts": verified["counts"],
            "manifest": str(manifest_path), "cache": str(cache), "candidates": str(candidates),
            "candidate_sha256": file_hash(candidate_path), "report_sha256": file_hash(report_path)},
        "limitations": contract["limitations"]}
    if source_fingerprint(root) != fingerprint or any(file_hash(p) != value for p, value in inputs.items()):
        raise ValueError("historical point sources changed during request preparation")
    output = Path(output)
    atomic_json(output, request)
    return {"status": "PREDICTION_REQUEST_PREPARED", "scope": SCOPE, "points": len(points),
            "production_ready": False, "physical_extraction": False, "historical_fit_ready": False,
            "experimental_data_admitted": False, "request_sha256": file_hash(output)}


def verify_request(path):
    """Rebuild the request from its original locked inputs, including cached use."""
    path = Path(path)
    request = read_json(path)
    source = request["source_preparation"]
    # Prepare to a sibling temporary request, leaving the original bytes intact.
    with tempfile.TemporaryDirectory(prefix="historical-request-", dir=path.parent) as folder:
        fresh = Path(folder)/"request.json"
        prepare_request(source["manifest"], source["cache"], source["candidates"], fresh)
        if read_json(fresh) != request:
            raise ValueError("changed source-derived historical point request")
    return request


def render_comparison(report_path, output):
    """Render completed Julia values; no Python theory, fitting or error combination."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    report = read_json(report_path)
    if report["scope"] != SCOPE or report["status"] != "PASS" or any(report.get(k) is not False for k in FALSE_FLAGS):
        raise ValueError("only the completed scoped diagnostic can produce comparison artifacts")
    output = Path(output)
    table = output/"point-comparison.csv"
    with table.open("w", encoding="utf-8", newline="") as stream:
        stream.write("# scope="+SCOPE+"; fixed rounded parameters; no fit or likelihood\n")
        writer = csv.writer(stream)
        writer.writerow(["point_id", "hadron", "x", "y", "z", "Q2_GeV2", "pT_GeV", "measured",
                         "statistical_error", "experimental_systematic", "model_systematic", "fixed_theory",
                         "unweighted_observed_minus_predicted"])
        for row in report["points"]:
            k, m = row["kinematics"], row["measured"]
            writer.writerow([row["id"], row["hadron"], *[k[name] for name in ("x", "y", "z", "Q2_GeV2", "pT_GeV")],
                m["central_value"], m["statistical_error"], m["systematic_components"]["experimental"],
                m["systematic_components"]["model"], row["amplitude"], row["unweighted_observed_minus_predicted"]])
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), sharex=True)
    for ax, hadron, label in zip(axes, ("pi+", "pi-"), (r"$\pi^+$", r"$\pi^-$")):
        rows = [row for row in report["points"] if row["hadron"] == hadron]
        xs = [row["kinematics"]["x"] for row in rows]
        ax.errorbar(xs, [row["measured"]["central_value"] for row in rows],
            yerr=[row["measured"]["statistical_error"] for row in rows],
            fmt="o", color="#1b5a85", capsize=3, label="Hall A neutron: reported statistical error")
        ax.scatter(xs, [row["amplitude"] for row in rows], marker="D", facecolors="none", edgecolors="#c25b12",
                   s=45, zorder=4, label="Fixed rounded KPSY parameters")
        ax.axhline(0, color="0.55", linewidth=.7)
        ax.set_title(label)
        ax.set_xlabel(r"Reported $\langle x_B\rangle$")
        ax.set_ylabel(r"Collins amplitude $A_C$")
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=8, loc="best")
    fig.suptitle("Hall A fixed-mean comparison\nReconstructed historical diagnostic", fontsize=12)
    fig.text(.5, .035, "Conditional massless fixed-y model; source acceptance/spin corrections remain unresolved.\n"
        "Statistical bars shown; experimental/model errors remain separate in the table. No fit or likelihood.",
        ha="center", va="bottom", fontsize=8)
    fig.tight_layout(rect=(0, .14, 1, .88))
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
    output.mkdir(parents=True, exist_ok=True)
    report_path, audit_path, invocation_path = [output/name for name in ("predictions.json", "audit.json", "invocation.json")]
    if any(path.exists() for path in (report_path, audit_path, invocation_path)):
        raise ValueError("use a fresh historical prediction output directory")
    fingerprint = source_fingerprint(root)
    hashes = {str(p): file_hash(p) for p in (request_path, native_build, parton_build)}
    invocation = {"status": "FAIL", "scope": SCOPE, **{key: False for key in FALSE_FLAGS},
        "G2PERF_pass": False, "source_fingerprint": fingerprint, "code_commit": git_commit(root),
        "input_hashes": hashes, "steps": [], "outputs": {}}
    commands = [("predictions", _julia_command(julia, root, "predict_historical_points.jl",
                   request_path, native_build, parton_build, report_path)),
                ("audit", _julia_command(julia, root, "audit_historical_points.jl", request_path, report_path, audit_path))]
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
        command = _julia_command(julia, root, "audit_historical_points.jl", request_path, output/"predictions.json", audit_path)
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
