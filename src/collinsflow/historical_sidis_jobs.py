"""Schedule complete Julia SIDIS jobs; Python performs no numerical predictions."""
from __future__ import annotations

import csv
import io
import json
from pathlib import Path
import subprocess
import time
import uuid

from .historical_sidis_pilot import _pilot_lock, _timed_job, _require_command_success, parse_resource_file
from .historical_sidis_points import verify_request
from .provenance import atomic_json, file_hash, git_commit, read_json, source_fingerprint, utc_now

ROOT = Path(__file__).resolve().parents[2]
POLICY = ROOT/"config/historical-sidis-job-policy.json"
POLICY_SHA = "7ce5d19f131c573caceaf64d24b2cc09fe589aaa4b628137d73d2b56df995911"
SCOPE = "historical_reconstructed_sidis_complete_point_jobs"
SCIENTIFIC_SCOPE = "historical_reconstructed_sidis_point_prediction_diagnostic"
FLAGS = ("production_ready", "physical_extraction", "historical_fit_ready", "experimental_data_admitted",
         "fit_performed", "accepted_prediction_diagnostic", "accepted_fitting_gate", "G2PERF_pass",
         "interval_coverage_verified", "all132_coverage_verified")
POINT_PHASES = ("native123", "native4", "tailref123", "reference4")
GLOBAL_PHASES = ("native-select3", "native-selection", "derivative-select3", "source-math", "diagnostic")
COMPLETED = "COMPLETE_JOB_RECORDED_NOT_PRODUCTION_GATE"
AUDITED = "SAVED_JOB_ARITHMETIC_REVIEWED_NOT_PRODUCTION_GATE"
TIME_BINARY = Path("/usr/bin/time")
FIGURES = ("comparison-hermes.png", "comparison-compass-proton.png", "comparison-compass-deuteron.png")


def load_policy():
    if file_hash(POLICY) != POLICY_SHA:
        raise ValueError("changed frozen complete SIDIS job policy")
    policy = read_json(POLICY)
    if (policy["scope"] != SCOPE or any(policy.get(k) is not False for k in FLAGS) or
            file_hash(ROOT/policy["parent_contract"]) != policy["parent_contract_sha256"]):
        raise ValueError("changed SIDIS scientific contract or self-admitting policy")
    return policy


def point_key(index):
    if isinstance(index, bool) or not isinstance(index, int) or not 1 <= index <= 132:
        raise ValueError("exact SIDIS source index1..132 required")
    return f"p{index:04d}"


def folder_for(run, phase, index=0):
    run = Path(run).resolve()
    if not run.is_relative_to(ROOT/"runs") or phase not in (*POINT_PHASES, *GLOBAL_PHASES):
        raise ValueError("use a canonical ignored full SIDIS run directory and declared phase")
    if phase in POINT_PHASES:
        return run/phase/point_key(index)
    if isinstance(index, bool) or not isinstance(index, int) or index != 0:
        raise ValueError("global jobs have point_index0")
    return run/phase


def _flags(item, final=False):
    return all(item.get(k) is (final and k in ("accepted_prediction_diagnostic", "all132_coverage_verified")) for k in FLAGS)


def _selection(run, phase, request, fingerprint):
    folder = folder_for(run, phase)
    report = read_json(folder/"report.json")
    audit = read_json(folder/"audit.json")
    if (report.get("phase") != phase or report.get("scope") != SCOPE or
            report.get("source_fingerprint") != fingerprint or not _flags(report) or
            report.get("source_ids") != [p["id"] for p in request["points"]] or
            audit.get("status") != AUDITED or audit.get("report_sha256") != file_hash(folder/"report.json")):
        raise ValueError("missing, stale or incomplete common SIDIS decision")
    allowed = {"native-select3": ("COMMON_NATIVE_LEVEL_SELECTED_NOT_DIAGNOSTIC_GATE", "NATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED"),
               "native-selection": ("COMMON_NATIVE_LEVEL_SELECTED_NOT_DIAGNOSTIC_GATE",),
               "derivative-select3": ("COMMON_DERIVATIVE_LEVEL_SELECTED_NOT_DIAGNOSTIC_GATE", "DERIVATIVE_LEVEL4_REQUIRED_NOT_ACCEPTED")}
    if report.get("status") not in allowed[phase]:
        raise ValueError("unaccepted or failed common SIDIS decision")
    return report


def input_paths(request_path, native, parton, run, phase, index, request, fingerprint):
    paths = [Path(request_path), Path(native), Path(parton), POLICY, TIME_BINARY, ROOT/"scripts/run_historical_sidis_job.py"]
    def add_points(name):
        paths.extend(folder_for(run, name, i)/"invocation.json" for i in range(1, 133))
    def add_global(name):
        paths.append(folder_for(run, name)/"invocation.json")
    if phase in POINT_PHASES:
        if index > 1:
            paths.append(folder_for(run, phase, index-1)/"invocation.json")
        if phase == "native4":
            paths.extend((folder_for(run, "native-select3")/"report.json", folder_for(run, "native123", index)/"report.json"))
        if phase in ("tailref123", "reference4"):
            chosen = _selection(run, "native-selection", request, fingerprint)
            paths.extend((folder_for(run, "native-selection")/"report.json",
                          folder_for(run, "native123" if chosen["selected_level"] == 3 else "native4", index)/"report.json"))
        if phase == "reference4":
            paths.extend((folder_for(run, "derivative-select3")/"report.json", folder_for(run, "tailref123", index)/"report.json"))
    elif phase != "source-math":
        add_points("native123")
        if phase != "native-select3":
            add_global("native-select3")
            if _selection(run, "native-select3", request, fingerprint)["selected_level"] == 0:
                add_points("native4")
        if phase in ("derivative-select3", "diagnostic"):
            add_global("native-selection")
            add_points("tailref123")
        if phase == "diagnostic":
            add_global("derivative-select3")
            add_global("source-math")
            if _selection(run, "derivative-select3", request, fingerprint)["selected_reference_level"] == 0:
                add_points("reference4")
    return sorted(set(p.resolve() for p in paths))


def _report_audit(folder, phase, index, request, fingerprint, request_sha):
    report, audit = read_json(folder/"report.json"), read_json(folder/"audit.json")
    final = phase == "diagnostic"
    if (report.get("source_fingerprint") != fingerprint or report.get("phase") != phase or
            report.get("scope") != (SCIENTIFIC_SCOPE if final else SCOPE) or not _flags(report, final) or
            report.get("source_points") != 132 or report.get("parameters") != request["parameters"] or
            report.get("request_sha256") != request_sha or report.get("policy_sha256") != POLICY_SHA or
            audit.get("status") != AUDITED or audit.get("source_fingerprint") != fingerprint or
            audit.get("scope") != report["scope"] or audit.get("request_sha256") != request_sha or audit.get("policy_sha256") != POLICY_SHA or
            audit.get("report_sha256") != file_hash(folder/"report.json") or not _flags(audit, final)):
        raise ValueError("incomplete, stale or self-admitting SIDIS report/audit")
    if phase in POINT_PHASES:
        if (report.get("status") != "POINT_JOB_RECORDED_NOT_DIAGNOSTIC_GATE" or report.get("point_index") != index or
                report.get("point") != request["points"][index-1]):
            raise ValueError("wrong exact source point or unfinished point job")
    elif phase == "diagnostic":
        if (report.get("status") != "CONDITIONAL_FULL132_POINT_DIAGNOSTIC_RECORDED_NOT_FITTING" or
                len(report.get("points", [])) != 132 or report.get("source_ids") != [p["id"] for p in request["points"]] or
                [p.get("id") for p in report["points"]] != report["source_ids"] or report.get("coordinates") != 13 or
                report.get("full132_native_tail_reference_derivative_policy_completed") is not True or
                type(report.get("native_level")) is not int or report["native_level"] not in (3, 4) or
                type(report.get("reference_level")) is not int or report["reference_level"] not in (3, 4) or
                any(p.get("native_level") != report["native_level"] or p.get("reference_level") != report["reference_level"] for p in report["points"])):
            raise ValueError("incomplete full132 diagnostic")
    elif phase == "source-math":
        if report.get("status") != "SYNTHETIC_SOURCE_MATH_RECORDED_NOT_DIAGNOSTIC_GATE":
            raise ValueError("missing explicitly synthetic source-math controls")
    else:
        _selection(folder.parent, phase, request, fingerprint)
    return report, audit


def _products(folder, final=False):
    names = ["report.json", "audit.json", "producer.log", "producer.resources.txt", "audit.log", "audit.resources.txt"]
    if final:
        names.extend(("point-comparison.csv", *FIGURES))
    return {str(folder/name) for name in names}


def verify_metadata(request_path, native, parton, run, phase, index=0):
    request_path, native, parton, run = map(lambda p: Path(p).resolve(), (request_path, native, parton, run))
    request = verify_request(request_path)
    fingerprint = source_fingerprint(ROOT)
    if (Path(request["source_preparation"]["manifest"]).parent.parent != ROOT or
            request["source_fingerprint"] != fingerprint):
        raise ValueError("current trusted full SIDIS source required")
    policy = load_policy()
    folder = folder_for(run, phase, index)
    invocation = read_json(folder/"invocation.json")
    final = phase == "diagnostic"
    if (invocation.get("status") != COMPLETED or invocation.get("source_fingerprint") != fingerprint or
            invocation.get("phase") != phase or type(invocation.get("point_index")) is not int or invocation.get("point_index") != index or not _flags(invocation, final) or
            invocation.get("policy_sha256") != POLICY_SHA or invocation.get("time_binary") != str(TIME_BINARY)):
        raise ValueError("stale, missing or failed complete SIDIS job")
    expected = {str(p) for p in input_paths(request_path, native, parton, run, phase, index, request, fingerprint)}
    if (set(invocation["input_hashes"]) != expected or
            any(file_hash(p) != h for p, h in invocation["input_hashes"].items()) or
            set(invocation["outputs"]) != _products(folder, final) or
            any(file_hash(p) != h for p, h in invocation["outputs"].items())):
        raise ValueError("changed or incomplete complete-job input/output bindings")
    if (file_hash(folder/"initial-invocation.json") != invocation["initial_record_sha256"] or
            read_json(folder/"initial-invocation.json")["record_role"] != "initial_sentinel_not_an_executed_failure"):
        raise ValueError("changed initial SIDIS invocation sentinel")
    if [s["name"] for s in invocation["steps"]] != ["producer", "audit"]:
        raise ValueError("incomplete SIDIS command ledger")
    budget = {"resources": policy["resources"][phase if phase in POINT_PHASES else "global"]}
    for step in invocation["steps"]:
        _require_command_success(step, budget)
        if (file_hash(folder/(step["name"]+".log")) != step["log_sha256"] or
                file_hash(folder/(step["name"]+".resources.txt")) != step["resource_sha256"] or
                parse_resource_file(folder/(step["name"]+".resources.txt")) != step["resource"]):
            raise ValueError("changed original command/resource records")
    if sum(s["wrapper_elapsed_seconds"] for s in invocation["steps"]) > budget["resources"]["whole_Julia_jobs_time_budget_seconds"]:
        raise ValueError("exhausted whole Julia-job budget")
    report, audit = _report_audit(folder, phase, index, request, fingerprint, file_hash(request_path))
    return request, policy, folder, invocation, report, audit


def _command(julia, program, *arguments):
    return [julia, "--startup-file=no", "--project="+str(ROOT/"numerics"), str(ROOT/"numerics/bin"/program), *map(str, arguments)]


def _deadline(limit, selected):
    if selected is None:
        return time.monotonic()+limit
    if isinstance(selected, bool) or not isinstance(selected, (int, float)) or not 0 < selected <= limit:
        raise ValueError("caller time budget may only tighten the frozen engineering limit")
    return time.monotonic()+selected


def _previous_guard(request_path, native, parton, run, phase, index):
    if phase in POINT_PHASES and index > 1:
        verify_metadata(request_path, native, parton, run, phase, index-1)


def run_complete_job(request_path, native, parton, run, phase, index=0, *, julia="julia", timeout_seconds=None):
    request_path, native, parton, run = map(lambda p: Path(p).resolve(), (request_path, native, parton, run))
    folder = folder_for(run, phase, index)
    if folder.exists() and any(folder.iterdir()):
        raise ValueError("preserve earlier complete/failed jobs; use a fresh run")
    folder.mkdir(parents=True, exist_ok=True)
    invocation = {"status": "FAIL", "scope": SCOPE, "phase": phase, "point_index": index,
                  "started_at": utc_now(), "steps": [], "outputs": {}, **dict.fromkeys(FLAGS, False)}
    initial = folder/"initial-invocation.json"
    atomic_json(initial, {**invocation, "record_role": "initial_sentinel_not_an_executed_failure"})
    invocation["initial_record_sha256"] = file_hash(initial)
    atomic_json(folder/"invocation.json", invocation)
    try:
        request = verify_request(request_path)
        fingerprint = source_fingerprint(ROOT)
        if Path(request["source_preparation"]["manifest"]).parent.parent != ROOT or request["source_fingerprint"] != fingerprint:
            raise ValueError("current trusted full132 request required")
        policy = load_policy()
        _previous_guard(request_path, native, parton, run, phase, index)
        inputs = input_paths(request_path, native, parton, run, phase, index, request, fingerprint)
        invocation.update(source_fingerprint=fingerprint, code_commit=git_commit(ROOT), policy_sha256=POLICY_SHA,
            input_hashes={str(p): file_hash(p) for p in inputs}, time_binary=str(TIME_BINARY),
            time_version=subprocess.run([str(TIME_BINARY), "--version"], capture_output=True, text=True, check=True, timeout=10).stdout)
        budget = {"resources": policy["resources"][phase if phase in POINT_PHASES else "global"]}
        report, audit = folder/"report.json", folder/"audit.json"
        producer = (_command(julia, "run_historical_sidis_point_job.jl", request_path, phase, index, native, parton, run, report)
                    if phase in POINT_PHASES else
                    _command(julia, "run_historical_sidis_global_job.jl", request_path, phase, native, parton, run, report))
        auditor = _command(julia, "audit_historical_sidis_job.jl", request_path, phase, index, native, parton, run, report, audit)
        with _pilot_lock(ROOT):
            deadline = _deadline(budget["resources"]["whole_Julia_jobs_time_budget_seconds"], timeout_seconds)
            for name, command in (("producer", producer), ("audit", auditor)):
                row = _timed_job(name, command, ROOT, folder, TIME_BINARY, budget, deadline)
                invocation["steps"].append(row)
                atomic_json(folder/"invocation.json", invocation)
                _require_command_success(row, budget)
        _report_audit(folder, phase, index, request, fingerprint, file_hash(request_path))
        if phase == "diagnostic":
            render_comparison(read_json(report), folder)
        verify_request(request_path)
        if source_fingerprint(ROOT) != fingerprint or any(file_hash(p) != h for p, h in invocation["input_hashes"].items()):
            raise ValueError("source or inputs changed during SIDIS job")
        invocation["outputs"] = {p: file_hash(p) for p in _products(folder, phase == "diagnostic")}
        invocation["status"] = COMPLETED
        if phase == "diagnostic":
            invocation["scope"] = SCIENTIFIC_SCOPE
            invocation["accepted_prediction_diagnostic"] = invocation["all132_coverage_verified"] = True
    except BaseException as error:
        invocation["error"] = type(error).__name__+": "+str(error)
        raise
    finally:
        invocation["completed_at"] = utc_now()
        atomic_json(folder/"invocation.json", invocation)
    return {k: invocation[k] for k in ("status", "phase", "point_index", *FLAGS)}


def verify_complete_job(request_path, native, parton, run, phase="diagnostic", index=0, *, julia="julia", timeout_seconds=None):
    request, policy, folder, invocation, report, _ = verify_metadata(request_path, native, parton, run, phase, index)
    review = Path(run).resolve()/"cache-reviews"/(phase+"-"+str(index)+"-"+uuid.uuid4().hex)
    review.mkdir(parents=True, exist_ok=False)
    record = {"status": "FAIL", "scope": "current_cached_SIDIS_saved_arithmetic_review_not_numerical_reexecution",
              "phase": phase, "point_index": index, "original_invocation_sha256": file_hash(folder/"invocation.json"),
              "source_fingerprint": report["source_fingerprint"], "started_at": utc_now(), **dict.fromkeys(FLAGS, False)}
    try:
        budget = {"resources": policy["resources"]["global"]}
        command = _command(julia, "audit_historical_sidis_job.jl", Path(request_path).resolve(), phase, index,
                           Path(native).resolve(), Path(parton).resolve(), Path(run).resolve(), folder/"report.json", review/"audit.json")
        with _pilot_lock(ROOT):
            row = _timed_job("cache-audit", command, ROOT, review, TIME_BINARY, budget,
                             _deadline(budget["resources"]["whole_Julia_jobs_time_budget_seconds"], timeout_seconds))
        record["step"] = row
        _require_command_success(row, budget)
        audit = read_json(review/"audit.json")
        if (audit["status"] != AUDITED or audit["report_sha256"] != file_hash(folder/"report.json") or
                audit["source_fingerprint"] != report["source_fingerprint"] or not _flags(audit, phase == "diagnostic")):
            raise ValueError("failed or stale current independent cached audit")
        if phase == "diagnostic" and (folder/"point-comparison.csv").read_text(encoding="utf-8") != comparison_csv(report):
            raise ValueError("changed current comparison CSV fields")
        verify_metadata(request_path, native, parton, run, phase, index)
        record["status"] = "CURRENT_CACHED_SIDIS_ARITHMETIC_VERIFIED_NOT_FITTING_GATE"
    except BaseException as error:
        record["error"] = type(error).__name__+": "+str(error)
        raise
    finally:
        record["completed_at"] = utc_now()
        atomic_json(review/"review.json", record)
    return record


def comparison_csv(report):
    fields = ("id", "dataset", "publication", "target", "hadron", "projection", "x", "z", "Q2_GeV2", "pT_GeV", "y",
              "measured", "statistical_error", "systematic_components", "reported_error_u", "uncertainty_status",
              "UU", "FUT", "amplitude", "unweighted_observed_minus_predicted", "native_level", "reference_level", "source_formula")
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for p in report["points"]:
        row = {"id": p["id"], "dataset": p["dataset_id"], "publication": p["publication"], "target": p["target"],
               "hadron": p["hadron"], "projection": p["projection"], **p["kinematics"],
               "measured": p["measured"]["central_value"], "statistical_error": p["measured"]["statistical_error"],
               "systematic_components": json.dumps(p["measured"]["systematic_components"], sort_keys=True),
               "reported_error_u": p["source_row"].get("reported_error_u"), "uncertainty_status": p["source_row"].get("uncertainty_status", "source_components_kept_separate"),
               "source_formula": p["model_convention"]["formula"],
               **{k: p[k] for k in ("UU", "FUT", "amplitude", "unweighted_observed_minus_predicted", "native_level", "reference_level")}}
        writer.writerow(row)
    return stream.getvalue()


def render_comparison(report, folder):
    if (report.get("status") != "CONDITIONAL_FULL132_POINT_DIAGNOSTIC_RECORDED_NOT_FITTING" or
            report.get("scope") != SCIENTIFIC_SCOPE or not _flags(report, True) or len(report.get("points", [])) != 132):
        raise ValueError("comparison products require the complete conditional132-point diagnostic")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    folder = Path(folder)
    csv_path = folder/"point-comparison.csv"
    if csv_path.exists() or any((folder/name).exists() for name in FIGURES):
        raise ValueError("preserve earlier comparison products")
    temp = csv_path.with_suffix(".csv.tmp")
    temp.write_text(comparison_csv(report), encoding="utf-8", newline="")
    temp.replace(csv_path)
    experiments = (("1006.4221v1", "HERMES proton", FIGURES[0]), ("1408.4405", "COMPASS proton", FIGURES[1]),
                   ("0802.2160", "COMPASS deuteron (conditional free isoscalar)", FIGURES[2]))
    for pub, title, filename in experiments:
        points = [p for p in report["points"] if p["publication"] == pub]
        fig, axes = plt.subplots(2, 3, figsize=(12, 7), squeeze=False)
        for row, hadron in enumerate(("pi+", "pi-")):
            for col, (projection, coordinate, xlabel) in enumerate((("x", "x", "reported mean x"), ("z", "z", "reported mean z"), ("pt", "pT_GeV", "reported mean P_h-perp [GeV]"))):
                group = [p for p in points if p["hadron"] == hadron and p["projection"] == projection]
                ax = axes[row, col]
                xs = [p["kinematics"][coordinate] for p in group]
                errors = [p["measured"]["statistical_error"] if p["measured"]["statistical_error"] is not None else p["source_row"]["reported_error_u"] for p in group]
                label = "reported Sigma (composition unresolved)" if pub == "1408.4405" else "reported statistical bars"
                ax.errorbar(xs, [p["measured"]["central_value"] for p in group], yerr=errors, fmt="o", ms=4, label=label)
                ax.plot(xs, [p["amplitude"] for p in group], "s", ms=5, mfc="none", label="fixed rounded KPSY parameters")
                ax.set_title(hadron+" / "+projection)
                ax.set_xlabel(xlabel)
                ax.set_ylabel("source Collins amplitude [fraction]")
                ax.axhline(0, color="gray", lw=.5)
                ax.legend(fontsize=6)
        fig.suptitle(title+" — reconstructed fixed-mean diagnostic")
        caption = {"1006.4221v1": "Statistical bars; experimental systematics remain separate in CSV.",
                   "1408.4405": "Reported Sigma bars; statistical/systematic composition remains unresolved.",
                   "0802.2160": "Statistical bars; unquantified systematics and the conditional free p+n target remain limitations."}[pub]
        fig.text(.5, .012, caption+" Fixed means, source phase, covariance/overlap and acceptance remain conditional. No fit interval.", ha="center", fontsize=6.5)
        fig.tight_layout(rect=(0, .04, 1, .95))
        temporary = folder/(filename+".tmp")
        fig.savefig(temporary, format="png", dpi=150)
        plt.close(fig)
        temporary.replace(folder/filename)
