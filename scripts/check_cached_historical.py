"""Check cached candidate integrity without authorizing a fit."""
import sys
import os
from pathlib import Path

from collinsflow.historical import verify_historical_intake
from collinsflow.provenance import atomic_json, file_hash, read_json

manifest, cache, output = sys.argv[1:]
verify_historical_intake(manifest, cache, output)
path = Path(output) / "candidates.json"
report_path = Path(output) / "report.json"
timestamps = {p: (p.stat().st_atime_ns, p.stat().st_mtime_ns) for p in (path, report_path)}
original = path.read_bytes()
original_report = report_path.read_bytes()
for mutation in ("value", "rehashed_value", "rehashed_angular_definition", "primary_decision",
                 "rehashed_proton_mean_provenance", "rehashed_proton_unbounded_bin", "proton_primary_decision"):
    rejected = False
    try:
        candidate = read_json(path)
        report = read_json(report_path)
        if mutation == "primary_decision":
            report["primary_BaBar_numeric_cells_crosschecked"] -= 1
        elif mutation == "proton_primary_decision":
            report["primary_COMPASS_proton_exact_numeric_cells"] -= 1
        elif mutation == "rehashed_proton_mean_provenance":
            view = next(v for v in candidate["views"] if
                        v.get("primary_numeric_crosscheck", {}).get("crosscheck_kind") == "COMPASS_proton_partial")
            point = view["points"][0]
            point["primary_reported_means"]["x"]["primary_literal"] = "0.125"
        elif mutation == "rehashed_proton_unbounded_bin":
            view = next(v for v in candidate["views"] if v["projection"] == "pt" and
                        v.get("primary_numeric_crosscheck", {}).get("crosscheck_kind") == "COMPASS_proton_partial")
            view["excluded"][-1]["primary_projection_bin"].update(high=2., upper_unbounded=False)
        elif mutation == "rehashed_angular_definition":
            view = next(v for v in candidate["views"] if v["process"] == "SIA" and v.get("primary_numeric_crosscheck"))
            view["points"][0]["angular_source"]["average_type"] = "ratio_of_angular_means"
        else:
            candidate["views"][0]["points"][0]["central_value"] += .001
        atomic_json(path, candidate)
        if mutation.startswith("rehashed_"):
            report["candidate_sha256"] = file_hash(path)
        atomic_json(report_path, report)
        try:
            verify_historical_intake(manifest, cache, output)
        except ValueError:
            rejected = True
    finally:
        path.write_bytes(original)
        report_path.write_bytes(original_report)
        for p in (path, report_path):
            os.utime(p, ns=timestamps[p])
    if not rejected:
        raise ValueError("cached historical intake accepted changed "+mutation)
    verify_historical_intake(manifest, cache, output)
verify_historical_intake(manifest, cache, output)
print("Historical candidate cache rejected tampering and accepted restored artifacts")
print("Rehashed values, angular/primary-mean provenance, excluded unbounded bins and primary counts were rejected")
