"""Check cached candidate integrity without authorizing a fit."""
import sys
from pathlib import Path

from collinsflow.historical import verify_historical_intake
from collinsflow.provenance import atomic_json, file_hash, read_json

manifest, cache, output = sys.argv[1:]
verify_historical_intake(manifest, cache, output)
path = Path(output) / "candidates.json"
original = path.read_bytes()
report_path = Path(output) / "report.json"
original_report = report_path.read_bytes()
for mutation in ("value", "rehashed_value", "rehashed_angular_definition", "primary_decision"):
    rejected = False
    try:
        candidate = read_json(path)
        report = read_json(report_path)
        if mutation == "primary_decision":
            report["primary_BaBar_numeric_cells_crosschecked"] -= 1
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
    if not rejected:
        raise ValueError("cached historical intake accepted changed "+mutation)
    verify_historical_intake(manifest, cache, output)
verify_historical_intake(manifest, cache, output)
print("Historical candidate cache rejected tampering and accepted restored artifacts")
print("Rehashed values, angular definitions and primary decision counts were also rejected")
