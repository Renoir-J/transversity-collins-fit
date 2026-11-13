"""Exercise experimental cache negatives without authorizing a physical fit."""
import sys
from pathlib import Path

from collinsflow.experimental import verify_intake
from collinsflow.provenance import atomic_json, file_hash, read_json


def main():
    manifest, cache, output = sys.argv[1:]
    output = Path(output)
    candidate_path, report_path = output / "candidates.json", output / "intake-report.json"
    receipt_path = output / "cache-integrity-check.json"
    atomic_json(receipt_path, {"scope": "candidate_cache_integrity", "status": "IN_PROGRESS",
                              "production_ready": False, "physical_extraction": False})
    original_candidate, original_report = candidate_path.read_bytes(), report_path.read_bytes()
    verify_intake(manifest, cache, output)
    checks = []
    mutations = ["value", "rehashed_value", "rehashed_harmonic", "rehashed_overlap", "rehashed_admission_type",
                 "report_count", "report_decision", "report_admission", "report_flag_type"]
    if "COMPASS2022_source_review" in read_json(report_path):
        mutations += ["rehashed_compass_total", "rehashed_compass_decimal", "rehashed_compass_sample",
                      "report_compass_flag"]
    for mutation in mutations:
        rejected = False
        try:
            candidate, report = read_json(candidate_path), read_json(report_path)
            if mutation == "report_count":
                report["points_across_overlapping_views"] += 1
            elif mutation == "report_decision":
                report["STAR_source_review"]["hepdata_views"] -= 1
            elif mutation == "report_admission":
                report["physical_extraction"] = True
            elif mutation == "report_flag_type":
                report["HERMES_3D_primary_binning"][0]["passed"] = 1
            elif mutation == "report_compass_flag":
                report["COMPASS2022_source_review"]["partial_2022_sample"] = 1
            elif mutation == "rehashed_admission_type":
                candidate["production_ready"] = 0
            elif mutation.startswith("rehashed_compass_"):
                dataset = next(v for v in candidate["views"]
                               if v.get("scientific_contract") == "COMPASS2022-DATA-001")
                item = dataset["points"][0]
                if mutation == "rehashed_compass_total":
                    item["reported_total_error"] *= 2
                elif mutation == "rehashed_compass_decimal":
                    item["source_decimal_fields"]["er_tot"] = str(item["reported_total_error"] * 2)
                else:
                    dataset["sample"] = "synthetic_wrong_complete_2022_sample"
            elif mutation in {"rehashed_harmonic", "rehashed_overlap"}:
                dataset = next(v for v in candidate["views"] if v["process"] == "STAR_510")
                if mutation == "rehashed_harmonic":
                    dataset["convention"] = "synthetic_wrong_sin(phiS+phiH)"
                else:
                    dataset["overlap_group"] = "synthetic_wrong_disjoint_sample"
            else:
                candidate["views"][0]["points"][0]["value"] += 0.001
            atomic_json(candidate_path, candidate)
            if mutation.startswith("rehashed_"):
                report["candidate_file_sha256"] = file_hash(candidate_path)
            atomic_json(report_path, report)
            try:
                verify_intake(manifest, cache, output)
            except ValueError as exc:
                expected = ("artifact was changed" if mutation == "value" else
                            "candidate differs from locked source replay" if mutation.startswith("rehashed_") else
                            "report differs from locked source replay")
                if expected not in str(exc):
                    raise
                rejected = True
                checks.append({"name": mutation, "rejected": True, "reason": str(exc)})
        finally:
            candidate_path.write_bytes(original_candidate)
            report_path.write_bytes(original_report)
        if not rejected:
            atomic_json(receipt_path, {"scope": "candidate_cache_integrity", "status": "FAIL",
                "production_ready": False, "physical_extraction": False, "accepted_mutation": mutation})
            raise ValueError("experimental cache accepted changed " + mutation)
    verification = verify_intake(manifest, cache, output)
    assert candidate_path.read_bytes() == original_candidate and report_path.read_bytes() == original_report
    report = read_json(report_path)
    atomic_json(receipt_path, {
        "scope": "candidate_cache_integrity", "status": "PASS", "production_ready": False,
        "physical_extraction": False, "source_fingerprint": report["source_fingerprint"],
        "manifest_sha256": file_hash(manifest), "candidate_sha256": file_hash(candidate_path),
        "intake_report_sha256": file_hash(report_path), "checker_sha256": file_hash(__file__),
        "checks": checks, "restored_byte_identical": True, "restored_verification": verification})
    print("Experimental candidate cache rejected tampering and accepted restored artifacts")
    print("Rehashed values, harmonics, overlap, scalar types and source decisions were also rejected")


if __name__ == "__main__":
    main()
