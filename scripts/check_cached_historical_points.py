"""Exercise source/provenance and saved-numerical-operand cache controls."""
import argparse
import os
from pathlib import Path

from collinsflow.historical_points import verify_point_job
from collinsflow.provenance import atomic_json, file_hash, read_json

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("request")
parser.add_argument("output")
parser.add_argument("--julia", default="julia")
args = parser.parse_args()
request, output = Path(args.request).resolve(), Path(args.output).resolve()
verify_point_job(request, output, julia=args.julia)
paths = [request, output/"predictions.json", output/"audit.json", output/"invocation.json", output/"point-comparison.png"]
timestamps = {path: (path.stat().st_atime_ns, path.stat().st_mtime_ns) for path in paths}
original = {path: path.read_bytes() for path in paths}
for mutation in ("figure", "rehashed_request_geometry", "rehashed_source_mapping", "rehashed_prediction", "production_scope"):
    rejected = False
    try:
        invocation = read_json(output/"invocation.json")
        if mutation == "figure":
            (output/"point-comparison.png").write_bytes(b"changed fixture control; original is restored")
        elif mutation in ("rehashed_request_geometry", "rehashed_source_mapping"):
            data = read_json(request)
            if mutation == "rehashed_request_geometry":
                data["points"][0]["kinematics"]["x"] += .001
            else:
                data["points"][0]["original_depolarization_removed"] = None
            atomic_json(request, data)
            invocation["input_hashes"][str(request)] = file_hash(request)
        elif mutation == "rehashed_prediction":
            report = read_json(output/"predictions.json")
            report["points"][0]["amplitude"] += .001
            atomic_json(output/"predictions.json", report)
            audit = read_json(output/"audit.json")
            audit["report_sha256"] = file_hash(output/"predictions.json")
            atomic_json(output/"audit.json", audit)
            for name in ("predictions.json", "audit.json"):
                invocation["outputs"][str(output/name)] = file_hash(output/name)
        else:
            invocation["scope"] = "production"
        atomic_json(output/"invocation.json", invocation)
        try:
            verify_point_job(request, output, julia=args.julia)
        except ValueError:
            rejected = True
    finally:
        for path, data in original.items():
            path.write_bytes(data)
            os.utime(path, ns=timestamps[path])
    if not rejected:
        raise ValueError("cached historical point diagnostic accepted "+mutation)
    verify_point_job(request, output, julia=args.julia)
print("Historical point diagnostic cache rejected five controls and accepted restored source/numerical artifacts")
