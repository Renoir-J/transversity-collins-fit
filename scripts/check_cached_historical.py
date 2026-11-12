"""Check cached candidate integrity without authorizing a fit."""
import sys
from pathlib import Path

from collinsflow.historical import verify_historical_intake
from collinsflow.provenance import atomic_json, read_json

manifest, cache, output = sys.argv[1:]
verify_historical_intake(manifest, cache, output)
path = Path(output) / "candidates.json"
original = path.read_bytes()
rejected = False
try:
    candidate = read_json(path)
    candidate["views"][0]["points"][0]["central_value"] += .001
    atomic_json(path, candidate)
    try:
        verify_historical_intake(manifest, cache, output)
    except ValueError:
        rejected = True
finally:
    path.write_bytes(original)
if not rejected:
    raise ValueError("cached historical intake accepted changed candidates")
verify_historical_intake(manifest, cache, output)
print("Historical candidate cache rejected tampering and accepted restored artifacts")
