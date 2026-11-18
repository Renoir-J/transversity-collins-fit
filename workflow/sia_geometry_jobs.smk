"""Serialize complete Julia pairs and always perform a fresh final cache review."""
from pathlib import Path

ROOT = Path(workflow.basedir).parent
MANIFEST = config["manifest"]
MANIFEST_SHA = config["manifest_sha256"]

rule all:
    retries: 0
    run:
        shell("python {ROOT}/scripts/run_sia_geometry_jobs.py all --manifest {MANIFEST:q} --manifest-sha {MANIFEST_SHA:q}")

rule cached_review:
    retries: 0
    run:
        shell("python {ROOT}/scripts/run_sia_geometry_jobs.py verify --manifest {MANIFEST:q} --manifest-sha {MANIFEST_SHA:q}")
