"""Schedule immutable continuation of the supported frozen SIDIS attempt."""
from pathlib import Path
import argparse
import importlib.util
import json
import sys

ADAPTER = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("sidis_continuation_transport", ADAPTER/"src/collinsflow/sidis_continuation.py")
transport = importlib.util.module_from_spec(spec)
spec.loader.exec_module(transport)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "plan", "prefix", "run", "verify", "all"))
    p.add_argument("--frozen-checkout", required=True)
    p.add_argument("--manifest")
    for name in ("original", "review", "integrity", "loss", "output", "native-build", "parton-build"):
        p.add_argument("--"+name)
    p.add_argument("--phase", choices=(*transport.POINT_PHASES, *transport.GLOBAL_PHASES), default="diagnostic")
    p.add_argument("--point-index", type=int, default=0)
    p.add_argument("--julia", default="julia")
    p.add_argument("--timeout-seconds", type=float)
    a = p.parse_args()
    frozen = Path(a.frozen_checkout).resolve()
    transport.contract(ADAPTER)
    transport.verify_git(frozen, transport.FROZEN_HEAD)
    # This process has not imported the adapter's collinsflow package. Import
    # the explicit, Git-verified frozen package, without replacing its globals.
    sys.path.insert(0, str(frozen/"src"))
    from collinsflow import historical_sidis_jobs as api
    if api.ROOT != frozen:
        raise ValueError("wrong frozen collinsflow import")
    if a.action == "prepare":
        names = ("original", "review", "integrity", "loss", "output", "native_build", "parton_build")
        if any(getattr(a, name) is None for name in names):
            p.error("prepare requires original/review/integrity/loss/output/native-build/parton-build")
        result = transport.prepare(frozen, ADAPTER, a.original, a.review, a.integrity, a.loss,
                                   a.output, a.native_build, a.parton_build, api)
        print(json.dumps({"manifest": str(result), "manifest_sha256": transport.file_hash(result),
                          "status": "PREPARED_NOT_EXECUTED_NOT_A_GATE"}), flush=True)
        return
    if a.manifest is None:
        p.error("an explicit prepared manifest is required")
    path = Path(a.manifest).resolve()
    m = transport.verify_manifest(path, frozen, ADAPTER, api)
    if a.action == "plan":
        result = {"status": "PLAN_ONLY_NO_JULIA_DISPATCH", "original_completed_prefix": 13,
                  "remaining_native123_points": 119, "all132_native4_if_required": True,
                  "all132_all13_tail_and_reference_required": True, "original_files_unchanged": 141,
                  "production_ready": False, "fit_performed": False}
    elif a.action == "prefix":
        result = transport.audit_prefix(m, path, api, a.julia)
    elif a.action == "run":
        result = transport.run_job(m, path, a.phase, a.point_index, api, a.julia, a.timeout_seconds)
    elif a.action == "verify":
        result = transport.audit_cached(m, path, a.phase, a.point_index, api, a.julia)
    else:
        result = transport.run_all(m, path, api, a.julia)
    print(json.dumps({k: result[k] for k in ("status", *transport.FLAGS) if k in result}), flush=True)


if __name__ == "__main__":
    main()
