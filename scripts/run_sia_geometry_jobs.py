"""Prepare, execute or freshly review the complete conditional SIA diagnostic."""
import argparse
import importlib.util
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("sia_complete_jobs_adapter", root / "src/collinsflow/sia_geometry_execution.py")
jobs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(jobs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)
    prepare = actions.add_parser("prepare")
    for name in ("plan", "plan-sha", "native", "parton", "source-controls", "sidis-owner", "retained-inventory", "output"):
        prepare.add_argument("--" + name, required=True)
    prepare.add_argument("--sidis-run", action="append", required=True)
    prepare.add_argument("--julia", default="julia")
    for name in ("all", "verify", "plan"):
        action = actions.add_parser(name)
        action.add_argument("--manifest", required=True)
        action.add_argument("--manifest-sha", required=True)
    args = parser.parse_args()
    if args.action == "prepare":
        result = jobs.prepare(args.plan, args.plan_sha, args.native, args.parton, args.source_controls,
            args.sidis_owner, args.sidis_run, args.retained_inventory, args.output, julia=args.julia)
    elif args.action == "all":
        result = jobs.run_all(args.manifest, args.manifest_sha)
    elif args.action == "verify":
        result = jobs.cached_products(args.manifest, args.manifest_sha)
    else:
        m = jobs.load(args.manifest, args.manifest_sha)
        result = dict(status="COMPLETE_JOB_PLAN_ONLY_NO_PHYSICAL_DISPATCH", native123_jobs=36, conditional_native4_jobs=36,
            tailref_jobs=36, derivative_geometry_indices=[1, 15, 36], global_jobs=list(jobs.GLOBALS),
            host=jobs.host.observe_host(m["sidis_owner"], m["sidis_run_roots"]), fit_performed=False)
    print(json.dumps({key: result[key] for key in result if key in ("status", "scope", "fit_performed", "production_ready",
        "accepted_fitting_gate", "host", "native123_jobs", "conditional_native4_jobs", "tailref_jobs", "derivative_geometry_indices", "global_jobs")}))


if __name__ == "__main__":
    main()
