"""Source-only planner; this CLI has no numerical execution action."""
import argparse
import importlib.util
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("sia_job_planner", root / "src/collinsflow/sia_geometry_jobs.py")
planner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(planner)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)
    plan = actions.add_parser("plan")
    plan.add_argument("--frozen", required=True)
    plan.add_argument("--request", required=True)
    plan.add_argument("--output", required=True)
    verify = actions.add_parser("verify")
    verify.add_argument("plan")
    verify.add_argument("--expected-sha256", required=True)
    args = parser.parse_args()
    if args.action == "plan":
        result = planner.prepare_plan(args.frozen, args.request, args.output)
    else:
        result = planner.verify_plan(args.plan, expected_sha256=args.expected_sha256)
    print(json.dumps({key: result[key] for key in ("status", "scope", "numerical_dispatch_ready", "fit_performed")}))


if __name__ == "__main__":
    main()
