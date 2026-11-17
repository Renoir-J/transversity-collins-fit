"""Complete synthetic source-control jobs only; no physical SIA action."""
import argparse
import importlib.util
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("sia_source_control_jobs", root / "src/collinsflow/sia_source_controls.py")
jobs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(jobs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)
    run = actions.add_parser("run")
    run.add_argument("plan")
    run.add_argument("--expected-sha256", required=True)
    run.add_argument("--output", required=True)
    run.add_argument("--julia", default="julia")
    verify = actions.add_parser("verify")
    verify.add_argument("plan")
    verify.add_argument("--expected-sha256", required=True)
    verify.add_argument("output")
    verify.add_argument("--julia", default="julia")
    args = parser.parse_args()
    if args.action == "run":
        result = jobs.run_pair(args.plan, args.expected_sha256, args.output, julia=args.julia)
    else:
        result = jobs.verify_pair(args.output, args.plan, args.expected_sha256, julia=args.julia)
    print(json.dumps({key: result[key] for key in ("status", "scope", "fit_performed", "production_ready")}))


if __name__ == "__main__":
    main()
