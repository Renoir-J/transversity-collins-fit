"""Command-line entry point for complete, reproducible jobs."""

import argparse
import json
import sys

from .data import normalize_table, plan_update, validate_dataset
from .experimental import run_intake
from .pipeline import Blocked, execute, production_readiness, root_path, run_stage
from .provenance import read_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("readiness")
    intake = commands.add_parser("intake")
    intake.add_argument("--manifest", default="config/experimental-sources.json")
    intake.add_argument("--cache", default="data/raw/intake")
    intake.add_argument("--output", default="data/normalized/intake")
    intake.add_argument("--seed-root")
    intake.add_argument("--offline", action="store_true")
    validate = commands.add_parser("validate-data")
    validate.add_argument("data")
    normalize = commands.add_parser("normalize")
    normalize.add_argument("source")
    normalize.add_argument("mapping")
    normalize.add_argument("output")
    update = commands.add_parser("plan-update")
    update.add_argument("current")
    update.add_argument("incoming")
    update.add_argument("--sample-kind", choices=["replicas", "posterior"], default="replicas")
    for command in ("run", "stage"):
        sub = commands.add_parser(command)
        sub.add_argument("--mode", choices=["closure", "production"], default="closure")
        sub.add_argument("--run-dir", required=True)
        sub.add_argument("--julia", default="julia")
        sub.add_argument("--target", default="products") if command == "run" else sub.add_argument("stage")
    args = parser.parse_args()
    try:
        if args.command == "readiness":
            missing = production_readiness(root_path())
            print(json.dumps({"production_ready": not missing, "missing": missing}, indent=2))
            return 2 if missing else 0
        if args.command == "intake":
            result = run_intake(args.manifest, args.cache, args.output,
                                seed_root=args.seed_root, offline=args.offline)
        elif args.command == "validate-data":
            result = validate_dataset(read_json(args.data))
        elif args.command == "normalize":
            result = normalize_table(args.source, args.mapping, args.output)
        elif args.command == "plan-update":
            result = plan_update(read_json(args.current), read_json(args.incoming), args.sample_kind)
        elif args.command == "stage":
            result = run_stage(args.stage, args.run_dir, args.mode, args.julia)
            result = {"stage": result["stage"], "status": result["status"], "scope": result["scope"]}
        else:
            execute(args.run_dir, args.mode, args.julia, args.target)
            result = {"completed": args.target, "mode": args.mode, "run_dir": args.run_dir}
        print(json.dumps(result, indent=2))
    except Blocked as exc:
        print(str(exc), file=sys.stderr)
        return 2
    except Exception as exc:
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

