"""Run or verify complete Julia fixed-context storage cases."""
from __future__ import annotations
import argparse
import json
from collinsflow.fixed_context_storage import run_case, verify_case, build_or_verify_ledger


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "verify", "ledger"):
        command = commands.add_parser(name)
        for key in ("sidis-request", "sia-request", "native-build", "parton-build", "run-dir"):
            command.add_argument("--" + key, required=True)
        command.add_argument("--julia", default="julia")
        command.add_argument("--timeout-seconds", type=float)
        if name == "ledger":
            command.add_argument("--output", required=True)
        else:
            command.add_argument("--case-index", type=int, required=True)
    args = parser.parse_args()
    common = (args.sidis_request, args.sia_request, args.native_build, args.parton_build, args.run_dir)
    kwargs = {"julia": args.julia, "timeout_seconds": args.timeout_seconds}
    if args.command == "run":
        result = run_case(*common, args.case_index, **kwargs)
    elif args.command == "verify":
        result = verify_case(*common, args.case_index, **kwargs)
    else:
        result = build_or_verify_ledger(*common, args.output, **kwargs)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
