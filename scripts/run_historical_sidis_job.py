"""Run or independently recheck one canonical complete Julia SIDIS job."""
import argparse
import json

from collinsflow.historical_sidis_jobs import GLOBAL_PHASES, POINT_PHASES, run_complete_job, verify_complete_job


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "verify"))
    parser.add_argument("--request", required=True)
    parser.add_argument("--native-build", required=True)
    parser.add_argument("--parton-build", required=True)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--phase", choices=(*POINT_PHASES, *GLOBAL_PHASES), required=True)
    parser.add_argument("--point-index", type=int, default=0)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--timeout-seconds", type=float)
    args = parser.parse_args()
    function = run_complete_job if args.action == "run" else verify_complete_job
    result = function(args.request, args.native_build, args.parton_build, args.run_dir, args.phase,
                      args.point_index, julia=args.julia, timeout_seconds=args.timeout_seconds)
    print(json.dumps({k: result[k] for k in ("status", "phase", "point_index")}))


if __name__ == "__main__":
    main()
