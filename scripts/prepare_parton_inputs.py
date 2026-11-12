"""Build pinned inputs and isolated reference programs; this does not authorize a fit."""
import argparse
from pathlib import Path

from collinsflow.parton import prepare


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", default=".runtime/parton-inputs")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--fc", default="gfortran")
    parser.add_argument("--fc-flag", action="append", default=[])
    parser.add_argument("--fortran-runtime")
    parser.add_argument("--jobs", type=int, default=2, choices=range(1, 9))
    args = parser.parse_args()
    print(prepare(root, root / args.destination, offline=args.offline,
                  fc=args.fc, fc_flags=args.fc_flag,
                  fortran_runtime=args.fortran_runtime, jobs=args.jobs))


if __name__ == "__main__":
    main()
