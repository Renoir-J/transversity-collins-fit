"""Build the pinned HOPPET adapter; this receipt does not authorize production."""
import argparse
from pathlib import Path

from collinsflow.hoppet import prepare


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", default=".runtime/hoppet")
    parser.add_argument("--source")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--fc", default="gfortran")
    parser.add_argument("--fc-flag", action="append", default=[])
    parser.add_argument("--jobs", type=int, choices=range(1, 9), default=2)
    args = parser.parse_args()
    print(prepare(root, root / args.destination, source=args.source, offline=args.offline,
                  cmake=args.cmake, fc=args.fc, fc_flags=args.fc_flag, jobs=args.jobs))


if __name__ == "__main__":
    main()
