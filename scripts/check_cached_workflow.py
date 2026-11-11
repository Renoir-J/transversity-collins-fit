"""Integration check: cached Snakemake runs must reject a changed product."""
import argparse
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    target = Path(args.run_dir) / "products/per-bin-predictions.csv"
    original = target.read_bytes()
    command = ["snakemake", "--snakefile", "workflow/Snakefile", "--cores", "2", "--quiet", "--config",
               "mode=closure", f"run_dir={args.run_dir}"]
    try:
        target.write_bytes(original + b"TAMPER_TEST\n")
        rejected = subprocess.run(command, cwd=root, capture_output=True, text=True)
    finally:
        target.write_bytes(original)
    if rejected.returncode == 0 or "changed receipt output" not in rejected.stdout + rejected.stderr:
        raise RuntimeError("cached workflow did not reject the intentionally changed product")
    accepted = subprocess.run(command, cwd=root, capture_output=True, text=True)
    if accepted.returncode != 0:
        raise RuntimeError("restored cached workflow did not validate: " + accepted.stderr[-3000:])
    print("Cached workflow rejected tampering and accepted restored artifacts")


if __name__ == "__main__":
    main()

