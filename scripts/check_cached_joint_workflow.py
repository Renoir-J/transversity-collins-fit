"""Check that cached joint validation rejects a changed report and accepts restoration."""
import argparse
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir")
    parser.add_argument("--native-build", default=".runtime/hoppet/build.json")
    parser.add_argument("--parton-build", default=".runtime/parton-inputs/build.json")
    parser.add_argument("--kind", choices=("all", "analytic", "inputs"), default="all")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    run = Path(args.run_dir)
    kinds = ("analytic", "inputs") if args.kind == "all" else (args.kind,)
    command = ["snakemake", "--snakefile", "workflow/joint-likelihood.smk", "--cores", "1", "--quiet",
               "--config", f"run_dir={run}", f"kind={args.kind}", f"native_build={args.native_build}",
               f"parton_build={args.parton_build}"]
    for kind in kinds:
        target = run / kind / "validation.json"
        original = target.read_bytes()
        changed = json.loads(original)
        changed["profile"]["chi2"] += 0.25
        try:
            target.write_text(json.dumps(changed) + "\n", encoding="utf-8", newline="\n")
            rejected = subprocess.run(command, cwd=root, capture_output=True, text=True)
        finally:
            target.write_bytes(original)
        (run / kind / "cache-rejected.log").write_text(rejected.stdout + rejected.stderr, encoding="utf-8")
        if rejected.returncode == 0 or "numerical operand differs" not in rejected.stdout + rejected.stderr:
            raise RuntimeError("cached joint workflow did not reject the changed report")
    restored = subprocess.run(command, cwd=root, capture_output=True, text=True)
    (run / "cache-restored.log").write_text(restored.stdout + restored.stderr, encoding="utf-8")
    if restored.returncode != 0:
        raise RuntimeError("restored joint workflow did not validate: " + restored.stderr[-3000:])
    print("Cached joint workflow rejected changed reports and accepted restored artifacts")


if __name__ == "__main__":
    main()
