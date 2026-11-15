"""Check cached fitting rejects report tampering and accepts the original bytes."""
import argparse
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir")
    parser.add_argument("--algorithm", choices=("LN_SBPLX", "LN_BOBYQA"), required=True)
    parser.add_argument("--maxeval", type=int, choices=(4000, 20000, 80000), required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    run = Path(args.run_dir)
    report = run / "validation.json"
    original = report.read_bytes()
    command = ["snakemake", "--snakefile", "workflow/bounded-joint-fit.smk", "--cores", "1", "--quiet",
               "--config", f"run_dir={run}", f"algorithm={args.algorithm}", f"maxeval={args.maxeval}"]
    changed = json.loads(original)
    changed["fit"]["chi2"] += 0.25
    try:
        report.write_text(json.dumps(changed) + "\n", encoding="utf-8", newline="\n")
        rejected = subprocess.run(command, cwd=root, capture_output=True, text=True)
    finally:
        report.write_bytes(original)
    (run / "cache-rejected.log").write_text(rejected.stdout + rejected.stderr, encoding="utf-8")
    if rejected.returncode == 0 or "changed numerical operand" not in rejected.stdout + rejected.stderr:
        raise RuntimeError("cached bounded workflow did not reject changed report")
    other = "LN_SBPLX" if args.algorithm == "LN_BOBYQA" else "LN_BOBYQA"
    different_budget = 4000 if args.maxeval != 4000 else 20000
    for name, key, value, message in (
        ("wrong-algorithm", "algorithm", other, "stale requested bounded configuration"),
        ("wrong-budget", "maxeval", str(different_budget), "stale requested bounded configuration"),
        ("production-mode", "mode", "production", "synthetic closure only"),
        ("production-scope", "scope", "production", "synthetic closure only"),
    ):
        altered = [arg if not arg.startswith(key + "=") else f"{key}={value}" for arg in command]
        if not any(arg.startswith(key + "=") for arg in command):
            altered.append(f"{key}={value}")
        altered.extend(["--rerun-triggers", "mtime"])
        failure = subprocess.run(altered, cwd=root, capture_output=True, text=True)
        (run / ("cache-" + name + ".log")).write_text(failure.stdout + failure.stderr, encoding="utf-8")
        if failure.returncode == 0 or message not in failure.stdout + failure.stderr:
            raise RuntimeError("cached bounded workflow failed control: " + name)
    restored = subprocess.run(command, cwd=root, capture_output=True, text=True)
    (run / "cache-restored.log").write_text(restored.stdout + restored.stderr, encoding="utf-8")
    if restored.returncode:
        raise RuntimeError("restored bounded workflow did not validate: " + restored.stderr[-3000:])
    print("Cached bounded workflow rejected tampering and accepted restored artifacts")


if __name__ == "__main__":
    main()
