"""Reject accidental secrets, private chat exports and generated bulk artifacts."""
import argparse
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DENIED = ("external/", ".runtime/", ".venv/", "runs/", "data/raw/", "artifacts/research-page.json")
SECRET = re.compile(rb"(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}|-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staged", action="store_true")
    args = parser.parse_args()
    command = ["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR", "-z"] if args.staged else ["git", "ls-files", "-z"]
    result = subprocess.run(command, cwd=ROOT, capture_output=True, check=True)
    failures = []
    for name in result.stdout.decode().split("\0"):
        if not name:
            continue
        path = ROOT / name
        if name.startswith(DENIED) or name.endswith((".pem", ".key")):
            failures.append(f"excluded material: {name}")
        if path.is_file():
            content = subprocess.run(["git", "show", ":" + name], cwd=ROOT, capture_output=True, check=True).stdout if args.staged else path.read_bytes()
            if len(content) > 2_000_000:
                failures.append(f"large artifact: {name}")
            if SECRET.search(content):
                failures.append(f"credential-like content: {name}")
            private_link = b"qcdphysics.slack.com" + b"/archives/"
            if private_link in content:
                failures.append(f"private chat export/link: {name}")
    if failures:
        raise SystemExit("\n".join(failures))
    print("Repository content checks passed")


if __name__ == "__main__":
    main()

