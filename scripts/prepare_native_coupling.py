"""Lock the primary paper and build an isolated pinned native coupling comparator."""
from __future__ import annotations

import argparse
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.request import Request, urlopen

from collinsflow.hoppet import prepare

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "config/native-coupling-contract.json"


class MathAnchors(HTMLParser):
    def __init__(self):
        super().__init__()
        self.values = {}

    def handle_starttag(self, tag, attrs):
        values = dict(attrs)
        if tag == "math" and "id" in values:
            if values["id"] in self.values:
                raise ValueError("duplicate source math anchor")
            self.values[values["id"]] = values.get("alttext", "")


def verify_paper(payload, policy):
    expected = policy["primary_source"]
    if len(payload) != expected["bytes"] or hashlib.sha256(payload).hexdigest() != expected["sha256"]:
        raise ValueError("primary paper hash/size mismatch")
    parser = MathAnchors()
    parser.feed(payload.decode("utf-8"))
    if any(parser.values.get(key) != value for key, value in expected["math_anchors"].items()):
        raise ValueError("primary paper math anchor mismatch")
    return {"status": "SOURCE_VERIFIED_NOT_A_SCIENCE_PASS", "production_ready": False,
            "sha256": expected["sha256"], "checked_anchors": len(expected["math_anchors"])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", default=".runtime/native-coupling")
    parser.add_argument("--source")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--fc", default="gfortran")
    parser.add_argument("--fc-flag", action="append", default=[])
    parser.add_argument("--jobs", type=int, choices=range(1, 9), default=2)
    args = parser.parse_args()
    destination = (ROOT / args.destination).resolve()
    if not destination.is_relative_to(ROOT / ".runtime"):
        raise ValueError("native coupling build must stay inside ignored .runtime")
    policy = json.loads(CONTRACT.read_text(encoding="utf-8"))
    paper = destination / "paper.html"
    if paper.is_file():
        payload = paper.read_bytes()
    else:
        if args.offline:
            raise FileNotFoundError("locked primary paper is unavailable offline")
        with urlopen(Request(policy["primary_source"]["url"], headers={"User-Agent": "CollinsFlow-source-review/1.0"}), timeout=60) as response:
            payload = response.read(policy["primary_source"]["bytes"] + 1)
    receipt = verify_paper(payload, policy)
    destination.mkdir(parents=True, exist_ok=True)
    if not paper.exists():
        paper.write_bytes(payload)
    paper.with_suffix(".receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(prepare(ROOT, destination, source=args.source, offline=args.offline, cmake=args.cmake,
                  fc=args.fc, fc_flags=args.fc_flag, jobs=args.jobs, with_coupling_reference=True))


if __name__ == "__main__":
    main()
