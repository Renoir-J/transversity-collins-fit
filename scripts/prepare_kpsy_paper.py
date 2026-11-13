"""Acquire and verify the exact primary paper used by the partial-moment check."""

from __future__ import annotations

import argparse
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "config/kpsy-paper-moment-contract.json"


class MathAnchors(HTMLParser):
    def __init__(self):
        super().__init__()
        self.values = {}

    def handle_starttag(self, tag, attrs):
        values = dict(attrs)
        if tag == "math" and "id" in values:
            key = values["id"]
            if key in self.values:
                raise ValueError(f"duplicate source math anchor: {key}")
            self.values[key] = values.get("alttext", "")


def verify_source(payload, contract):
    expected = contract["primary_source"]
    if len(payload) != expected["bytes"] or hashlib.sha256(payload).hexdigest() != expected["sha256"]:
        raise ValueError("primary paper byte hash/size mismatch; no replacement source is admitted")
    parser = MathAnchors()
    parser.feed(payload.decode("utf-8"))
    for key, value in expected["math_anchors"].items():
        if parser.values.get(key) != value:
            raise ValueError(f"primary paper math anchor mismatch: {key}")
    return {"status": "SOURCE_VERIFIED_NOT_A_PHYSICS_PASS", "scope": "source_acquisition",
            "production_ready": False, "version": expected["version"],
            "url": expected["url"], "sha256": expected["sha256"], "bytes": len(payload),
            "checked_anchors": len(expected["math_anchors"])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / ".runtime/paper-moments/source/1505.05589v1.html")
    args = parser.parse_args()
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if args.output.exists():
        payload = args.output.read_bytes()
    else:
        request = Request(contract["primary_source"]["url"], headers={"User-Agent": "CollinsFlow-source-review/1.0"})
        with urlopen(request, timeout=60) as response:
            payload = response.read(contract["primary_source"]["bytes"] + 1)
    receipt = verify_source(payload, contract)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if not args.output.exists():
        args.output.write_bytes(payload)
    args.output.with_suffix(".receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(receipt))


if __name__ == "__main__":
    main()
