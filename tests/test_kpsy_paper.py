"""Synthetic source-verifier controls; no paper download or physics PASS fixture."""

from copy import deepcopy
import hashlib
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("prepare_kpsy_paper", ROOT / "scripts/prepare_kpsy_paper.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def synthetic_source():
    payload = b'<html><math id="test.eq" alttext="synthetic h1"></math></html>'
    contract = {"primary_source": {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest(),
                "math_anchors": {"test.eq": "synthetic h1"}, "version": "synthetic", "url": "https://example.invalid"}}
    return payload, contract


def test_source_receipt_does_not_claim_physics_pass():
    payload, contract = synthetic_source()
    result = MODULE.verify_source(payload, contract)
    assert result["status"] == "SOURCE_VERIFIED_NOT_A_PHYSICS_PASS"
    assert result["production_ready"] is False
    assert result["checked_anchors"] == 1


def test_different_bytes_or_anchors_fail():
    payload, contract = synthetic_source()
    with pytest.raises(ValueError, match="byte hash/size mismatch"):
        MODULE.verify_source(payload.replace(b"h1", b"xh"), contract)
    altered = deepcopy(contract)
    altered["primary_source"]["math_anchors"]["test.eq"] = "wrong convention"
    with pytest.raises(ValueError, match="math anchor mismatch"):
        MODULE.verify_source(payload, altered)


def test_duplicate_source_anchor_fails_even_with_matching_bytes():
    payload, contract = synthetic_source()
    payload += b'<math id="test.eq" alttext="synthetic h1"></math>'
    contract["primary_source"].update(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
    with pytest.raises(ValueError, match="duplicate source math anchor"):
        MODULE.verify_source(payload, contract)
