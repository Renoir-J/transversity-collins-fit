"""Synthetic source-verifier controls; no physical coupling is generated in Python."""
from copy import deepcopy
import hashlib
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("prepare_native_coupling", ROOT / "scripts/prepare_native_coupling.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def synthetic_source():
    data = b'<html><math id="synthetic.alpha" alttext="synthetic alpha"></math></html>'
    return data, {"primary_source": {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(),
                                    "math_anchors": {"synthetic.alpha": "synthetic alpha"}}}


def test_acquisition_is_not_science_pass():
    data, policy = synthetic_source()
    result = MODULE.verify_paper(data, policy)
    assert result["status"] == "SOURCE_VERIFIED_NOT_A_SCIENCE_PASS"
    assert result["production_ready"] is False


def test_altered_bytes_and_mathematical_anchor_reject():
    data, policy = synthetic_source()
    with pytest.raises(ValueError, match="hash/size mismatch"):
        MODULE.verify_paper(data + b" ", policy)
    wrong = deepcopy(policy)
    wrong["primary_source"]["math_anchors"]["synthetic.alpha"] = "different alpha"
    with pytest.raises(ValueError, match="math anchor mismatch"):
        MODULE.verify_paper(data, wrong)


def test_duplicate_anchor_rejects_even_when_bytes_match():
    data, policy = synthetic_source()
    data += b'<math id="synthetic.alpha" alttext="other"></math>'
    policy["primary_source"].update(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
    with pytest.raises(ValueError, match="duplicate"):
        MODULE.verify_paper(data, policy)
