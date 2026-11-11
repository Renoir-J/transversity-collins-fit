from pathlib import Path

import pytest

from collinsflow.pipeline import production_readiness
from collinsflow.provenance import atomic_json, digest, file_hash, verify_receipt


def test_modified_or_wrong_scope_receipts_rejected(tmp_path):
    artifact = tmp_path / "artifact.json"
    atomic_json(artifact, {"value": 1})
    receipt = tmp_path / "receipt.json"
    payload = {"status": "PASS", "scope": "synthetic_closure", "source_fingerprint": "abc",
               "checks": [{"passed": True}], "input_hashes": {},
               "outputs": [{"path": str(artifact), "sha256": file_hash(artifact)}]}
    atomic_json(receipt, payload)
    verify_receipt(receipt, "synthetic_closure", "abc")
    with pytest.raises(ValueError, match="wrong-scope"):
        verify_receipt(receipt, "production", "abc")
    with pytest.raises(ValueError, match="stale"):
        verify_receipt(receipt, "synthetic_closure", "changed-code")
    atomic_json(artifact, {"value": 2})
    with pytest.raises(ValueError, match="changed receipt output"):
        verify_receipt(receipt, "synthetic_closure", "abc")


def test_cache_keys_change_with_every_dependency():
    dependencies = dict(theory="a", cuts="b", scales="c", spin="d", grid="e", parameters="f")
    original = digest(dependencies)
    for key in dependencies:
        changed = dict(dependencies)
        changed[key] += "changed"
        assert digest(changed) != original


def test_production_is_explicitly_blocked():
    root = Path(__file__).resolve().parents[1]
    assert any("production_ready" in item for item in production_readiness(root))

