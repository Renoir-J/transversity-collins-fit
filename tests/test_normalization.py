import csv
from pathlib import Path

import pytest

from collinsflow.data import normalize_table, validate_dataset
from collinsflow.provenance import atomic_json, file_hash, read_json


def test_percent_sign_and_stat_covariance_conversion(tmp_path):
    raw = tmp_path / "measurement.csv"
    with raw.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerows([["A", "stat", "xlow", "xhigh"], [2.0, 0.5, 0.1, 0.2]])
    mapping = tmp_path / "mapping.json"
    atomic_json(mapping, {"source_sha256": file_hash(raw), "input_unit": "percent", "sign": -1,
        "dataset": {"id": "normalization_test", "publication": "synthetic normalization test",
            "table": "test table", "source_url": "https://example.invalid/test", "version": "test-v1",
            "source_sha256": file_hash(raw), "covariance_status": "published_uncorrelated",
            "selection_reason": "unit test only"}, "process": "SIDIS", "hadron": "pi+",
        "convention": "explicit-test-sign", "columns": {"value": "A", "stat": "stat"},
        "bins": {"x": ["xlow", "xhigh"]}})
    output = tmp_path / "normalized.json"
    normalize_table(raw, mapping, output)
    data = read_json(output)
    assert data["points"][0]["value"] == pytest.approx(-0.02)
    assert data["covariance"][0][0] == pytest.approx(0.005**2)
    raw.write_text(raw.read_text() + "3,1,0.2,0.3\n")
    with pytest.raises(ValueError, match="checksum"):
        normalize_table(raw, mapping, output)


def test_real_lattice_constraint_is_single_nonduplicated_isovector():
    root = Path(__file__).resolve().parents[1]
    data = read_json(root / "data/constraints/lattice-alexandrou2020.json")
    assert validate_dataset(data)["points"] == 1
    assert data["points"][0]["value"] == 0.936
    assert data["covariance"][0][0]**0.5 == 0.025

