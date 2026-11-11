from pathlib import Path

import numpy as np
import pytest

from collinsflow.data import plan_update, validate_dataset
from collinsflow.provenance import read_json

ROOT = Path(__file__).resolve().parents[1]


def fixture():
    return read_json(ROOT / "fixtures/closure/data.json")


def test_fixture_covariance_and_scope():
    assert validate_dataset(fixture())["points"] == 66


@pytest.mark.parametrize("mutation", ["percent", "overlap", "negative_cov", "hadron", "replaces"])
def test_invalid_data_stops_before_fitting(mutation):
    data = fixture()
    if mutation == "percent":
        data["points"][0]["unit"] = "percent"
    elif mutation == "overlap":
        data["datasets"][1]["overlap_group"] = data["datasets"][0]["overlap_group"]
    elif mutation == "negative_cov":
        data["covariance"][0][0] = -1
    elif mutation == "hadron":
        data["points"][0]["hadron"] = "h+"
    else:
        data["datasets"][1]["replaces"] = [data["datasets"][0]["id"]]
    with pytest.raises(ValueError):
        validate_dataset(data)


def test_update_replaces_old_likelihood_and_replicas_require_refit():
    old = {"datasets": [{"id": "old", "overlap_group": "sample"}]}
    incoming = {"datasets": [{"id": "new", "overlap_group": "sample", "replaces": ["old"]}]}
    result = plan_update(old, incoming)
    assert result["active"] == ["new"]
    assert result["remove"] == ["old"]
    assert result["action"] == "refit"
    new = {"datasets": [{"id": "other", "overlap_group": "other"}]}
    assert plan_update(old, new, "posterior")["action"] == "posterior_reweight_candidate"
    assert plan_update(old, new, "replicas")["action"] == "refit"
    assert plan_update(old, new, "posterior", False)["action"] == "refit"


def test_profiled_nuisance_matches_effective_covariance():
    data = fixture()
    C = np.asarray(data["covariance"])
    B = np.asarray(data["nuisance_design"])
    r = np.linspace(-0.01, 0.01, len(C))
    a = np.linalg.solve(np.eye(B.shape[1]) + B.T @ np.linalg.solve(C, B), B.T @ np.linalg.solve(C, r))
    shifted = r - B @ a
    profiled = shifted @ np.linalg.solve(C, shifted) + a @ a
    assert profiled == pytest.approx(r @ np.linalg.solve(C + B @ B.T, r), rel=1e-12)

