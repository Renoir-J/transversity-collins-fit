"""Explicit units/conventions, covariance and overlap validation."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from .provenance import atomic_json, file_hash, read_json

PROCESSES = {"SIDIS", "SIA", "STAR_200", "STAR_510", "lattice"}


def validate_dataset(data):
    if data.get("schema_version") != 1 or data.get("scope") not in {"synthetic_closure", "production"}:
        raise ValueError("unsupported data schema or scope")
    points = data["points"]
    if not points or len({p["id"] for p in points}) != len(points):
        raise ValueError("empty data or duplicate point ID")
    datasets = {d["id"]: d for d in data["datasets"]}
    if len(datasets) != len(data["datasets"]):
        raise ValueError("duplicate dataset ID")
    active = set(datasets)
    overlap = {}
    for dataset in datasets.values():
        if set(dataset.get("replaces", [])) & active:
            raise ValueError("a replacement and its superseded data are both active")
        sample = dataset.get("overlap_group")
        if sample:
            overlap.setdefault(sample, []).append(dataset)
    for group in overlap.values():
        if len(group) > 1 and not all(d.get("joint_covariance_resolved") for d in group):
            raise ValueError("overlapping projections require a joint covariance or disjoint selection")
    for point in points:
        if point["dataset_id"] not in datasets or point["process"] not in PROCESSES:
            raise ValueError("unknown dataset or process")
        if not math.isfinite(point["value"]) or point.get("unit") != "fraction":
            raise ValueError("values must be finite and normalized to fractions")
        if not point.get("convention") or not point.get("hadron"):
            raise ValueError("explicit convention and hadron are required")
        if point["hadron"] in {"h+", "h-"} and not point.get("mixture_adapter"):
            raise ValueError("unidentified hadrons require a mixture adapter")
        for bounds in point.get("bins", {}).values():
            if len(bounds) != 2 or not all(math.isfinite(x) for x in bounds) or bounds[0] > bounds[1]:
                raise ValueError("invalid bin bounds")
    covariance = np.asarray(data["covariance"], dtype=float)
    if covariance.shape != (len(points), len(points)) or not np.isfinite(covariance).all():
        raise ValueError("invalid covariance dimensions or entries")
    if not np.allclose(covariance, covariance.T, atol=1e-14, rtol=1e-12):
        raise ValueError("covariance is not symmetric")
    try:
        np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError("covariance must be positive definite; no silent diagonal jitter") from exc
    nuisance = np.asarray(data.get("nuisance_design", [[] for _ in points]), dtype=float)
    if nuisance.ndim != 2 or nuisance.shape[0] != len(points) or not np.isfinite(nuisance).all():
        raise ValueError("invalid additive nuisance design")
    if data["scope"] == "production":
        for dataset in datasets.values():
            for key in ("publication", "table", "source_url", "source_sha256", "covariance_status",
                        "selection_reason", "version"):
                if not dataset.get(key):
                    raise ValueError(f"production data missing {key}")
            if dataset["covariance_status"] not in {"published_full", "published_uncorrelated", "approved_model"}:
                raise ValueError("unresolved production covariance")
    return {"points": len(points), "datasets": len(datasets), "scope": data["scope"],
            "covariance_min_eigenvalue": float(np.linalg.eigvalsh(covariance).min())}


def normalize_table(source, mapping, output):
    """Normalize a reviewed CSV/XLSX mapping. No inferred signs or column meanings."""
    import csv

    source = Path(source)
    spec = read_json(mapping)
    if file_hash(source) != spec["source_sha256"]:
        raise ValueError("raw source checksum mismatch")
    if source.suffix.lower() == ".xlsx":
        from openpyxl import load_workbook
        workbook = load_workbook(source, read_only=True, data_only=True)
        sheet = workbook[spec["sheet"]]
        rows = list(sheet.iter_rows(values_only=True))
        header = rows[spec.get("header_row", 1) - 1]
        records = [dict(zip(header, row)) for row in rows[spec.get("header_row", 1):]]
    else:
        with source.open(encoding="utf-8-sig", newline="") as stream:
            records = list(csv.DictReader(stream))
    selected = []
    multiplier = {"fraction": 1.0, "percent": 0.01}[spec["input_unit"]]
    if spec["sign"] not in (-1, 1):
        raise ValueError("sign must be explicitly +1 or -1")
    for index, record in enumerate(records):
        if not all(str(record.get(key)) == str(value) for key, value in spec.get("select", {}).items()):
            continue
        point = {"id": f"{spec['dataset']['id']}:{index}", "dataset_id": spec["dataset"]["id"],
                 "process": spec["process"], "hadron": spec["hadron"], "unit": "fraction",
                 "convention": spec["convention"],
                 "value": float(record[spec["columns"]["value"]]) * multiplier * spec["sign"],
                 "bins": {key: [float(record[low]), float(record[high])]
                          for key, (low, high) in spec.get("bins", {}).items()}}
        sigma = float(record[spec["columns"]["stat"]]) * multiplier
        if not math.isfinite(sigma) or sigma <= 0:
            raise ValueError("nonpositive statistical error")
        selected.append((point, sigma))
    if not selected:
        raise ValueError("table selection produced no points")
    # Correlated systematics are supplied separately, never inferred as diagonal.
    n = len(selected)
    covariance = np.diag([sigma**2 for _, sigma in selected])
    if spec.get("covariance_file"):
        covariance = np.asarray(read_json(Path(mapping).parent / spec["covariance_file"])) * multiplier**2
    result = {"schema_version": 1, "scope": "production", "datasets": [spec["dataset"]],
              "points": [point for point, _ in selected], "covariance": covariance.tolist(),
              "nuisance_design": spec.get("nuisance_design", [[] for _ in range(n)])}
    summary = validate_dataset(result)
    atomic_json(output, result)
    return summary


def plan_update(current, incoming, sample_kind="replicas", unchanged_contract=True):
    old = {d["id"] for d in current["datasets"]}
    new = {d["id"] for d in incoming["datasets"]}
    if old & new:
        raise ValueError("dataset IDs already active; use a versioned replacement")
    replaced = {key for d in incoming["datasets"] for key in d.get("replaces", [])}
    if not replaced <= old:
        raise ValueError("replacement target is not active")
    combined = [d for d in current["datasets"] if d["id"] not in replaced] + incoming["datasets"]
    groups = [d.get("overlap_group") for d in combined if d.get("overlap_group")]
    if len(groups) != len(set(groups)) and not all(d.get("joint_covariance_resolved") for d in combined):
        raise ValueError("update introduces unresolved overlapping datasets")
    reweight = sample_kind == "posterior" and unchanged_contract and not replaced
    return {"remove": sorted(replaced), "add": sorted(new), "active": [d["id"] for d in combined],
            "action": "posterior_reweight_candidate" if reweight else "refit",
            "must_refit_on_low_ess": True, "invalidate": ["data", "benchmark", "kernels",
                "gradient_validation", "performance_validation", "closure", "fit", "ablate", "products"]}

