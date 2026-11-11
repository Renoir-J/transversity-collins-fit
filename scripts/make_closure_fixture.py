"""Generate explicitly synthetic fixtures using an independent Python expression.

These are algebra/optimizer fixtures, not a phenomenological model or experimental data.
"""
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    target = ROOT / "fixtures/closure"
    target.mkdir(parents=True, exist_ok=True)
    truth = [0.72, -0.28, 0.24, -0.19]
    theta = [math.log(truth[0]), *truth[1:]]
    rng = np.random.default_rng(20261111)
    points, bins, golden, datasets = [], [], [], []
    for process in ("SIDIS", "SIA", "STAR_200", "STAR_510", "lattice"):
        datasets.append({"id": "synthetic_" + process, "publication": "NONE: synthetic closure fixture",
                         "overlap_group": "synthetic_" + process})
        for j in range(16 if process != "lattice" else 2):
            x = (j + 1) / 18
            # Exact algebra known independently of the JSON kernel evaluator.
            if process == "SIA":
                value = (1 + x) * truth[2]**2 - 0.6 * x * truth[2] * truth[3] + x**2 * truth[3]**2
                terms = [(3, 3, 1 + x), (3, 4, -0.6 * x), (4, 4, x*x)]
                linear = []
            elif process == "lattice":
                value = truth[0] - (1 + j) * truth[1]
                terms, linear = [], [{"index": 1, "weight": 1}, {"index": 2, "weight": -(1+j)}]
            else:
                value = (1+x) * truth[0] * truth[2] + (0.5-x) * truth[1] * truth[3] + x*x * truth[0] * truth[3]
                terms = [(1, 3, 1+x), (2, 4, 0.5-x), (1, 4, x*x)]
                linear = []
            denominator = 1.0 + x
            value /= denominator
            point_id = f"synthetic_{process}:{j}"
            points.append({"id": point_id, "dataset_id": "synthetic_" + process, "process": process,
                           "hadron": "pi+" if process != "lattice" else "none", "unit": "fraction",
                           "convention": "SYNTHETIC_ALGEBRA_V1", "value": value,
                           "bins": {"fixture_coordinate": [x, x]}})
            bins.append({"point_id": point_id, "sign": 1, "denominator": denominator,
                         "linear": linear,
                         "bilinear": [{"left": a, "right": b, "weight": w} for a, b, w in terms]})
            golden.append(value)
    n = len(points)
    sigma = np.array([0.01 if p["process"] != "lattice" else 0.035 for p in points])
    covariance = np.diag(sigma*sigma)
    nuisance = np.zeros((n, 2))
    for i, point in enumerate(points):
        nuisance[i, 0] = 0.003 if point["process"] == "SIDIS" else 0
        nuisance[i, 1] = 0.002 if point["process"] == "STAR_510" else 0
    noise = rng.multivariate_normal(np.zeros(n), covariance + nuisance @ nuisance.T)
    for point, shift in zip(points, noise):
        point["value"] += float(shift)
    data = {"schema_version": 1, "scope": "synthetic_closure", "datasets": datasets,
            "points": points, "covariance": covariance.tolist(), "nuisance_design": nuisance.tolist()}
    kernel = {"schema_version": 1, "scope": "synthetic_closure", "theory_version": "synthetic-algebra-v1",
              "physics_qualification": "No real cross sections; workflow closure only",
              "parameters": [{"name": name, "transform": "exp" if i == 0 else "identity"}
                             for i, name in enumerate(["synthetic_hu", "synthetic_hd", "synthetic_Cfav", "synthetic_Cunf"])],
              "initial": [math.log(0.6), -0.2, 0.2, -0.15], "truth_theta": theta,
              "independent_predictions": golden, "bins": bins}
    for name, value in (("data.json", data), ("kernel.json", kernel)):
        (target / name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Wrote {n} explicitly synthetic points")


if __name__ == "__main__":
    main()

