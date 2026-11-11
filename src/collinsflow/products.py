"""Export labeled closure predictions and uncertainty figures."""

import csv

import matplotlib
import numpy as np

from .provenance import atomic_json, read_json

matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402


def build_products(run_dir, fit, data):
    target = run_dir / "products"
    target.mkdir(parents=True, exist_ok=True)
    paths = []
    csv_path = target / "per-bin-predictions.csv"
    sigma = np.sqrt(np.diag(data["covariance"]))
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["scope", "point_id", "process", "synthetic_observation", "prediction", "uncorrelated_sigma"])
        for point, prediction, error in zip(data["points"], fit["predictions"], sigma):
            writer.writerow([fit["scope"], point["id"], point["process"], point["value"], prediction, error])
    paths.append(csv_path)
    fig, ax = plt.subplots(figsize=(10, 4.6), layout="constrained")
    index = np.arange(len(sigma))
    ax.errorbar(index, [p["value"] for p in data["points"]], yerr=sigma, fmt=".", label="synthetic pseudo-data")
    ax.plot(index, fit["predictions"], ".", label="Julia fit")
    ax.set(xlabel="Synthetic bin index", ylabel="Fixture observable (fraction)",
           title="Synthetic closure only — no extracted physical distributions")
    ax.legend()
    figure = target / "closure-predictions.png"
    fig.savefig(figure, dpi=150)
    plt.close(fig)
    paths.append(figure)
    replicas = np.asarray(fit["replicas"])
    covariance = np.cov(replicas, rowvar=False).tolist() if len(replicas) > 1 else None
    summary = {"scope": fit["scope"], "physical_extraction": False, "chi2": fit["chi2"],
               "fit_parameters": fit["coefficients"], "replicas": len(replicas),
               "sample_kind": fit["sample_kind"], "unconstrained_parameter_replica_covariance": covariance,
               "production_h1_Collins_grids": "BLOCKED: process-accurate kernels not admitted",
               "tensor_charge": "NOT_EXTRACTED_FROM_REAL_DATA",
               "closure": read_json(run_dir / "closure/closure-report.json"),
               "receipts": "../receipts", "limitations": fit["limitations"]}
    report = target / "summary.json"
    atomic_json(report, summary)
    paths.append(report)
    return paths

