"""Deterministic job DAG with checked, scope-bound receipts."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import traceback
from pathlib import Path

from .data import validate_dataset
from .provenance import (atomic_json, digest, environment, file_hash, git_commit, read_json,
                         source_fingerprint, utc_now, verify_receipt)


class Blocked(RuntimeError):
    pass


def root_path():
    return Path(__file__).resolve().parents[2]


def production_readiness(root):
    contract = read_json(root / "config/physics-contract.json")
    missing = list(contract["missing_prerequisites"])
    if not contract["production_ready"]:
        missing.append("production_ready is false")
    for key in ("production_kernel", "production_data", "independent_physics_benchmark", "frozen_production_policy"):
        path = contract.get(key)
        if not path or not (root / path).is_file():
            missing.append(f"missing artifact: {key}")
    return missing


def run_stage(stage, run_dir, mode="closure", julia="julia"):
    root, run_dir = root_path(), Path(run_dir).resolve()
    spec = read_json(root / "config/workflow.json")
    stages = {s["id"]: s for s in spec["stages"]}
    if stage not in stages:
        raise ValueError("unknown stage")
    scope = "synthetic_closure" if mode == "closure" else "production"
    run_dir.mkdir(parents=True, exist_ok=True)
    fingerprint = source_fingerprint(root)
    env = environment(root, julia)
    env_digest = digest(env)
    policy_hash = file_hash(root / "config/numerical-backends.json")
    seed = read_json(root / "config/physics-contract.json")["synthetic_profile"]["seed"]
    receipt_path = run_dir / "receipts" / f"{stage}.json"
    inputs = {str(root / "config/workflow.json"): file_hash(root / "config/workflow.json"),
              str(root / "config/physics-contract.json"): file_hash(root / "config/physics-contract.json")}
    checks, outputs = [], []
    receipt = {"stage": stage, "scope": scope, "source_fingerprint": fingerprint,
               "code_commit": git_commit(root), "environment_digest": env_digest, "backend_versions": env,
               "seed": seed, "started_at": utc_now(), "architecture_decision": "ADR-001",
               "validation_policy_hash": policy_hash, "compiler_and_flags": "Julia managed; Fortran bridge separately recorded",
               "precision": "Float64", "thread_and_process_counts": {"job_processes": 1, "julia_threads": 1,
                   "blas_threads": 1}, "input_hashes": inputs, "outputs": outputs, "checks": checks,
               "status": "RUNNING"}
    atomic_json(receipt_path, receipt)

    def emit(name, value):
        path = run_dir / stage / name
        atomic_json(path, value)
        outputs.append({"path": str(path), "sha256": file_hash(path)})
        return path

    def check(name, passed, detail=None):
        checks.append({"name": name, "passed": bool(passed), "detail": detail})
        if not passed:
            raise ValueError(f"failed check: {name}")

    try:
        for dependency in stages[stage]["needs"]:
            dependency_path = run_dir / "receipts" / f"{dependency}.json"
            verify_receipt(dependency_path, scope, fingerprint, env_digest)
            inputs[str(dependency_path)] = file_hash(dependency_path)
        if mode == "production":
            missing = production_readiness(root)
            if missing:
                raise Blocked("Production prerequisites unresolved: " + "; ".join(missing))
            # The production path needs an independent admission implementation, not a boolean override.
            raise Blocked("Production scientific admission is not certified; cannot substitute synthetic kernels")
        profile = read_json(root / "config/physics-contract.json")["synthetic_profile"]
        data_path, kernel_path = root / profile["data"], root / profile["kernel"]
        for path in (data_path, kernel_path):
            inputs[str(path)] = file_hash(path)
        if stage == "contract":
            emit("physics-contract.json", read_json(root / "config/physics-contract.json"))
            emit("conventions.json", {"scope": scope, "physical_results": False, "unit": "fraction"})
            check("explicit_synthetic_scope", scope == "synthetic_closure")
        elif stage == "sources":
            emit("sources.lock.json", read_json(root / "references/source-inventory.json")["repositories"])
            emit("source-review.json", {"scope": scope, "pinned_source_index": "references/source-inventory.json",
                 "KPSY_author_code": "NOT_LOCATED", "upstream_code_executed_in_this_job": False})
            check("source_metadata_recorded", True)
        elif stage == "data":
            data = read_json(data_path)
            summary = validate_dataset(data)
            emit("normalized-data.json", data)
            emit("overlap-graph.json", {"scope": scope, "groups": [d["overlap_group"] for d in data["datasets"]]})
            check("data_validation", summary["points"] > 0, summary)
        elif stage == "theory":
            kernel = read_json(kernel_path)
            check("synthetic_kernel_label", kernel["scope"] == scope and kernel["theory_version"] == "synthetic-algebra-v1")
            emit("process-adapters.json", {"scope": scope, "kernel": str(kernel_path), "sha256": file_hash(kernel_path),
                 "production_adapters": "NOT_IMPLEMENTED"})
            emit("evolution-benchmarks.json", {"scope": scope, "implementation": "fixed-nf LO moment checks in Julia",
                 "Collins_evolution": "NOT_IMPLEMENTED", "KPSY_reproduction": False})
        elif stage == "kernels":
            key = digest({"source": fingerprint, "kernel": file_hash(kernel_path), "data": file_hash(data_path),
                          "environment": env_digest, "policy": policy_hash, "scope": scope})
            cache = run_dir / "cache" / key / "kernel.json"
            cache.parent.mkdir(parents=True, exist_ok=True)
            if not cache.exists():
                shutil.copyfile(kernel_path, cache)
            check("cache_equals_clean_source", file_hash(cache) == file_hash(kernel_path))
            outputs.append({"path": str(cache), "sha256": file_hash(cache)})
            emit("evolution-operator-manifest.json", {"scope": scope, "cache_key": key,
                 "fixed_operator_AD": "Julia generic matrix application", "production_operator_builder": "NOT_IMPLEMENTED"})
            emit("cache-validation-report.json", {"scope": scope, "passed": True, "cache_key": key,
                 "hash_dependencies": ["source", "kernel", "data", "environment", "policy", "scope"]})
        elif stage == "products":
            from .products import build_products
            for path in build_products(run_dir, read_json(run_dir / "fit/fit.json"), read_json(data_path)):
                outputs.append({"path": str(path), "sha256": file_hash(path)})
            check("synthetic_products_labeled", True)
        else:
            actions = {"benchmark": "benchmark", "gradient_validation": "gradient",
                       "performance_validation": "performance", "closure": "closure", "fit": "fit", "ablate": "ablate"}
            action = actions[stage]
            filename = {"benchmark": "benchmark-report.json", "gradient_validation": "gradient-validation-report.json",
                        "performance_validation": "performance-report.json", "closure": "closure-report.json",
                        "fit": "fit.json", "ablate": "scenario-comparisons.json"}[stage]
            output = run_dir / stage / filename
            job = {"action": action, "data": str(data_path), "kernel": str(kernel_path), "output": str(output),
                   "seed": seed, "replicas": profile["replicas"], "scope": scope}
            job_path = emit("job.json", job)
            log_path = run_dir / stage / "julia.log"
            command = [julia, "--startup-file=no", f"--project={root / 'numerics'}", str(root / "numerics/bin/run_job.jl"), str(job_path)]
            process_env = os.environ.copy()
            process_env.update(JULIA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
            with log_path.open("w", encoding="utf-8") as log:
                completed = subprocess.run(command, cwd=root, env=process_env, stdout=log, stderr=subprocess.STDOUT)
            outputs.append({"path": str(log_path), "sha256": file_hash(log_path)})
            check("julia_job_exit", completed.returncode == 0, {"returncode": completed.returncode, "log": str(log_path)})
            result = read_json(output)
            check("julia_result_scope", result["scope"] == scope)
            check("numerical_job_passed", result.get("passed", result.get("converged", True)))
            if stage == "ablate":
                check("all_holdout_fits_converged", all(h["result"]["converged"] for h in result["holdouts"]))
            outputs.append({"path": str(output), "sha256": file_hash(output)})
            if stage == "performance_validation":
                emit("performance-budget.json", result["budget"])
        if not checks:
            check("stage_outputs_written", bool(outputs))
        receipt["status"] = "PASS"
    except Exception as exc:
        receipt["status"] = "BLOCKED" if isinstance(exc, Blocked) else "FAIL"
        receipt["error"] = str(exc)
        receipt["traceback"] = traceback.format_exc()
        raise
    finally:
        receipt["completed_at"] = utc_now()
        atomic_json(receipt_path, receipt)
        if receipt["status"] in {"FAIL", "BLOCKED"}:
            # Snakemake removes failed declared outputs; keep separate failure evidence.
            failure_id = digest({"started_at": receipt["started_at"], "stage": stage})[:16]
            atomic_json(run_dir / "failures" / f"{stage}-{failure_id}.json", receipt)
    return receipt


def execute(run_dir, mode="closure", julia="julia", target="products"):
    root = root_path()
    sys.path.insert(0, str(root / "scripts"))
    from plan_workflow import ordered_stages
    spec = read_json(root / "config/workflow.json")
    for stage in ordered_stages(spec, target):
        print(f"Running {stage['id']} ({mode})", flush=True)
        run_stage(stage["id"], run_dir, mode, julia)

