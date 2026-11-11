"""Validate the research DAG and distinguish closure implementation from physics readiness."""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def ordered_stages(spec, target=None):
    stages = {stage["id"]: stage for stage in spec["stages"]}
    if len(stages) != len(spec["stages"]):
        raise ValueError("duplicate stage ID")
    order, active, done = [], set(), set()

    def visit(key):
        if key not in stages:
            raise ValueError(f"unknown dependency or target {key}")
        if key in active:
            raise ValueError(f"cycle at {key}")
        if key in done:
            return
        active.add(key)
        for dependency in stages[key]["needs"]:
            visit(dependency)
        active.remove(key)
        done.add(key)
        order.append(stages[key])

    for key in [target] if target else stages:
        visit(key)
    return order


def validate(spec, registry, numerics):
    order = ordered_stages(spec)
    gates = {gate["id"] for gate in spec["gates"]}
    for stage in order:
        if stage.get("gate") and stage["gate"] not in gates:
            raise ValueError(f"missing gate for {stage['id']}")
    scenarios = {scenario["id"] for scenario in spec["scenarios"]}
    if len(scenarios) != len(spec["scenarios"]):
        raise ValueError("duplicate scenario ID")
    for scenario in spec["scenarios"]:
        if scenario.get("parent") and scenario["parent"] not in scenarios:
            raise ValueError(f"missing parent for {scenario['id']}")
    ids = {dataset["id"] for dataset in registry["datasets"]}
    if len(ids) != len(registry["datasets"]):
        raise ValueError("duplicate dataset ID")
    for dataset in registry["datasets"]:
        if dataset.get("replaces") and dataset["replaces"] not in ids:
            raise ValueError(f"unknown replaced dataset in {dataset['id']}")
    # Prevent the design scaffold from being misrepresented as an executed fit.
    if spec["status"] not in {"DESIGN_ONLY", "CLOSURE_IMPLEMENTED_PRODUCTION_BLOCKED"} or registry["normalized_for_production"]:
        raise ValueError("implementation/data readiness needs a separate audited change")
    if not spec["update_policy"]["replicas_are_not_posterior_samples"]:
        raise ValueError("replica and posterior distinction is required")
    if spec["execution"]["agent_calls_in_optimizer"] != 0:
        raise ValueError("optimizer numerical path must be deterministic")
    if not spec["receipts"]["fail_stops_descendants"]:
        raise ValueError("failed stages must stop their descendants")
    architecture = numerics["architecture"]
    if (spec["execution"]["architecture_decision"] != numerics["decision_id"]
            or spec["execution"]["numerical_backend"] != architecture["numerical_language"]
            or spec["execution"]["coordinator"] != architecture["orchestration_language"]
            or spec["execution"]["deterministic_runner"] != architecture["dag_runner"]):
        raise ValueError("workflow and numerical architecture disagree")
    if (architecture["numerical_language"] != "Julia"
            or architecture["orchestration_language"] != "Python"
            or architecture["dag_runner"] != "Snakemake"
            or architecture["numerical_python_callbacks"] != 0):
        raise ValueError("ADR-001 requires a Julia numerical loop and Python/Snakemake job orchestration")
    if numerics["implementation_status"] not in {"DESIGN_ONLY", "CLOSURE_IMPLEMENTED_PRODUCTION_BLOCKED"}:
        raise ValueError("numerical implementation readiness needs a separate audited change")
    required_gates = {"benchmark": "G2", "gradient_validation": "G2AD",
                      "performance_validation": "G2PERF", "closure": "G3"}
    stage_by_id = {stage["id"]: stage for stage in order}
    fit_ancestors = {stage["id"] for stage in ordered_stages(spec, "fit")}
    for stage_id, gate in required_gates.items():
        if stage_id not in fit_ancestors or stage_by_id[stage_id].get("gate") != gate:
            raise ValueError(f"fit must depend on {stage_id} and its {gate} gate")
    if "gradient_validation" not in {
            stage["id"] for stage in ordered_stages(spec, "performance_validation")}:
        raise ValueError("performance acceptance must follow gradient validation")
    produced = {output for stage in order for output in stage["outputs"]}
    missing = set(numerics["validation"]["required_reports"]) - produced
    if missing:
        raise ValueError(f"validation reports absent from DAG: {sorted(missing)}")
    if numerics["fortran_policy"]["upstream_checkout_mutation"]:
        raise ValueError("the pinned upstream Fortran reference must remain unchanged")
    inventory = json.loads((ROOT / "references/source-inventory.json").read_text(encoding="utf-8"))
    indexed_sources = {"external/" + entry["repo"] + "/" + entry["path"] for entry in inventory["source_files"]}
    for module in numerics["fortran_policy"]["module_actions"]:
        for source in module["sources"]:
            if not (ROOT / source).is_file() and source not in indexed_sources:
                raise ValueError(f"unknown Fortran reference source: {source}")
    for key in ("architecture_decision", "backend_versions", "compiler_and_flags",
                "precision", "thread_and_process_counts", "validation_policy_hash"):
        if key not in spec["receipts"]["required"]:
            raise ValueError(f"numerical provenance missing from receipts: {key}")
    return {"design_valid": True, "stages": len(order), "scenarios": len(scenarios),
            "dataset_candidates": len(ids), "architecture_decision": numerics["decision_id"],
            "numerical_backend": architecture["numerical_language"],
            "dag_runner": architecture["dag_runner"],
            "production_fit_implemented": False, "production_data_normalized": False,
            "physics_gradient_performance_checks_executed": False,
            "closure_engine_implemented": spec["status"] == "CLOSURE_IMPLEMENTED_PRODUCTION_BLOCKED"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["validate", "plan"], nargs="?", default="plan")
    parser.add_argument("--target", default="products")
    args = parser.parse_args()
    spec = json.loads((ROOT / "config/workflow.json").read_text(encoding="utf-8"))
    registry = json.loads((ROOT / "config/data-registry.json").read_text(encoding="utf-8"))
    numerics = json.loads((ROOT / spec["execution"]["numerical_specification"]).read_text(encoding="utf-8"))
    validation = validate(spec, registry, numerics)
    if args.action == "validate":
        print(json.dumps(validation, indent=2))
        return
    print("Planning only: this command launches no jobs; use collinsflow or Snakemake.\n")
    print(f"Architecture: {validation['numerical_backend']} numerical core; "
          f"Python/{validation['dag_runner']} orchestration.\n")
    for stage in ordered_stages(spec, args.target):
        print(f"{stage['id']:12} <- {', '.join(stage['needs']) or '(start)'}")
        print("  produces: " + ", ".join(stage["outputs"]))
        if stage.get("gate"):
            print("  requires: " + stage["gate"])


if __name__ == "__main__":
    main()
