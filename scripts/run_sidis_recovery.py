"""Run the declared full132 SIDIS recovery without rewriting any earlier namespace."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

ADAPTER = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("sidis_recovery_transport", ADAPTER / "src/collinsflow/sidis_recovery.py")
transport = importlib.util.module_from_spec(spec)
spec.loader.exec_module(transport)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "plan", "review", "all", "verify", "snapshot"))
    parser.add_argument("--frozen-checkout", required=True)
    parser.add_argument("--predecessor-checkout", required=True)
    for name in ("manifest", "manifest-sha", "parent-manifest", "parent-sha", "inventory", "inventory-sha", "output"):
        parser.add_argument("--" + name)
    parser.add_argument("--julia", default="/home/renoir_j/.juliaup/bin/julia")
    args = parser.parse_args()
    frozen = Path(args.frozen_checkout)
    predecessor = Path(args.predecessor_checkout)
    transport.lineage.source_bindings(frozen, predecessor, ADAPTER)
    sys.path.insert(0, str(frozen / "src"))
    from collinsflow import historical_sidis_jobs as api
    if api.ROOT != frozen:
        raise ValueError("wrong explicitly verified frozen package")
    if args.action == "snapshot":
        if any(getattr(args, key) is None for key in ("parent_manifest", "parent_sha", "output")):
            parser.error("snapshot requires parent-manifest/parent-sha/output")
        result = transport.snapshot(args.parent_manifest, args.parent_sha, args.output)
    elif args.action == "prepare":
        if any(getattr(args, key) is None for key in ("parent_manifest", "parent_sha", "inventory", "inventory_sha", "output")):
            parser.error("prepare requires parent-manifest/parent-sha/inventory/inventory-sha/output")
        with transport.host.family_lock():
            path = transport.prepare(frozen, predecessor, ADAPTER, args.parent_manifest, args.parent_sha,
                args.inventory, args.inventory_sha, args.output, api, args.julia)
        result = dict(status="PREPARED_NOT_EXECUTED_NOT_A_GATE", manifest=str(path), manifest_sha256=transport.file_hash(path))
    else:
        if args.manifest is None or args.manifest_sha is None:
            parser.error("explicit manifest and externally bound manifest-sha required")
        path = Path(args.manifest)
        m = transport.verify_manifest(path, frozen, predecessor, ADAPTER, api, args.manifest_sha)
        if args.julia != m["runtime"]["launcher"]:
            raise ValueError("caller Julia launcher differs from the prepared runtime binding")
        if args.action == "plan":
            result = dict(status="PLAN_ONLY_NO_JULIA_DISPATCH", retained_jobs=len(m["retained_jobs"]),
                retained_native123=sum(row["phase"] == "native123" for row in m["retained_jobs"]),
                all132_common_native_tail_and_all13_reference_policy_still_required=True,
                **dict.fromkeys(transport.FLAGS, False))
        elif args.action == "all":
            with transport.driver(m, path, api, args.action) as lease:
                result = transport.run_all(m, path, api, lease)
        else:
            with transport.host.family_lock() as lease:
                result = transport.audit_retained(m, path, api, lease) if args.action == "review" else transport.audit_cached(m, path, "diagnostic", 0, api, lease)
    print(json.dumps({key: result[key] for key in ("status", "manifest", "manifest_sha256", "retained_jobs", "retained_native123", *transport.FLAGS) if key in result}), flush=True)


if __name__ == "__main__":
    main()
