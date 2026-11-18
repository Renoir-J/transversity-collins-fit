"""Read the pinned24 predecessor pairs using source blobs verified at this job boundary."""
from pathlib import Path


def input_hashes(parent, manifest, index, verified_source_files, old, api):
    if (type(index) is not int or not 14 <= index <= 37 or len(verified_source_files) != 348 or
            dict(verified_source_files) != parent["adapter_files"]):
        raise ValueError("exact current348-blob predecessor proof and declared pair14..37 required")
    # These exact canonical source paths and bytes were all checked against raw
    # Git blobs in the same whole-job boundary. Do not resolve/read them348 times
    # for each old pair; they are not physical query or parameter-response data.
    hashes = {str(Path(parent["adapter"]) / name): sha for name, sha in dict(verified_source_files).items()}
    frozen = Path(parent["frozen"])
    extra = [Path(parent[key]) for key in ("request", "native", "parton")]
    extra += [frozen / "config/historical-sidis-job-policy.json", api.TIME_BINARY,
              frozen / "scripts/run_historical_sidis_job.py", Path(manifest),
              old.point_folder(parent, "native123", index - 1) / "invocation.json",
              Path(parent["run"]) / "prefix-review/invocation.json"]
    hashes.update({str(path): old.file_hash(path) for path in extra})
    return hashes


def verify(parent, manifest, index, verified_source_files, old, api):
    """Equivalent complete-pair fields/resources; all old packets remain at their actual paths."""
    folder = old.point_folder(parent, "native123", index)
    inv = old.read_json(folder / "invocation.json")
    expected_provenance = old.provenance(parent, manifest)
    if (inv.get("status") != api.COMPLETED or inv.get("phase") != "native123" or
            type(inv.get("point_index")) is not int or inv["point_index"] != index or
            inv.get("continuation") != expected_provenance or inv.get("source_fingerprint") != parent["numerical_source_fingerprint"] or
            inv.get("policy_sha256") != api.POLICY_SHA or not api._flags(inv) or
            inv.get("time_binary") != str(api.TIME_BINARY) or not inv.get("completed_at")):
        raise ValueError("missing, failed, stale or self-admitting predecessor job")
    if (inv["input_hashes"] != input_hashes(parent, manifest, index, verified_source_files, old, api) or
            set(inv["outputs"]) != api._products(folder) or any(old.file_hash(path) != sha for path, sha in inv["outputs"].items()) or
            old.file_hash(folder / "initial-invocation.json") != inv["initial_record_sha256"] or
            old.read_json(folder / "initial-invocation.json")["record_role"] != "initial_sentinel_not_an_executed_failure" or
            [step["name"] for step in inv["steps"]] != ["producer", "audit"]):
        raise ValueError("changed or incomplete predecessor inputs, outputs, sentinel or commands")
    budget = {"resources": api.load_policy()["resources"]["native123"]}
    for step in inv["steps"]:
        if any(type(step.get(key)) is not int or step[key] != 0 for key in ("returncode", "process_returncode")):
            raise ValueError("inexact predecessor command outcomes")
        api._require_command_success(step, budget)
        name = step["name"]
        if (old.file_hash(folder / (name + ".log")) != step["log_sha256"] or
                old.file_hash(folder / (name + ".resources.txt")) != step["resource_sha256"] or
                api.parse_resource_file(folder / (name + ".resources.txt")) != step["resource"]):
            raise ValueError("changed complete predecessor command measurements")
    if sum(step["wrapper_elapsed_seconds"] for step in inv["steps"]) > budget["resources"]["whole_Julia_jobs_time_budget_seconds"]:
        raise ValueError("predecessor whole-job resource budget exhausted")
    report, audit = api._report_audit(folder, "native123", index, old.read_json(parent["request"]),
                                     parent["numerical_source_fingerprint"], parent["request_sha256"])
    if any(item.get("continuation") != expected_provenance for item in (report, audit)):
        raise ValueError("unbound predecessor producer or auditor")
    marker = folder / "job.completed"
    if not marker.is_file() or marker.is_symlink():
        raise ValueError("missing canonical predecessor scheduling marker")
    return inv
