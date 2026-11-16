"""Explicitly synthetic full-Snakemake failure/interruption/forgery controls."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from collinsflow.historical_sidis_jobs import FLAGS, verify_complete_job
from collinsflow.provenance import atomic_json, file_hash, read_json, source_fingerprint, utc_now

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--primary-cache", required=True)
    parser.add_argument("--native-build", required=True)
    parser.add_argument("--parton-build", required=True)
    parser.add_argument("--seed-root")
    parser.add_argument("--mode", choices=("failure", "timeout", "signal"), default="failure")
    args = parser.parse_args()
    if os.name != "posix":
        raise ValueError("run complete controller controls on Linux")
    folder = (ROOT/args.run_dir).resolve()
    if not folder.is_relative_to(ROOT/"runs"):
        raise ValueError("fresh ignored control run directory required")
    folder.mkdir(parents=True, exist_ok=False)
    full = folder/"full"
    runtime = ROOT/".runtime"
    runtime.mkdir(exist_ok=True)
    fake_dir = Path(tempfile.mkdtemp(prefix="sidis-complete-control-", dir=runtime))
    fake = fake_dir/"intentional-numerical-nonexecution"
    payload = {"status": "FAIL", "scope": "explicitly_synthetic_controller_control_not_a_physical_prediction",
               "numerical_predictions_executed": False, "control_mode": args.mode, **dict.fromkeys(FLAGS, False)}
    fake.write_text(f'''#!{sys.executable}
from pathlib import Path
import json,os,signal,sys,time
payload={payload!r}
output=Path(sys.argv[-1])
output.write_text(json.dumps(payload,indent=2)+'\\n',encoding='utf-8')
Path(str(output)+'.progress.json').write_text(json.dumps(payload,indent=2)+'\\n',encoding='utf-8')
with Path(__file__).with_suffix('.calls').open('a',encoding='utf-8') as stream: stream.write('called\\n')
if payload['control_mode']!='failure':
 def stopped(number,frame): raise SystemExit(0)
 signal.signal(signal.SIGINT,stopped)
 if payload['control_mode']=='signal':
  parent=Path('/proc')/str(os.getppid())
  owner=int(next(line.split(':',1)[1] for line in (parent/'status').read_text().splitlines() if line.startswith('PPid:')))
  command=(Path('/proc')/str(owner)/'cmdline').read_bytes()
  if b'run_historical_sidis_job.py' not in command: raise RuntimeError('refuse signal outside owned CLI controller')
  os.kill(owner,signal.SIGTERM)
 time.sleep(30)
raise SystemExit(2)
''', encoding="utf-8")
    fake.chmod(0o700)
    config = ["run_dir="+str(full), "cache="+str((ROOT/args.cache).resolve()),
              "primary_cache="+str((ROOT/args.primary_cache).resolve()), "native_build="+str((ROOT/args.native_build).resolve()),
              "parton_build="+str((ROOT/args.parton_build).resolve()), "offline=True", "julia="+str(fake)]
    if args.seed_root:
        config.append("seed_root="+str(Path(args.seed_root).resolve()))
    if args.mode == "timeout":
        config.append("timeout_seconds=2")
    command = ["snakemake", "--snakefile", "workflow/historical-sidis-jobs.smk", "--cores", "2", "--config", *config]
    record = {"status": "FAIL", "scope": "explicitly_synthetic_complete_scheduler_interruption_and_forgery_control",
              "started_at": utc_now(), "control_mode": args.mode, "source_fingerprint": source_fingerprint(ROOT),
              "control_script_sha256": file_hash(__file__), "fake_executable_sha256": file_hash(fake),
              "command": command, "numerical_predictions_executed": False, **dict.fromkeys(FLAGS, False)}
    try:
        with (folder/"workflow.log").open("wb") as stream:
            result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
        record["workflow_returncode"] = result.returncode
        if result.returncode == 0:
            raise ValueError("synthetic failed complete worker unexpectedly succeeded")
        first = full/"native123/p0001"
        report, progress, invocation = first/"report.json", first/"report.json.progress.json", first/"invocation.json"
        if read_json(report) != payload or read_json(progress) != payload:
            raise ValueError("scheduler lost or changed failed report/progress")
        inv = read_json(invocation)
        expected = {"failure": (2, 2, None), "timeout": (124, 0, "resource_timeout"), "signal": (143, 0, "parent_signal")}[args.mode]
        if (inv["status"] != "FAIL" or "completed_at" not in inv or len(inv["steps"]) != 1 or
                [(s["returncode"], s["process_returncode"], s["interruption"]) for s in inv["steps"]] != [expected]):
            raise ValueError("failed/interrupted complete-command ledger was not retained")
        if inv["steps"][0]["resource"]["child_exit_code"] != (2 if args.mode == "failure" else 0):
            raise ValueError("raw GNU-time child outcome differs from owned control")
        if (first/"job.completed").exists() or (first/"audit.json").exists():
            raise ValueError("failed producer generated a completion marker or dispatched audit")
        originals = [report, progress, invocation, first/"initial-invocation.json", first/"producer.log", first/"producer.resources.txt"]
        hashes = {str(p): file_hash(p) for p in originals}
        rejected = False
        (first/"job.completed").touch()
        try:
            try:
                verify_complete_job(full/"request.json", args.native_build, args.parton_build, full, "native123", 1, julia=str(fake))
            except (ValueError, FileNotFoundError, KeyError):
                rejected = True
        finally:
            (first/"job.completed").unlink()
        if not rejected:
            raise ValueError("failed predecessor accepted a forged marker")
        forged = [full/phase/f"p{i:04d}"/"job.completed" for phase in ("native123", "tailref123") for i in range(1, 133)]
        forged.extend(full/phase/"job.completed" for phase in ("native-select3", "native-selection", "source-math", "derivative-select3", "diagnostic"))
        for marker in forged:
            marker.parent.mkdir(parents=True, exist_ok=True)
            if marker.exists():
                raise ValueError("refuse to replace an original marker during forgery control")
            marker.touch()
        forged_report = full/"diagnostic/report.json"
        claim = {"status": "CONDITIONAL_FULL132_POINT_DIAGNOSTIC_RECORDED_NOT_FITTING",
                 "scope": "explicitly_synthetic_forged_completion_claim", "source_fingerprint": record["source_fingerprint"],
                 "source_ids": [p["id"] for p in read_json(full/"request.json")["points"]], **dict.fromkeys(FLAGS, False)}
        claim["accepted_prediction_diagnostic"] = claim["all132_coverage_verified"] = True
        atomic_json(forged_report, claim)
        atomic_json(folder/"forged-report-retained.json", claim)
        try:
            with (folder/"forged-workflow.log").open("wb") as stream:
                forged_result = subprocess.run([*command, "--ignore-incomplete"], cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
            record["forged_workflow_returncode"] = forged_result.returncode
            if forged_result.returncode == 0:
                raise ValueError("full scheduler accepted forged all132 markers")
            try:
                verify_complete_job(full/"request.json", args.native_build, args.parton_build, full, "diagnostic", julia=str(fake))
            except (ValueError, FileNotFoundError, KeyError):
                pass
            else:
                raise ValueError("direct full132 verifier accepted a forged claim")
        finally:
            for p in [*forged, forged_report]:
                assert p.resolve().is_relative_to(full.resolve())
                p.unlink(missing_ok=True)
        if fake.with_suffix(".calls").read_text(encoding="utf-8").splitlines() != ["called"]:
            raise ValueError("forgery caused additional numerical dispatch")
        if any(file_hash(p) != h for p, h in hashes.items()) or source_fingerprint(ROOT) != record["source_fingerprint"]:
            raise ValueError("original failure evidence or source changed during controls")
        record.update(status="COMPLETE_SCHEDULER_FAILURE_INTERRUPTION_AND_FORGERY_CONTROL_VERIFIED",
                      original_failed_artifact_hashes=hashes, original_failed_artifacts_retained=True,
                      all132_marker_forgery_rejected=True, direct_cached_verifiers_rejected=True,
                      numerical_audit_dispatched=False, additional_numerical_dispatch_after_forgery=False,
                      raw_child_exit_code=inv["steps"][0]["resource"]["child_exit_code"],
                      effective_wrapper_exit_code=inv["steps"][0]["returncode"])
    except BaseException as error:
        record["error"] = type(error).__name__+": "+str(error)
        raise
    finally:
        record["completed_at"] = utc_now()
        atomic_json(folder/"control.json", record)
    print(json.dumps({k: record[k] for k in ("status", "control_mode", "numerical_predictions_executed")}))


if __name__ == "__main__":
    main()
