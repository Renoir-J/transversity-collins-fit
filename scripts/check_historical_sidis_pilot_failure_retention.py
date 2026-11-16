"""Owned synthetic scheduler failure and forged-marker controls; no predictions."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from collinsflow.historical_sidis_pilot import FALSE_FLAGS, load_policy, verify_pilot_job
from collinsflow.provenance import atomic_json, file_hash, git_commit, read_json, source_fingerprint, utc_now

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', required=True)
    parser.add_argument('--cache', required=True)
    parser.add_argument('--primary-cache', required=True)
    parser.add_argument('--native-build', required=True)
    parser.add_argument('--parton-build', required=True)
    parser.add_argument('--seed-root')
    args = parser.parse_args()
    if os.name != 'posix':
        raise ValueError('run the owned scheduler control on Linux')
    folder = (ROOT/args.run_dir).resolve()
    if not folder.is_relative_to(ROOT/'runs'):
        raise ValueError('use a fresh ignored runs directory')
    folder.mkdir(parents=True, exist_ok=False)
    full = folder/'full'
    runtime = ROOT/'.runtime'
    runtime.mkdir(exist_ok=True)
    fake_dir = Path(tempfile.mkdtemp(prefix='sidis-pilot-failure-control-', dir=runtime))
    fake = fake_dir/'intentional-failing-julia'
    payload = dict(status='FAIL', scope='explicitly_synthetic_pilot_failure_retention_control',
                   numerical_predictions_executed=False, **{k: False for k in FALSE_FLAGS})
    fake.write_text(f'''#!{sys.executable}
from pathlib import Path
import json,sys
payload={payload!r}
Path(sys.argv[-1]).write_text(json.dumps(payload)+"\\n")
Path(sys.argv[-1]+".progress.json").write_text(json.dumps(payload)+"\\n")
with Path(__file__).with_suffix(".calls").open("a") as stream:stream.write("called\\n")
raise SystemExit(2)
''', encoding='utf-8')
    fake.chmod(0o700)
    config = ['run_dir='+str(full), 'cache='+str(Path(args.cache).resolve()),
              'primary_cache='+str(Path(args.primary_cache).resolve()),
              'native_build='+str(Path(args.native_build).resolve()),
              'parton_build='+str(Path(args.parton_build).resolve()), 'julia='+str(fake), 'offline=True']
    if args.seed_root:
        config.append('seed_root='+str(Path(args.seed_root).resolve()))
    command = ['snakemake', '--snakefile', str(ROOT/'workflow/historical-sidis-pilot.smk'), '--cores', '1', '--config', *config]
    record = dict(status='FAIL', scope='explicitly_synthetic_scheduler_failure_and_marker_control',
        **{k: False for k in FALSE_FLAGS}, numerical_predictions_executed=False, started_at=utc_now(),
        source_fingerprint=source_fingerprint(ROOT), code_commit=git_commit(ROOT),
        control_script_sha256=file_hash(__file__), fake_executable_sha256=file_hash(fake), command=command)
    try:
        with (folder/'workflow.log').open('w', encoding='utf-8') as stream:
            result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
        record['workflow_returncode'] = result.returncode
        if result.returncode == 0:
            raise ValueError('intentional failed producer unexpectedly succeeded')
        cases = load_policy(ROOT)['selection']['cases']
        diagnostic = full/'cases'/cases[0]['id']
        report, progress = diagnostic/'pilot.json', diagnostic/'pilot.json.progress.json'
        invocation = diagnostic/'invocation.json'
        if read_json(report) != payload or read_json(progress) != payload:
            raise ValueError('failed synthetic report/progress lost')
        inv = read_json(invocation)
        if (inv['status'] != 'FAIL' or len(inv['steps']) != 1 or inv['steps'][0]['name'] != 'producer' or
                inv['steps'][0]['returncode'] != 2 or inv['steps'][0]['process_returncode'] != 2 or
                inv['steps'][0]['resource']['child_exit_code'] != 2 or inv['steps'][0]['interruption'] is not None):
            raise ValueError('complete-command failure ledger not preserved')
        if (diagnostic/'job.completed').exists() or (diagnostic/'audit.json').exists() or (full/'pilot-ledger.json').exists():
            raise ValueError('failed worker produced an acceptance artifact')
        retained = [report, progress, invocation, diagnostic/'initial-invocation.json',
                    diagnostic/'producer.log', diagnostic/'producer.resources.txt']
        hashes = {str(p): file_hash(p) for p in retained}
        stamps = [full/'cases'/case['id']/'job.completed' for case in cases]
        forged = full/'pilot-ledger.json'
        forgery = dict(status='PILOT_COMPLETED_NOT_DIAGNOSTIC_GATE', scope='explicitly_synthetic_forged_marker_control',
                       **{k: False for k in FALSE_FLAGS})
        atomic_json(folder/'forged-ledger.json', forgery)
        rejected = False
        try:
            for stamp in stamps:
                stamp.parent.mkdir(parents=True, exist_ok=True)
                stamp.touch()
            atomic_json(forged, forgery)
            try:
                verify_pilot_job(full/'request.json', cases[0]['id'], diagnostic, julia=str(fake))
            except ValueError:
                rejected = True
            # Ignore scheduler incomplete metadata deliberately: even an operator
            # override plus forged markers must reach and fail the final verifier.
            forged_command = [*command, '--ignore-incomplete']
            record['forged_command'] = forged_command
            with (folder/'workflow-forged.log').open('w', encoding='utf-8') as stream:
                replay = subprocess.run(forged_command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
            record['forged_workflow_returncode'] = replay.returncode
            if not rejected or replay.returncode == 0:
                raise ValueError('forged marker or ledger was accepted')
            if 'stale, missing or self-admitting cached pilot' not in (folder/'workflow-forged.log').read_text(encoding='utf-8'):
                raise ValueError('forged workflow did not reach the final cached verifier')
        finally:
            for stamp in stamps:
                stamp.unlink(missing_ok=True)
            forged.unlink(missing_ok=True)
        if fake.with_suffix('.calls').read_text().splitlines() != ['called']:
            raise ValueError('forged marker caused another numerical dispatch')
        if any(file_hash(p) != h for p, h in hashes.items()) or source_fingerprint(ROOT) != record['source_fingerprint']:
            raise ValueError('original failed evidence or source changed')
        record.update(status='SYNTHETIC_SCHEDULER_FAILURE_AND_FORGERY_CONTROLS_VERIFIED',
            raw_failure_hashes=hashes, original_failure_retained=True, audit_dispatched=False,
            forged_marker_rejected_before_dispatch=True,
            timeout_and_parent_signal_scope='Separate owned GNU-time-wrapped unit controls, not this Snakemake control')
    except BaseException as error:
        record['error'] = type(error).__name__+': '+str(error)
        raise
    finally:
        record['completed_at'] = utc_now()
        atomic_json(folder/'control.json', record)
    print(json.dumps({k: record[k] for k in ('status', 'scope', 'numerical_predictions_executed')}))


if __name__ == '__main__':
    main()
