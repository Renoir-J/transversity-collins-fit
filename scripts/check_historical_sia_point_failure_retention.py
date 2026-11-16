"""Verify raw failure retention with an explicitly synthetic failing whole job."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from collinsflow.historical_sia_points import verify_point_job
from collinsflow.provenance import atomic_json, file_hash, git_commit, read_json, source_fingerprint, utc_now

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', required=True)
    parser.add_argument('--cache', default='data/raw/historical')
    parser.add_argument('--native-build', default='.runtime/hoppet/build.json')
    parser.add_argument('--parton-build', default='.runtime/parton-inputs/build.json')
    parser.add_argument('--seed-root')
    args = parser.parse_args()
    if os.name != 'posix':
        raise ValueError('run the scheduler retention control on Linux/WSL')
    folder = (ROOT/args.run_dir).resolve()
    if not folder.is_relative_to(ROOT/'runs'):
        raise ValueError('use a fresh ignored runs directory for this failure control')
    folder.mkdir(parents=True, exist_ok=False)
    full = folder/'full'
    runtime = ROOT/'.runtime'
    runtime.mkdir(exist_ok=True)
    fake_dir = Path(tempfile.mkdtemp(prefix='historical-sia-failure-control-', dir=runtime))
    fake = fake_dir/'intentional-failing-julia'
    payload = dict(status='FAIL', scope='explicitly_synthetic_failure_retention_control',
        production_ready=False, physical_extraction=False, experimental_data_admitted=False,
        fit_performed=False, numerical_predictions_executed=False,
        reason='Intentional complete-job failure before Julia or numerical predictions')
    fake.write_text(f'''#!{sys.executable}
from pathlib import Path
import json
import sys
payload = {payload!r}
Path(sys.argv[-1]).write_text(json.dumps(payload, indent=2)+'\\n', encoding='utf-8')
with Path(__file__).with_suffix('.calls').open('a', encoding='utf-8') as stream:
    stream.write('called\\n')
raise SystemExit(2)
''', encoding='utf-8')
    fake.chmod(0o700)
    config = ['run_dir='+str(full), 'cache='+str((ROOT/args.cache).resolve()),
        'native_build='+str((ROOT/args.native_build).resolve()),
        'parton_build='+str((ROOT/args.parton_build).resolve()), 'julia='+str(fake), 'offline=True']
    if args.seed_root:
        config.append('seed_root='+str(Path(args.seed_root).resolve()))
    command = ['snakemake', '--snakefile', str(ROOT/'workflow/historical-sia-points.smk'),
               '--cores', '1', '--config', *config]
    record = dict(status='FAIL', scope='explicitly_synthetic_failure_retention_control',
        production_ready=False, physical_extraction=False, experimental_data_admitted=False,
        accepted_prediction_diagnostic=False, accepted_fitting_gate=False,
        numerical_predictions_executed=False, started_at=utc_now(),
        source_fingerprint=source_fingerprint(ROOT), code_commit=git_commit(ROOT),
        control_script_sha256=file_hash(__file__), command=command,
        fake_executable_sha256=file_hash(fake))
    try:
        with (folder/'workflow.log').open('w', encoding='utf-8') as stream:
            result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
        record['workflow_returncode'] = result.returncode
        if result.returncode == 0:
            raise ValueError('deliberate failed whole job unexpectedly succeeded')
        diagnostic = full/'diagnostic'
        report, invocation = diagnostic/'predictions.json', diagnostic/'invocation.json'
        if read_json(report) != payload:
            raise ValueError('synthetic failed report was lost or changed by scheduler cleanup')
        inv = read_json(invocation)
        if (inv['status'] != 'FAIL' or 'completed_at' not in inv or
                [(step['name'], step['returncode']) for step in inv['steps']] != [('predictions', 2)]):
            raise ValueError('expected explicit failed invocation was lost or changed')
        stamp = diagnostic/'job.completed'
        if stamp.exists():
            raise ValueError('failed job created a completion marker')
        hashes = {str(path): file_hash(path) for path in (report, invocation)}
        rejected = False
        stamp.touch()
        try:
            try:
                verify_point_job(full/'request.json', diagnostic, julia=str(fake))
            except ValueError:
                rejected = True
        finally:
            stamp.unlink()
        if not rejected or fake.with_suffix('.calls').read_text(encoding='utf-8').splitlines() != ['called']:
            raise ValueError('cached verifier failed to reject a forged scheduling marker before numerical dispatch')
        if any(file_hash(path) != sha for path, sha in hashes.items()):
            raise ValueError('failed original artifacts changed during retention control')
        if source_fingerprint(ROOT) != record['source_fingerprint']:
            raise ValueError('source changed during failure-retention control')
        record.update(status='FAILURE_RETENTION_CONTROL_VERIFIED', raw_failure_hashes=hashes,
            original_failed_report_and_invocation_retained=True, completion_marker_produced=False,
            direct_cached_verifier_rejected_forged_marker=True, numerical_audit_dispatched=False)
    except Exception as exc:
        record['error'] = type(exc).__name__+': '+str(exc)
        raise
    finally:
        record['completed_at'] = utc_now()
        atomic_json(folder/'control.json', record)
    print(json.dumps({key: record[key] for key in ('status', 'scope', 'numerical_predictions_executed')}))


if __name__ == '__main__':
    main()
