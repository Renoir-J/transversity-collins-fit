"""Synthetic job/resource controls. No physical source loader or Julia prediction runs."""
import copy
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest

from collinsflow import historical_sidis_pilot as pilot
from collinsflow.provenance import file_hash, read_json, source_fingerprint

ROOT = Path(__file__).resolve().parents[1]
LINUX_TIME = os.name == 'posix' and Path('/usr/bin/time').is_file()


def invented_request():
    case = pilot.load_policy(ROOT)['selection']['cases'][0]
    points = [dict(id=case['point_id'], publication=case['publication'], target=case['target'])]
    points += [dict(id=f'explicitly-synthetic:{i}', publication='invented', target='invented') for i in range(131)]
    return dict(scope='explicitly_synthetic_test_only_not_source_replayed',
        source_preparation=dict(manifest=str(ROOT/'config/historical-sources.json')),
        source_fingerprint=source_fingerprint(ROOT), points=points)


def fake_worker(tmp_path, *, cooperative=False):
    path = tmp_path/'synthetic-worker'
    path.write_text(f'#!{sys.executable}\n'+'''import json,signal,sys,time
from pathlib import Path
out=Path(sys.argv[-1]);out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(dict(status="FAIL",scope="explicitly_synthetic_worker_control",production_ready=False)))
Path(str(out)+".progress.json").write_text(json.dumps(dict(status="IN_PROGRESS_NOT_A_GATE",scope="explicitly_synthetic_worker_control")))
print("explicitly synthetic partial worker",flush=True)
'''+('''
def stop(number,frame):
    print("explicitly synthetic cooperative exit0",flush=True)
    raise SystemExit(0)
signal.signal(signal.SIGINT,stop)
(out.parent/"ready").write_text("synthetic worker ready")
time.sleep(30)
''' if cooperative else 'raise SystemExit(7)\n'), encoding='utf-8')
    path.chmod(0o755)
    return path


@pytest.fixture
def fake_source_job(tmp_path, monkeypatch):
    # Only the unit supervisor receives this invented source representation.
    # The actual raw-source/Julia loaders reject it and are never replaced in production.
    request = invented_request()
    requestfile = tmp_path/'invented-request.json'
    requestfile.write_text(json.dumps(request), encoding='utf-8')
    native, parton = tmp_path/'invented-native.json', tmp_path/'invented-parton.json'
    native.write_text('{}')
    parton.write_text('{}')
    monkeypatch.setattr(pilot, 'verify_request', lambda path: copy.deepcopy(request))
    monkeypatch.setattr(pilot, 'available_memory_bytes', lambda: 8*1024**3)
    return requestfile, native, parton


def test_fixed_policy_keeps_partial_measurements_outside_every_admission():
    policy = pilot.load_policy(ROOT)
    assert all(policy[k] is False for k in pilot.FALSE_FLAGS)
    assert [c['point_id'] for c in policy['selection']['cases']] == [
        'KPSY15_SIDIS_3000:r2', 'KPSY15_SIDIS_3025:r2', 'KPSY15_SIDIS_4000:r2']
    assert policy['numerics']['native_level'] == 4
    assert policy['numerics']['derivative_evaluations']['derivative_accuracy_accepted'] is False


def test_resource_units_and_measurement_scope_are_explicit(tmp_path):
    path = tmp_path/'resource.txt'
    path.write_text('COLLINSFLOW_TIME_V1\t65536\t1.25\t0.8\t0.1\t0\n')
    actual = pilot.parse_resource_file(path)
    assert actual['maximum_RSS_KiB'] == 65536
    assert actual['wall_seconds'] == 1.25 and actual['child_exit_code'] == 0
    assert 'not_simultaneous_process_tree_peak' in actual['scope']


@pytest.mark.skipif(not LINUX_TIME, reason='Linux GNU time required')
def test_measured_RSS_has_KiB_scale_for_an_owned_known_allocation(tmp_path, monkeypatch):
    monkeypatch.setattr(pilot, 'available_memory_bytes', lambda: 8*1024**3)
    command = [sys.executable, '-c', 'b=bytearray(64*1024*1024);b[::4096]=b"x"*(len(b)//4096);print(len(b))']
    row = pilot._timed_job('known-allocation', command, ROOT, tmp_path, Path('/usr/bin/time'),
                           pilot.load_policy(ROOT), time.monotonic()+15)
    pilot._require_command_success(row, pilot.load_policy(ROOT))
    assert 64*1024 <= row['resource']['maximum_RSS_KiB'] < 300*1024


@pytest.mark.parametrize('line', ['', 'COLLINSFLOW_TIME_V1\t0\t1\t1\t1\t0',
    'COLLINSFLOW_TIME_V1\t1\tnan\t1\t1\t0', 'COLLINSFLOW_TIME_V1\t1\t-1\t1\t1\t0',
    'COLLINSFLOW_TIME_V1\t1\t1\t1\t1\t-1'])
def test_missing_or_invalid_resource_measurement_is_rejected(tmp_path, line):
    path = tmp_path/'resource.txt'
    path.write_text(line)
    with pytest.raises(ValueError):
        pilot.parse_resource_file(path)


def test_absent_and_corrupt_primary_files_fail_without_fallback(tmp_path):
    manifest = ROOT/'config/historical-sidis-primary-sources.json'
    cache = tmp_path/'cache'
    with pytest.raises(FileNotFoundError):
        pilot.prepare_primary_cache(manifest, cache, tmp_path/'missing.json', offline=True)
    assert read_json(tmp_path/'missing.json')['status'] == 'FAIL'
    source = read_json(manifest)['sources'][0]
    wrong = cache/source['filename']
    wrong.write_bytes(b'explicitly synthetic wrong bytes')
    before = file_hash(wrong)
    with pytest.raises(ValueError, match='checksum'):
        pilot.prepare_primary_cache(manifest, cache, tmp_path/'corrupt.json', offline=True)
    assert read_json(tmp_path/'corrupt.json')['status'] == 'FAIL' and file_hash(wrong) == before


@pytest.mark.skipif(not LINUX_TIME, reason='Linux GNU time and process groups required')
def test_failed_worker_preserves_logs_and_a_forged_marker_cannot_approve_it(tmp_path, monkeypatch, fake_source_job):
    request, native, parton = fake_source_job
    worker = fake_worker(tmp_path)
    out = tmp_path/'job'
    with pytest.raises(ValueError, match='failed complete pilot'):
        pilot.run_pilot_job(request, 'hermes-first', native, parton, out, julia=str(worker))
    invocation = read_json(out/'invocation.json')
    assert invocation['status'] == 'FAIL' and invocation['steps'][0]['returncode'] == 7
    assert invocation['steps'][0]['resource']['child_exit_code'] == 7
    assert len(invocation['steps']) == 1 and not (out/'audit.json').exists()
    assert read_json(out/'pilot.json')['status'] == 'FAIL'
    assert file_hash(out/'initial-invocation.json') == invocation['initial_record_sha256']
    (out/'job.completed').write_text('explicitly synthetic forged marker')
    monkeypatch.setattr(pilot, '_timed_job', lambda *a, **k: pytest.fail('forged marker must not dispatch a worker'))
    with pytest.raises(ValueError, match='cached pilot'):
        pilot.verify_pilot_job(request, 'hermes-first', out)


@pytest.mark.skipif(not LINUX_TIME, reason='Linux GNU time and process groups required')
def test_timeout_retains_failure_even_when_timed_child_and_wrapper_exit_zero(tmp_path, monkeypatch, fake_source_job):
    request, native, parton = fake_source_job
    frozen = pilot.load_policy(ROOT)
    synthetic = copy.deepcopy(frozen)
    synthetic['resources']['whole_Julia_jobs_time_budget_seconds'] = 2.
    monkeypatch.setattr(pilot, 'load_policy', lambda root: synthetic)
    worker = fake_worker(tmp_path, cooperative=True)
    out = tmp_path/'timeout'
    with pytest.raises(ValueError, match='interrupted'):
        pilot.run_pilot_job(request, 'hermes-first', native, parton, out, julia=str(worker))
    inv = read_json(out/'invocation.json')
    step = inv['steps'][0]
    assert inv['status'] == 'FAIL' and step['returncode'] == 124
    assert step['process_returncode'] == step['resource']['child_exit_code'] == 0
    assert step['interruption'] == 'resource_timeout'
    assert read_json(out/'pilot.json')['status'] == 'FAIL' and not (out/'audit.json').exists()
    assert 'cooperative exit0' in (out/'producer.log').read_text()


@pytest.mark.skipif(not LINUX_TIME, reason='Linux GNU time and process groups required')
def test_parent_signal_preserves_nonzero_status_and_cooperative_child_evidence(tmp_path):
    request = tmp_path/'invented-request.json'
    request.write_text(json.dumps(invented_request()))
    native, parton = tmp_path/'invented-native.json', tmp_path/'invented-parton.json'
    native.write_text('{}')
    parton.write_text('{}')
    worker = fake_worker(tmp_path, cooperative=True)
    out = tmp_path/'signal'
    harness = tmp_path/'harness.py'
    harness.write_text('''import json,sys
from pathlib import Path
from collinsflow import historical_sidis_pilot as p
request,native,parton,out,worker=sys.argv[1:]
p.verify_request=lambda path: json.loads(Path(request).read_text())
p.available_memory_bytes=lambda: 8*1024**3
try:p.run_pilot_job(request,"hermes-first",native,parton,out,julia=worker)
except ValueError:raise SystemExit(2)
''')
    with (tmp_path/'harness.log').open('w') as stream:
        process = subprocess.Popen([sys.executable, str(harness), str(request), str(native), str(parton), str(out), str(worker)],
            cwd=ROOT, env=dict(os.environ), stdout=stream, stderr=subprocess.STDOUT)
        try:
            deadline = time.monotonic()+15
            while not (out/'ready').exists() and process.poll() is None and time.monotonic()<deadline:
                time.sleep(.05)
            assert (out/'ready').exists(), (tmp_path/'harness.log').read_text()
            os.kill(process.pid, signal.SIGTERM)
            assert process.wait(timeout=8) != 0
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait(timeout=8)
    inv = read_json(out/'invocation.json')
    step = inv['steps'][0]
    assert inv['status'] == 'FAIL' and step['returncode'] == 143 and step['process_returncode'] == 0
    assert step['interruption'] == 'parent_signal' and step['termination_signal'] == signal.SIGTERM
    assert step['resource']['child_exit_code'] == 0


def test_insufficient_memory_never_dispatches_a_command(tmp_path, monkeypatch):
    monkeypatch.setattr(pilot, 'available_memory_bytes', lambda: 1)
    monkeypatch.setattr(pilot, 'run_complete_process', lambda *a, **k: pytest.fail('must not start'))
    with pytest.raises(ValueError, match='insufficient available memory'):
        pilot._timed_job('synthetic', ['not-executed'], ROOT, tmp_path, Path('/usr/bin/time'),
                         pilot.load_policy(ROOT), time.monotonic()+1)


@pytest.mark.skipif(os.name != 'posix', reason='Linux advisory lock')
def test_only_one_pilot_can_own_a_checkout():
    with pilot._pilot_lock(ROOT), pytest.raises(ValueError, match='another complete'):
        with pilot._pilot_lock(ROOT):
            pytest.fail('second pilot must not acquire lock')


@pytest.mark.parametrize('damage', ['missing_resource', 'rss_budget', 'signal'])
def test_zero_process_exit_does_not_override_missing_or_failed_resource_evidence(damage):
    policy = pilot.load_policy(ROOT)
    row = dict(name='synthetic', returncode=0, process_returncode=0, interruption=None, termination_signal=None,
               resource=dict(child_exit_code=0, maximum_RSS_KiB=1))
    if damage == 'missing_resource':
        del row['resource']
    elif damage == 'rss_budget':
        row['resource']['maximum_RSS_KiB'] = policy['resources']['maximum_command_RSS_KiB']+1
    else:
        row['termination_signal'] = signal.SIGTERM
    with pytest.raises(ValueError):
        pilot._require_command_success(row, policy)


def test_forged_ledger_and_all_markers_cannot_replace_absent_job_evidence(tmp_path, monkeypatch, fake_source_job):
    request, _, _ = fake_source_job
    run = tmp_path/'forged-run'
    for case in pilot.load_policy(ROOT)['selection']['cases']:
        directory = run/'cases'/case['id']
        directory.mkdir(parents=True)
        (directory/'job.completed').write_text('synthetic forged marker')
    ledger = run/'pilot-ledger.json'
    ledger.write_text(json.dumps(dict(status=pilot.COMPLETED, scope=pilot.SCOPE)))
    monkeypatch.setattr(pilot, '_timed_job', lambda *a, **k: pytest.fail('absent evidence must not dispatch'))
    with pytest.raises(FileNotFoundError):
        pilot.verify_pilot_ledger(request, run, ledger)
