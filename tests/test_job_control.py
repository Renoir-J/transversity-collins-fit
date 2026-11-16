"""Owned synthetic subprocesses; no numerical worker or production admission."""
import json
import os
import signal
import subprocess
import sys
import time

import pytest

from collinsflow.job_control import run_complete_process


@pytest.mark.parametrize('code', [0, 7])
def test_retains_worker_exit_and_log(tmp_path, code):
    with (tmp_path/'worker.log').open('w') as log:
        result = run_complete_process([sys.executable, '-c', f'print("synthetic log"); raise SystemExit({code})'],
            cwd=tmp_path, env=dict(os.environ), stdout=log)
    assert result.returncode == result.process_returncode == code
    assert result.interruption is result.termination_signal is None
    assert 'synthetic log' in (tmp_path/'worker.log').read_text()


@pytest.mark.parametrize('budget', [True, 0, -1, float('nan'), float('inf')])
def test_invalid_budget_never_starts_a_worker(tmp_path, monkeypatch, budget):
    def unexpected(*a, **k):
        raise AssertionError('worker must not start')
    monkeypatch.setattr(subprocess, 'Popen', unexpected)
    with (tmp_path/'worker.log').open('w') as log, pytest.raises(ValueError):
        run_complete_process(['not-started'], cwd=tmp_path, env=dict(os.environ), stdout=log, timeout_seconds=budget)


COOPERATIVE = '''import signal,time
def stop(number,frame):
    print("explicitly synthetic cooperative stop",flush=True)
    raise SystemExit(0)
signal.signal(signal.SIGINT,stop)
print("explicitly synthetic partial progress",flush=True)
time.sleep(10)
'''


@pytest.mark.skipif(os.name != 'posix', reason='Linux process-group interruption contract')
def test_timeout_never_accepts_cooperative_exit_zero(tmp_path):
    before = {sig: signal.getsignal(sig) for sig in (signal.SIGINT, signal.SIGTERM)}
    with (tmp_path/'worker.log').open('w') as log:
        result = run_complete_process([sys.executable, '-c', COOPERATIVE], cwd=tmp_path,
            env=dict(os.environ), stdout=log, timeout_seconds=2.)
    assert result.returncode == 124 and result.process_returncode == 0
    assert result.interruption == 'resource_timeout'
    assert 'partial progress' in (tmp_path/'worker.log').read_text()
    assert 'cooperative stop' in (tmp_path/'worker.log').read_text()
    assert all(signal.getsignal(sig) == old for sig, old in before.items())


@pytest.mark.skipif(os.name != 'posix', reason='Linux process-group interruption contract')
def test_parent_SIGTERM_is_not_a_successful_job(tmp_path):
    child = tmp_path/'child.py'
    child.write_text(COOPERATIVE.replace('print("explicitly synthetic partial progress",flush=True)',
        'from pathlib import Path\nPath("ready").write_text("synthetic child ready")\nprint("explicitly synthetic partial progress",flush=True)'))
    harness = tmp_path/'harness.py'
    harness.write_text('''import json,os,sys
from pathlib import Path
from collinsflow.job_control import run_complete_process
with open("worker.log","w") as log:
    r=run_complete_process([sys.executable,"child.py"],cwd=Path.cwd(),env=dict(os.environ),stdout=log)
Path("outcome.json").write_text(json.dumps(r._asdict()))
raise SystemExit(r.returncode)
''')
    with (tmp_path/'harness.log').open('w') as log:
        parent = subprocess.Popen([sys.executable, str(harness)], cwd=tmp_path, env=dict(os.environ), stdout=log, stderr=subprocess.STDOUT)
        try:
            deadline = time.monotonic()+10
            while not (tmp_path/'ready').exists() and parent.poll() is None and time.monotonic() < deadline:
                time.sleep(.05)
            assert (tmp_path/'ready').exists(), (tmp_path/'harness.log').read_text()
            os.kill(parent.pid, signal.SIGTERM)
            assert parent.wait(timeout=8) == 128+signal.SIGTERM
        finally:
            if parent.poll() is None:
                parent.terminate()
                parent.wait(timeout=8)
    result = json.loads((tmp_path/'outcome.json').read_text())
    assert result['returncode'] == 128+signal.SIGTERM and result['process_returncode'] == 0
    assert result['interruption'] == 'parent_signal' and result['termination_signal'] == signal.SIGTERM
    assert 'partial progress' in (tmp_path/'worker.log').read_text()
