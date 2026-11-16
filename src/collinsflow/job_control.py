"""Supervise complete numerical worker processes, outside parameter evaluation."""
from __future__ import annotations

import math
import os
import signal
import subprocess
import threading
import time
from typing import NamedTuple


class WholeJobResult(NamedTuple):
    returncode: int
    process_returncode: int
    interruption: str | None
    termination_signal: int | None


def run_complete_process(command, *, cwd, env, stdout, timeout_seconds=None):
    """Retain a nonzero result on interruption, including a child that exits0."""
    if timeout_seconds is not None and (isinstance(timeout_seconds, bool) or
            not math.isfinite(timeout_seconds) or timeout_seconds <= 0):
        raise ValueError("a positive finite complete-job time budget is required")
    received = []
    handlers = {}
    if threading.current_thread() is threading.main_thread():
        for sig in (signal.SIGINT, signal.SIGTERM):
            handlers[sig] = signal.getsignal(sig)
            signal.signal(sig, lambda number, frame: received.append(number))
    process = None
    started = time.monotonic()
    stop_started = None
    sent = 0
    interruption = None
    kill_signal = getattr(signal, "SIGKILL", signal.SIGTERM)

    def send(sig):
        if process is None or process.poll() is not None:
            return
        try:
            if os.name == "posix":
                os.killpg(process.pid, sig)
            elif sig == kill_signal:
                process.kill()
            else:
                process.terminate()
        except ProcessLookupError:
            pass

    try:
        process = subprocess.Popen(command, cwd=cwd, env=env, stdout=stdout,
            stderr=subprocess.STDOUT, start_new_session=os.name == "posix")
        while True:
            now = time.monotonic()
            if received and interruption is None:
                interruption = "parent_signal"
            if timeout_seconds is not None and now-started >= timeout_seconds and interruption is None:
                interruption = "resource_timeout"
            if interruption is not None and stop_started is None:
                stop_started = now
                send(signal.SIGINT)
                sent = 1
            if stop_started is not None and now-stop_started >= 2 and sent == 1:
                send(signal.SIGTERM)
                sent = 2
            if stop_started is not None and now-stop_started >= 3 and sent == 2:
                send(kill_signal)
                sent = 3
            try:
                code = process.wait(timeout=.05)
                if received and interruption is None:
                    interruption = "parent_signal"
                result_code = 128+received[0] if received else 124 if interruption else code
                return WholeJobResult(result_code, code, interruption, received[0] if received else None)
            except subprocess.TimeoutExpired:
                pass
    except BaseException:
        send(signal.SIGTERM)
        if process is not None:
            try:
                process.wait(timeout=1)
            except subprocess.TimeoutExpired:
                send(kill_signal)
                process.wait()
        raise
    finally:
        for sig, handler in handlers.items():
            signal.signal(sig, handler)
