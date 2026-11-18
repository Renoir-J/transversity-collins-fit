"""Real Linux descriptor observations on private test locks; no production lease."""
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from collinsflow import sia_geometry_host as host

ROOT=Path(__file__).resolve().parents[1]


@pytest.mark.skipif(not sys.platform.startswith("linux") or not shutil.which("julia"), reason="Linux and Julia required")
def test_julia_reads_real_private_flock_and_rejects_stale_or_substituted_proofs(tmp_path):
    private=tmp_path/"private-fixture.lock"
    other=tmp_path/"other-private.lock"
    # Neither helper acquires the standard SIDIS/SIA lock.
    assert private != host.standard_lock_path()
    with host._lock_at(private) as proof, host._lock_at(other) as other_proof:
        cases=[dict(name="positive-live-private-lock",proof=proof,expected_path=str(private))]
        for field,value in (("pid",True),("start_ticks",False),("start_ticks",proof["start_ticks"]+1),
            ("boot_id","11111111-2222-4333-8444-555555555555"),("file_descriptor",-1),
            ("inode",proof["inode"]+1),("device",proof["device"]+1),("lock_path",str(other)),
            ("file_descriptor",other_proof["file_descriptor"])):
            changed=copy.deepcopy(proof)
            changed[field]=value
            cases.append(dict(name=f"changed-{field}-{len(cases)}",proof=changed,expected_path=str(private)))
        cases.append(dict(name="changed-expected-path",proof=proof,expected_path=str(other)))
        # A live open descriptor with matching metadata is insufficient without
        # an exclusive flock. Observe it while it is genuinely open.
        unheld=tmp_path/"unheld.lock"
        descriptor=os.open(unheld,os.O_RDWR|os.O_CREAT,0o600)
        try:
            info=os.fstat(descriptor)
            no_lock={**proof,"file_descriptor":descriptor,"lock_path":str(unheld),"device":info.st_dev,"inode":info.st_ino}
            cases.append(dict(name="live-open-but-unlocked",proof=no_lock,expected_path=str(unheld)))
            fixture=tmp_path/"cases.json"
            fixture.write_text(json.dumps(dict(scope="explicitly_synthetic_private_resource_cases_not_dispatch",cases=cases)))
            output=tmp_path/"result.json"
            command=[shutil.which("julia"),"--startup-file=no","--project="+str(ROOT/"numerics"),
                str(ROOT/"numerics/test/sia_resource_lease_probe.jl"),str(fixture),str(output)]
            result=subprocess.run(command,cwd=ROOT,env={**os.environ,"JULIA_NUM_THREADS":"1","OPENBLAS_NUM_THREADS":"1","OMP_NUM_THREADS":"1"},
                capture_output=True,text=True,timeout=180)
            assert result.returncode==0,result.stdout+result.stderr
            report=json.loads(output.read_text())
            assert [row["accepted"] for row in report["cases"]]==[True]+[False]*(len(cases)-1),report
            assert report["cases"][0]["owner_pid"]==os.getpid()
            assert report["production_locks_acquired"] is False
            assert all(report[key] is False for key in ("physical_Julia_jobs_started","fit_performed","production_ready","accepted_fitting_gate"))
        finally:
            os.close(descriptor)
