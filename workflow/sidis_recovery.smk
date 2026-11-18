"""One serial complete-job driver; retries never overwrite interrupted outputs."""
from pathlib import Path

ADAPTER = str(Path(workflow.basedir).parent)
FROZEN = config["frozen_checkout"]
PREDECESSOR = config["predecessor_checkout"]
MANIFEST = config["recovery_manifest"]
MANIFEST_SHA = config["recovery_manifest_sha256"]
RUN = str(Path(MANIFEST).parent)

rule all:
    input:
        f"{RUN}/diagnostic/job.completed"
    # A cached marker is never admission. Always launch a fresh Julia reader.
    shell:
        "python {ADAPTER:q}/scripts/run_sidis_recovery.py verify --frozen-checkout {FROZEN:q} "
        "--predecessor-checkout {PREDECESSOR:q} --manifest {MANIFEST:q} --manifest-sha {MANIFEST_SHA:q}"

rule complete_recovery:
    input:
        MANIFEST
    output:
        f"{RUN}/diagnostic/job.completed"
    threads: 1
    resources:
        mem_mb=6144
    retries: 0
    shell:
        "python {ADAPTER:q}/scripts/run_sidis_recovery.py all --frozen-checkout {FROZEN:q} "
        "--predecessor-checkout {PREDECESSOR:q} --manifest {MANIFEST:q} --manifest-sha {MANIFEST_SHA:q}"
