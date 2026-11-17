"""Whole-Julia-job continuation; raw evidence is never a cleanup output."""
from pathlib import Path
from collinsflow.provenance import computational_sources

ROOT = Path(workflow.basedir).parent.resolve()
MANIFEST = str(Path(config["manifest"]).resolve())
FROZEN = str(Path(config["frozen_checkout"]).resolve())
RUN = str(Path(MANIFEST).parent)
JULIA = str(config.get("julia", "julia"))
CODE = [str(p) for p in computational_sources(ROOT)] + [str(ROOT/"scripts/run_sidis_continuation.py")]

rule all:
    input:
        f"{RUN}/continuation.completed"
    params:
        root=str(ROOT), manifest=MANIFEST, frozen=FROZEN, julia=JULIA
    shell:
        "cd {params.root:q} && python scripts/run_sidis_continuation.py verify --manifest {params.manifest:q} --frozen-checkout {params.frozen:q} --julia {params.julia:q} --phase diagnostic"

rule complete_continuation:
    input:
        manifest=MANIFEST, code=CODE
    output:
        touch(f"{RUN}/continuation.completed")
    params:
        root=str(ROOT), frozen=FROZEN, julia=JULIA
    shell:
        "cd {params.root:q} && python scripts/run_sidis_continuation.py all --manifest {input.manifest:q} --frozen-checkout {params.frozen:q} --julia {params.julia:q}"
