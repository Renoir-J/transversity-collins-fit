"""Real-row fixed-parameter diagnostics; all fitting/production gates remain closed."""
from pathlib import Path
import shlex

from collinsflow.historical_sia_points import verify_point_job
from collinsflow.provenance import computational_sources

ROOT = Path(workflow.basedir).parent
RUN = config.get("run_dir", "runs/historical-sia-points")
MANIFEST = str(ROOT / config.get("manifest", "config/historical-sources.json"))
CACHE = config.get("cache", "data/raw/historical")
CANDIDATES = f"{RUN}/historical"
REQUEST = f"{RUN}/request.json"
NUMERICAL = f"{RUN}/diagnostic"
COMPLETED = f"{NUMERICAL}/job.completed"
NATIVE = str(ROOT / config.get("native_build", ".runtime/hoppet/build.json"))
PARTON = str(ROOT / config.get("parton_build", ".runtime/parton-inputs/build.json"))
JULIA = config.get("julia", "julia")
SEED = ("--seed-root " + shlex.quote(str(config["seed_root"]))) if config.get("seed_root") else ""
OFFLINE = "--offline" if config.get("offline", False) else ""
HOPPET_SOURCE = ("--source " + shlex.quote(str(config["hoppet_source"]))) if config.get("hoppet_source") else ""
CODE = [str(path) for path in computational_sources(ROOT)]

rule all:
    input:
        COMPLETED
    run:
        verify_point_job(REQUEST, NUMERICAL, julia=JULIA)

rule historical_candidates:
    input:
        manifest=MANIFEST,
        code=CODE
    output:
        candidates=f"{CANDIDATES}/candidates.json",
        report=f"{CANDIDATES}/report.json"
    params:
        root=str(ROOT), output=CANDIDATES, cache=CACHE, seed=SEED, offline=OFFLINE
    threads: 1
    resources:
        mem_mb=512
    shell:
        "cd {params.root:q} && python -m collinsflow.cli historical-intake --manifest {input.manifest:q} "
        "--cache {params.cache:q} --output {params.output:q} {params.seed} {params.offline}"

rule source_request:
    input:
        candidates=f"{CANDIDATES}/candidates.json",
        report=f"{CANDIDATES}/report.json",
        code=CODE
    output:
        REQUEST
    params:
        root=str(ROOT), manifest=MANIFEST, cache=CACHE, candidates=CANDIDATES
    threads: 1
    resources:
        mem_mb=512
    shell:
        "cd {params.root:q} && python -m collinsflow.cli historical-sia-point-request --manifest {params.manifest:q} "
        "--cache {params.cache:q} --candidates {params.candidates:q} --output {output:q}"

if not config.get("native_build"):
    rule retained_native_build:
        input:
            bridge=str(ROOT / "fortran/hoppet_bridge.f90"),
            contract=str(ROOT / "config/hoppet-contract.json"),
            prepare=str(ROOT / "src/collinsflow/hoppet.py")
        output:
            NATIVE
        params:
            root=str(ROOT), destination=str(Path(NATIVE).parent), source=HOPPET_SOURCE, offline=OFFLINE
        threads: 2
        resources:
            mem_mb=2048
        shell:
            "cd {params.root:q} && python scripts/prepare_hoppet.py --destination {params.destination:q} "
            "{params.source} {params.offline} --jobs {threads}"

if not config.get("parton_build"):
    rule pinned_parton_build:
        input:
            bridge=str(ROOT / "bridges/lhapdf.cpp"),
            sources=str(ROOT / "config/parton-sources.json"),
            contract=str(ROOT / "config/parton-contract.json"),
            prepare=str(ROOT / "src/collinsflow/parton.py")
        output:
            PARTON
        params:
            root=str(ROOT), destination=str(Path(PARTON).parent), offline=OFFLINE
        threads: 2
        resources:
            mem_mb=2048
        shell:
            "cd {params.root:q} && python scripts/prepare_parton_inputs.py --destination {params.destination:q} "
            "{params.offline} --jobs {threads}"

rule historical_sia_point_diagnostic:
    input:
        request=REQUEST,
        native=NATIVE,
        parton=PARTON,
        code=CODE
    output:
        completion=touch(COMPLETED)
    params:
        root=str(ROOT), output=NUMERICAL, julia=JULIA
    threads: 1
    resources:
        mem_mb=6144
    shell:
        "cd {params.root:q} && python -m collinsflow.cli historical-sia-point-predictions --request {input.request:q} "
        "--native-build {input.native:q} --parton-build {input.parton:q} --output {params.output:q} --julia {params.julia:q}"
