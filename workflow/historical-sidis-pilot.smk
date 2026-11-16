"""Three resource measurements; full132 physics/fitting gates remain closed."""
from pathlib import Path
import shlex

from collinsflow.historical_sidis_pilot import build_pilot_ledger, load_policy, verify_pilot_ledger
from collinsflow.provenance import computational_sources

ROOT = Path(workflow.basedir).parent.resolve()
def absolute(value):
    return str((ROOT / str(value)).resolve())

RUN = absolute(config.get("run_dir", "runs/historical-sidis-pilot"))
CACHE = absolute(config.get("cache", "data/raw/historical"))
PRIMARY_CACHE = absolute(config.get("primary_cache", "data/raw/historical-sidis-primary"))
CANDIDATES = f"{RUN}/historical"
REQUEST = f"{RUN}/request.json"
PRIMARY_RECEIPT = f"{RUN}/primary.json"
LEDGER = f"{RUN}/pilot-ledger.json"
NATIVE = absolute(config.get("native_build", ".runtime/hoppet/build.json"))
PARTON = absolute(config.get("parton_build", ".runtime/parton-inputs/build.json"))
JULIA = config.get("julia", "julia")
TIME_BINARY = config.get("time_binary", "/usr/bin/time")
POLICY = load_policy(ROOT)
CASES = [c["id"] for c in POLICY["selection"]["cases"]]
CASE_INDEX = {name: i for i, name in enumerate(CASES)}
COMPLETIONS = [f"{RUN}/cases/{case}/job.completed" for case in CASES]
CODE = [str(path) for path in computational_sources(ROOT)]
SEED = ("--seed-root " + shlex.quote(absolute(config["seed_root"]))) if config.get("seed_root") else ""
OFFLINE = "--offline" if config.get("offline", False) else ""
HOPPET_SOURCE = ("--source " + shlex.quote(absolute(config["hoppet_source"]))) if config.get("hoppet_source") else ""

rule all:
    input:
        LEDGER
    run:
        verify_pilot_ledger(REQUEST, RUN, LEDGER, julia=JULIA, time_binary=TIME_BINARY)

rule historical_candidates:
    input:
        manifest=str(ROOT / "config/historical-sources.json"), code=CODE
    output:
        candidates=f"{CANDIDATES}/candidates.json", report=f"{CANDIDATES}/report.json"
    params:
        root=str(ROOT), cache=CACHE, output=CANDIDATES, seed=SEED, offline=OFFLINE
    threads: 1
    resources:
        mem_mb=512
    shell:
        "cd {params.root:q} && python -m collinsflow.cli historical-intake --manifest {input.manifest:q} "
        "--cache {params.cache:q} --output {params.output:q} {params.seed} {params.offline}"

rule primary_convention_bytes:
    input:
        manifest=str(ROOT / "config/historical-sidis-primary-sources.json"), code=CODE
    output:
        PRIMARY_RECEIPT
    params:
        root=str(ROOT), cache=PRIMARY_CACHE, offline=OFFLINE
    threads: 1
    resources:
        mem_mb=256
    shell:
        "cd {params.root:q} && python -m collinsflow.cli historical-sidis-primary "
        "--manifest {input.manifest:q} --cache {params.cache:q} --output {output:q} {params.offline}"

rule source_request:
    input:
        candidates=f"{CANDIDATES}/candidates.json", report=f"{CANDIDATES}/report.json",
        primary=PRIMARY_RECEIPT, code=CODE
    output:
        REQUEST
    params:
        root=str(ROOT), manifest=str(ROOT / "config/historical-sources.json"), cache=CACHE,
        candidates=CANDIDATES, primary_cache=PRIMARY_CACHE
    threads: 1
    resources:
        mem_mb=512
    shell:
        "cd {params.root:q} && python -m collinsflow.cli historical-sidis-point-request --manifest {params.manifest:q} "
        "--cache {params.cache:q} --candidates {params.candidates:q} --primary-cache {params.primary_cache:q} --output {output:q}"

if not config.get("native_build"):
    rule retained_native_build:
        input:
            bridge=str(ROOT / "fortran/hoppet_bridge.f90"), contract=str(ROOT / "config/hoppet-contract.json"),
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
            bridge=str(ROOT / "bridges/lhapdf.cpp"), sources=str(ROOT / "config/parton-sources.json"),
            contract=str(ROOT / "config/parton-contract.json"), prepare=str(ROOT / "src/collinsflow/parton.py")
        output:
            PARTON
        params:
            root=str(ROOT), destination=str(Path(PARTON).parent), offline=OFFLINE
        threads: 2
        resources:
            mem_mb=2048
        shell:
            "cd {params.root:q} && python scripts/prepare_parton_inputs.py --destination {params.destination:q} {params.offline} --jobs {threads}"

wildcard_constraints:
    case="|".join(CASES)

rule measured_pilot_case:
    input:
        request=REQUEST, native=NATIVE, parton=PARTON, code=CODE,
        previous=lambda w: [] if CASE_INDEX[w.case] == 0 else [COMPLETIONS[CASE_INDEX[w.case]-1]]
    output:
        completion=touch(f"{RUN}/cases/{{case}}/job.completed")
    params:
        root=str(ROOT), output=lambda w: f"{RUN}/cases/{w.case}", julia=JULIA, time_binary=TIME_BINARY
    threads: 1
    resources:
        mem_mb=3072
    shell:
        "cd {params.root:q} && python -m collinsflow.cli historical-sidis-pilot --request {input.request:q} --case {wildcards.case:q} "
        "--native-build {input.native:q} --parton-build {input.parton:q} --output {params.output:q} "
        "--julia {params.julia:q} --time-binary {params.time_binary:q}"

rule measurement_ledger:
    input:
        completed=COMPLETIONS, request=REQUEST, code=CODE
    output:
        LEDGER
    run:
        build_pilot_ledger(REQUEST, RUN, str(output[0]), julia=JULIA, time_binary=TIME_BINARY)
