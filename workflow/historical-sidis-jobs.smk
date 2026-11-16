"""Complete132-point Julia validation; no fitting or production admission."""
from pathlib import Path
import shlex
from collinsflow.historical_sidis_jobs import load_policy, verify_complete_job, _selection
from collinsflow.provenance import computational_sources, read_json, source_fingerprint

ROOT = Path(workflow.basedir).parent.resolve()
def absolute(value):
    return str((ROOT/str(value)).resolve())
RUN = absolute(config.get("run_dir", "runs/historical-sidis-jobs"))
CACHE = absolute(config.get("cache", "data/raw/historical"))
PRIMARY_CACHE = absolute(config.get("primary_cache", "data/raw/historical-sidis-primary"))
CANDIDATES = f"{RUN}/historical"
REQUEST = f"{RUN}/request.json"
NATIVE = absolute(config.get("native_build", ".runtime/hoppet/build.json"))
PARTON = absolute(config.get("parton_build", ".runtime/parton-inputs/build.json"))
JULIA = config.get("julia", "julia")
TIMEOUT = "" if config.get("timeout_seconds") is None else " --timeout-seconds "+shlex.quote(str(config["timeout_seconds"]))
SEED = ("--seed-root "+shlex.quote(absolute(config["seed_root"]))) if config.get("seed_root") else ""
OFFLINE = "--offline" if config.get("offline", False) else ""
HOPPET_SOURCE = ("--source "+shlex.quote(absolute(config["hoppet_source"]))) if config.get("hoppet_source") else ""
POLICY = load_policy()
KEYS = [f"p{i:04d}" for i in range(1, 133)]
CODE = [str(p) for p in computational_sources(ROOT)]

def point_index(key):
    if key not in KEYS:
        raise ValueError("exact132 point index required")
    return KEYS.index(key)+1

def previous(phase, key):
    i=point_index(key)
    return [] if i==1 else [f"{RUN}/{phase}/{KEYS[i-2]}/job.completed"]

def native_final_inputs(wildcards):
    first=checkpoints.native_first_selection.get()
    decision=_selection(RUN,"native-select3",read_json(REQUEST),source_fingerprint(ROOT))
    return list(first.output) + ([] if decision["selected_level"]==3 else
        [f"{RUN}/native4/{key}/job.completed" for key in KEYS])

def final_reference_inputs(wildcards):
    first=checkpoints.derivative_first_selection.get()
    decision=_selection(RUN,"derivative-select3",read_json(REQUEST),source_fingerprint(ROOT))
    return list(first.output) + ([] if decision["selected_reference_level"]==3 else
        [f"{RUN}/reference4/{key}/job.completed" for key in KEYS])

COMMON = " --request "+shlex.quote(REQUEST)+" --native-build "+shlex.quote(NATIVE)+" --parton-build "+shlex.quote(PARTON)+" --run-dir "+shlex.quote(RUN)+" --julia "+shlex.quote(str(JULIA))+TIMEOUT

rule all:
    input:
        f"{RUN}/diagnostic/job.completed"
    run:
        verify_complete_job(REQUEST,NATIVE,PARTON,RUN,"diagnostic",julia=JULIA,timeout_seconds=config.get("timeout_seconds"))

rule historical_candidates:
    input:
        manifest=str(ROOT/"config/historical-sources.json"), code=CODE
    output:
        f"{CANDIDATES}/candidates.json", f"{CANDIDATES}/report.json"
    params:
        root=str(ROOT),cache=CACHE,output=CANDIDATES,seed=SEED,offline=OFFLINE
    shell:
        "cd {params.root:q} && python -m collinsflow.cli historical-intake --manifest {input.manifest:q} --cache {params.cache:q} --output {params.output:q} {params.seed} {params.offline}"

rule primary_conventions:
    input:
        manifest=str(ROOT/"config/historical-sidis-primary-sources.json"),code=CODE
    output:
        f"{RUN}/primary.json"
    params:
        root=str(ROOT),cache=PRIMARY_CACHE,offline=OFFLINE
    shell:
        "cd {params.root:q} && python -m collinsflow.cli historical-sidis-primary --manifest {input.manifest:q} --cache {params.cache:q} --output {output:q} {params.offline}"

rule source_request:
    input:
        candidates=f"{CANDIDATES}/candidates.json",report=f"{CANDIDATES}/report.json",primary=f"{RUN}/primary.json",code=CODE
    output:
        REQUEST
    params:
        root=str(ROOT),manifest=str(ROOT/"config/historical-sources.json"),cache=CACHE,candidates=CANDIDATES,primary=PRIMARY_CACHE
    shell:
        "cd {params.root:q} && python -m collinsflow.cli historical-sidis-point-request --manifest {params.manifest:q} --cache {params.cache:q} --candidates {params.candidates:q} --primary-cache {params.primary:q} --output {output:q}"

if not config.get("native_build"):
    rule native_build:
        input:
            bridge=str(ROOT/"fortran/hoppet_bridge.f90"),contract=str(ROOT/"config/hoppet-contract.json"),prepare=str(ROOT/"src/collinsflow/hoppet.py")
        output:
            NATIVE
        params:
            root=str(ROOT),destination=str(Path(NATIVE).parent),source=HOPPET_SOURCE,offline=OFFLINE
        threads: 2
        shell:
            "cd {params.root:q} && python scripts/prepare_hoppet.py --destination {params.destination:q} {params.source} {params.offline} --jobs {threads}"

if not config.get("parton_build"):
    rule parton_build:
        input:
            bridge=str(ROOT/"bridges/lhapdf.cpp"),sources=str(ROOT/"config/parton-sources.json"),contract=str(ROOT/"config/parton-contract.json"),prepare=str(ROOT/"src/collinsflow/parton.py")
        output:
            PARTON
        params:
            root=str(ROOT),destination=str(Path(PARTON).parent),offline=OFFLINE
        threads: 2
        shell:
            "cd {params.root:q} && python scripts/prepare_parton_inputs.py --destination {params.destination:q} {params.offline} --jobs {threads}"

wildcard_constraints:
    point="p[0-9]{4}"

rule native_first_points:
    input:
        request=REQUEST,native=NATIVE,parton=PARTON,previous=lambda w:previous("native123",w.point),code=CODE
    output:
        touch(f"{RUN}/native123/{{point}}/job.completed")
    params:
        root=str(ROOT),index=lambda w:point_index(w.point),common=COMMON
    shell:
        "cd {params.root:q} && python scripts/run_historical_sidis_job.py run --phase native123 --point-index {params.index} {params.common}"

checkpoint native_first_selection:
    input:
        completed=[f"{RUN}/native123/{key}/job.completed" for key in KEYS],request=REQUEST,native=NATIVE,parton=PARTON,code=CODE
    output:
        touch(f"{RUN}/native-select3/job.completed")
    params:
        root=str(ROOT),common=COMMON
    shell:
        "cd {params.root:q} && python scripts/run_historical_sidis_job.py run --phase native-select3 {params.common}"

rule native_fourth_points:
    input:
        selection=lambda w:checkpoints.native_first_selection.get().output,request=REQUEST,native=NATIVE,parton=PARTON,
        previous=lambda w:previous("native4",w.point),code=CODE
    output:
        touch(f"{RUN}/native4/{{point}}/job.completed")
    params:
        root=str(ROOT),index=lambda w:point_index(w.point),common=COMMON
    shell:
        "cd {params.root:q} && python scripts/run_historical_sidis_job.py run --phase native4 --point-index {params.index} {params.common}"

rule native_final_selection:
    input:
        completed=native_final_inputs,request=REQUEST,native=NATIVE,parton=PARTON,code=CODE
    output:
        touch(f"{RUN}/native-selection/job.completed")
    params:
        root=str(ROOT),common=COMMON
    shell:
        "cd {params.root:q} && python scripts/run_historical_sidis_job.py run --phase native-selection {params.common}"

rule source_math:
    input:
        selected=f"{RUN}/native-selection/job.completed",request=REQUEST,native=NATIVE,parton=PARTON,code=CODE
    output:
        touch(f"{RUN}/source-math/job.completed")
    params:
        root=str(ROOT),common=COMMON
    shell:
        "cd {params.root:q} && python scripts/run_historical_sidis_job.py run --phase source-math {params.common}"

rule first_reference_points:
    input:
        math=f"{RUN}/source-math/job.completed",request=REQUEST,native=NATIVE,parton=PARTON,
        previous=lambda w:previous("tailref123",w.point),code=CODE
    output:
        touch(f"{RUN}/tailref123/{{point}}/job.completed")
    params:
        root=str(ROOT),index=lambda w:point_index(w.point),common=COMMON
    shell:
        "cd {params.root:q} && python scripts/run_historical_sidis_job.py run --phase tailref123 --point-index {params.index} {params.common}"

checkpoint derivative_first_selection:
    input:
        completed=[f"{RUN}/tailref123/{key}/job.completed" for key in KEYS],request=REQUEST,native=NATIVE,parton=PARTON,code=CODE
    output:
        touch(f"{RUN}/derivative-select3/job.completed")
    params:
        root=str(ROOT),common=COMMON
    shell:
        "cd {params.root:q} && python scripts/run_historical_sidis_job.py run --phase derivative-select3 {params.common}"

rule fourth_reference_points:
    input:
        selection=lambda w:checkpoints.derivative_first_selection.get().output,request=REQUEST,native=NATIVE,parton=PARTON,
        previous=lambda w:previous("reference4",w.point),code=CODE
    output:
        touch(f"{RUN}/reference4/{{point}}/job.completed")
    params:
        root=str(ROOT),index=lambda w:point_index(w.point),common=COMMON
    shell:
        "cd {params.root:q} && python scripts/run_historical_sidis_job.py run --phase reference4 --point-index {params.index} {params.common}"

rule final_diagnostic:
    input:
        completed=final_reference_inputs,request=REQUEST,native=NATIVE,parton=PARTON,code=CODE
    output:
        touch(f"{RUN}/diagnostic/job.completed")
    params:
        root=str(ROOT),common=COMMON
    shell:
        "cd {params.root:q} && python scripts/run_historical_sidis_job.py run --phase diagnostic {params.common}"
