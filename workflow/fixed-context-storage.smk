"""Complete source-bound fixed-context storage jobs; no fitting admission."""
from pathlib import Path
import shlex
from collinsflow.provenance import computational_sources
from collinsflow.fixed_context_storage import build_or_verify_ledger

ROOT = Path(workflow.basedir).parent.resolve()
def absolute(value):
    return str((ROOT / str(value)).resolve())

RUN = absolute(config.get("run_dir", "runs/fixed-context-storage"))
CACHE = absolute(config.get("cache", "data/raw/historical"))
PRIMARY = absolute(config.get("primary_cache", "data/raw/historical-sidis-primary"))
NATIVE = absolute(config.get("native_build", ".runtime/hoppet/build.json"))
PARTON = absolute(config.get("parton_build", ".runtime/parton-inputs/build.json"))
JULIA = config.get("julia", "julia")
SIDIS = RUN + "/sidis-request.json"
SIA = RUN + "/sia-request.json"
CASES = [f"case{i:02d}" for i in range(1, 5)]
CODE = [str(p) for p in computational_sources(ROOT)]
OFFLINE = "--offline" if config.get("offline", False) else ""
SEED = "--seed-root " + shlex.quote(absolute(config["seed_root"])) if config.get("seed_root") else ""
TIMEOUT = "" if config.get("timeout_seconds") is None else " --timeout-seconds " + shlex.quote(str(config["timeout_seconds"]))
COMMON = " --sidis-request " + shlex.quote(SIDIS) + " --sia-request " + shlex.quote(SIA) + " --native-build " + shlex.quote(NATIVE) + " --parton-build " + shlex.quote(PARTON) + " --run-dir " + shlex.quote(RUN) + " --julia " + shlex.quote(str(JULIA)) + TIMEOUT

def index(key):
    if key not in CASES:
        raise ValueError("exact four storage cases required")
    return CASES.index(key) + 1

rule all:
    input:
        [RUN + "/" + c + "/job.completed" for c in CASES]
    run:
        build_or_verify_ledger(SIDIS, SIA, NATIVE, PARTON, RUN, RUN + "/storage-ledger.json", julia=JULIA, timeout_seconds=config.get("timeout_seconds"))

rule historical_candidates:
    input:
        manifest=str(ROOT / "config/historical-sources.json"), code=CODE
    output:
        RUN + "/historical/candidates.json", RUN + "/historical/report.json"
    params:
        root=str(ROOT), cache=CACHE, folder=RUN + "/historical", seed=SEED, offline=OFFLINE
    shell:
        "cd {params.root:q} && python -m collinsflow.cli historical-intake --manifest {input.manifest:q} --cache {params.cache:q} --output {params.folder:q} {params.seed} {params.offline}"

rule source_requests:
    input:
        candidates=RUN + "/historical/candidates.json", report=RUN + "/historical/report.json", code=CODE
    output:
        sidis=SIDIS, sia=SIA
    params:
        root=str(ROOT), manifest=str(ROOT / "config/historical-sources.json"), cache=CACHE, primary=PRIMARY,
        candidates=RUN + "/historical"
    shell:
        "cd {params.root:q} && python -m collinsflow.cli historical-sidis-point-request --manifest {params.manifest:q} --cache {params.cache:q} --candidates {params.candidates:q} --primary-cache {params.primary:q} --output {output.sidis:q} && python -m collinsflow.cli historical-sia-point-request --manifest {params.manifest:q} --cache {params.cache:q} --candidates {params.candidates:q} --output {output.sia:q}"

rule storage_case:
    input:
        sidis=SIDIS, sia=SIA, native=NATIVE, parton=PARTON, code=CODE,
        previous=lambda w: [] if index(w.case)==1 else [RUN + "/" + CASES[index(w.case)-2] + "/job.completed"]
    output:
        touch(RUN + "/{case}/job.completed")
    wildcard_constraints:
        case="case0[1-4]"
    params:
        root=str(ROOT), common=COMMON, case_index=lambda w:index(w.case)
    shell:
        "cd {params.root:q} && python scripts/run_fixed_context_storage.py run --case-index {params.case_index} {params.common}"
