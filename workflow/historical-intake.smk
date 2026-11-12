"""Complete historical candidate job; numerical fitting is a separate boundary."""
from pathlib import Path
import shlex
from collinsflow.historical import verify_historical_intake
from collinsflow.provenance import computational_sources

ROOT = Path(workflow.basedir).parent
MANIFEST = str(ROOT / config.get("manifest", "config/historical-sources.json"))
CACHE = config.get("cache", "data/raw/historical")
OUTPUT = config.get("output", "data/normalized/historical")
SEED_ROOT = config.get("seed_root")
OFFLINE = config.get("offline", False)

rule all:
    input:
        f"{OUTPUT}/candidates.json",
        f"{OUTPUT}/report.json"
    run:
        verify_historical_intake(MANIFEST, CACHE, OUTPUT)

rule historical_candidates:
    input:
        manifest=MANIFEST,
        code=[str(path) for path in computational_sources(ROOT)]
    output:
        candidates=f"{OUTPUT}/candidates.json",
        report=f"{OUTPUT}/report.json"
    params:
        cache=CACHE,
        output=OUTPUT,
        root=str(ROOT),
        seed=("--seed-root " + shlex.quote(str(SEED_ROOT))) if SEED_ROOT else "",
        offline="--offline" if OFFLINE else ""
    threads: 1
    resources:
        mem_mb=512
    shell:
        "cd {params.root:q} && python -m collinsflow.cli historical-intake --manifest {input.manifest:q} "
        "--cache {params.cache:q} --output {params.output:q} {params.seed} {params.offline}"
