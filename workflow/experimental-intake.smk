"""Pinned public-source acquisition and experimental candidates. No fitting is authorized."""
import re
from pathlib import Path
from collinsflow.experimental import acquire, verify_intake
from collinsflow.provenance import atomic_json, computational_sources, file_hash, read_json

ROOT = Path(workflow.basedir).parent
MANIFEST = str(ROOT / config.get("manifest", "config/experimental-sources.json"))
SOURCES = {source["id"]: source for source in read_json(MANIFEST)["sources"]}
CACHE = config.get("cache", "data/raw/intake")
OUTPUT = config.get("output", "data/normalized/intake")
SEED_ROOT = config.get("seed_root")

rule all:
    input:
        f"{OUTPUT}/candidates.json",
        f"{OUTPUT}/intake-report.json"
    run:
        verify_intake(MANIFEST, CACHE, OUTPUT)

rule locked_source:
    input:
        manifest=MANIFEST
    output:
        receipt=f"{OUTPUT}/sources/{{source_id}}.json"
    wildcard_constraints:
        source_id="|".join(re.escape(key) for key in SOURCES)
    threads: 1
    resources:
        mem_mb=256
    run:
        source = SOURCES[wildcards.source_id]
        seed = Path(SEED_ROOT) / source["seed_path"] if SEED_ROOT else None
        path = acquire(source, CACHE, seed=seed)
        atomic_json(output.receipt, {"id": source["id"], "path": str(path), "sha256": file_hash(path),
                                    "scope": "experimental_candidate"})

rule candidates:
    input:
        manifest=MANIFEST,
        locks=[f"{OUTPUT}/sources/{key}.json" for key in SOURCES],
        code=[str(path) for path in computational_sources(ROOT)]
    output:
        candidates=f"{OUTPUT}/candidates.json",
        report=f"{OUTPUT}/intake-report.json"
    params:
        cache=CACHE,
        output=OUTPUT,
        root=str(ROOT)
    threads: 1
    resources:
        mem_mb=1024
    shell:
        "cd {params.root:q} && python -m collinsflow.cli intake --manifest {input.manifest:q} "
        "--cache {params.cache:q} --output {params.output:q} --offline"
