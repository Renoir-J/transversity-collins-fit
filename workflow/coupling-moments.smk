"""Whole Julia jobs only. A failed scientific diagnostic stays a failed workflow."""
from pathlib import Path

ROOT = Path(workflow.basedir).parent
RUN = config.get("run_dir", "runs/coupling-moments")
NATIVE = config.get("native_build", ".runtime/native-coupling/build.json")
PARTON = config.get("parton_build", ".runtime/parton-inputs/build.json")
PAPER = config.get("paper", ".runtime/native-coupling/paper.html")
JULIA = config.get("julia", "julia")
CODE = [str(p) for folder in ("numerics/src", "numerics/bin", "config")
        for p in sorted((ROOT / folder).glob("*")) if p.suffix in (".jl", ".json")]
CODE += [str(ROOT / p) for p in ("numerics/Project.toml", "numerics/Manifest.toml",
          "workflow/coupling-moments.smk", "fortran/hoppet_bridge.f90",
          "fortran/native_coupling_reference.f90", "bridges/lhapdf.cpp",
          "src/collinsflow/hoppet.py", "src/collinsflow/parton.py",
          "src/collinsflow/provenance.py", "scripts/prepare_hoppet.py",
          "scripts/prepare_native_coupling.py", "scripts/prepare_kpsy_paper.py",
          "scripts/prepare_parton_inputs.py")]

rule all:
    input:
        f"{RUN}/sensitivity/accepted.json"

rule coupling_reference:
    input:
        build=NATIVE,
        code=CODE
    output:
        protected(f"{RUN}/coupling/accepted.json")
    log:
        job=f"{RUN}/coupling/job.log",
        report=f"{RUN}/coupling/validation.json"
    params:
        root=str(ROOT), julia=JULIA, directory=f"{RUN}/coupling"
    threads: 1
    shell:
        "cd {params.root:q} && test ! -e {log.report:q} && mkdir -p {params.directory:q} && "
        "JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "
        "{params.julia:q} --startup-file=no --project=numerics numerics/bin/check_native_coupling.jl "
        "{input.build:q} {log.report:q} >{log.job:q} 2>&1 && cp {log.report:q} {output:q}"

rule prepared_time_controls:
    input:
        build=NATIVE,
        coupling=f"{RUN}/coupling/accepted.json",
        code=CODE
    output:
        protected(f"{RUN}/prepared/accepted.json")
    log:
        job=f"{RUN}/prepared/job.log",
        report=f"{RUN}/prepared/validation.json"
    params:
        root=str(ROOT), julia=JULIA, directory=f"{RUN}/prepared"
    threads: 1
    shell:
        "cd {params.root:q} && test ! -e {log.report:q} && mkdir -p {params.directory:q} && "
        "JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "
        "{params.julia:q} --startup-file=no --project=numerics numerics/bin/check_prepared_times.jl "
        "{input.build:q} {input.coupling:q} {log.report:q} >{log.job:q} 2>&1 && cp {log.report:q} {output:q}"

rule physical_sensitivity:
    input:
        build=NATIVE,
        parton=PARTON,
        paper=PAPER,
        coupling=f"{RUN}/coupling/accepted.json",
        prepared=f"{RUN}/prepared/accepted.json",
        code=CODE
    # Validation is a log so Snakemake retains failed scientific evidence.
    # Only a successful process may create the separate accepted output.
    output:
        protected(f"{RUN}/sensitivity/accepted.json")
    log:
        job=f"{RUN}/sensitivity/job.log",
        report=f"{RUN}/sensitivity/validation.json"
    params:
        root=str(ROOT), julia=JULIA, directory=f"{RUN}/sensitivity"
    threads: 1
    resources:
        mem_mb=4096
    shell:
        "cd {params.root:q} && test ! -e {log.report:q} && mkdir -p {params.directory:q} && "
        "JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "
        "{params.julia:q} --startup-file=no --project=numerics numerics/bin/check_coupling_moments.jl "
        "{input.build:q} {input.parton:q} {input.paper:q} {input.coupling:q} {input.prepared:q} "
        "{log.report:q} >{log.job:q} 2>&1 && cp {log.report:q} {output:q}"
