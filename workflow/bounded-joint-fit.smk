"""Whole Julia bounded fitting jobs on the fixed explicitly synthetic problem."""
from pathlib import Path

ROOT = Path(workflow.basedir).parent
RUN = config.get("run_dir", "runs/bounded-joint-fit")
JULIA = config.get("julia", "julia")
ALGORITHM = config.get("algorithm")
BUDGET = config.get("maxeval")
if config.get("mode", "synthetic_closure") != "synthetic_closure" or config.get("scope", "synthetic_closure") != "synthetic_closure":
    raise ValueError("bounded fitting workflow implements synthetic closure only")
if ALGORITHM not in ("LN_SBPLX", "LN_BOBYQA"):
    raise ValueError("explicit reviewed algorithm required")
if isinstance(BUDGET, bool) or not isinstance(BUDGET, int) or BUDGET not in (4000, 20000, 80000):
    raise ValueError("explicit frozen resource-ladder maxeval required")
CODE = [str(path) for folder in ("numerics/src", "numerics/bin", "numerics/test", "config")
        for path in sorted((ROOT / folder).glob("*")) if path.suffix in (".jl", ".json")]
CODE += [str(ROOT / name) for name in ("numerics/Project.toml", "numerics/Manifest.toml",
         "workflow/bounded-joint-fit.smk", "scripts/check_cached_bounded_workflow.py")]

rule all:
    input:
        f"{RUN}/accepted.json"
    params:
        root=str(ROOT), julia=JULIA, run=RUN, algorithm=ALGORITHM, maxeval=BUDGET
    shell:
        "cd {params.root:q} && export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 && "
        "{params.julia:q} --startup-file=no --project=numerics numerics/bin/check_bounded_joint_cache.jl "
        "{params.run:q}/validation.json {params.run:q}/audit.json {input:q} {params.algorithm:q} {params.maxeval}"

rule fit_and_audit:
    input:
        CODE
    output:
        protected(f"{RUN}/accepted.json")
    log:
        job=f"{RUN}/job.log",
        report=f"{RUN}/validation.json",
        audit=f"{RUN}/audit.json",
        controls=f"{RUN}/report-controls.log"
    params:
        root=str(ROOT), julia=JULIA, run=RUN, algorithm=ALGORITHM, maxeval=BUDGET
    threads: 1
    resources:
        mem_mb=4096
    shell:
        "cd {params.root:q} && test ! -e {log.report:q} && test ! -e {log.audit:q} && "
        "mkdir -p {params.run:q} && export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 && "
        "{params.julia:q} --startup-file=no --project=numerics numerics/bin/check_bounded_joint_fit.jl "
        "{params.algorithm:q} {params.maxeval} {log.report:q} >{log.job:q} 2>&1 && "
        "{params.julia:q} --startup-file=no --project=numerics numerics/bin/audit_bounded_joint_fit.jl "
        "{log.report:q} {log.audit:q} >>{log.job:q} 2>&1 && "
        "{params.julia:q} --startup-file=no --project=numerics numerics/test/bounded_report_controls.jl "
        "{log.report:q} >{log.controls:q} 2>&1 && cp {log.audit:q} {output:q}"
