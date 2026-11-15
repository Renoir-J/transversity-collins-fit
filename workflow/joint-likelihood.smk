"""Complete Julia validation jobs; synthetic observations never admit production."""
from pathlib import Path

ROOT = Path(workflow.basedir).parent
RUN = config.get("run_dir", "runs/joint-likelihood")
NATIVE = config.get("native_build", ".runtime/hoppet/build.json")
PARTON = config.get("parton_build", ".runtime/parton-inputs/build.json")
JULIA = config.get("julia", "julia")
CODE = [str(p) for folder in ("numerics/src", "numerics/bin", "numerics/test", "config")
        for p in sorted((ROOT / folder).glob("*")) if p.suffix in (".jl", ".json")]
CODE += [str(ROOT / p) for p in ("numerics/Project.toml", "numerics/Manifest.toml",
          "workflow/joint-likelihood.smk")]

rule all:
    input:
        expand(f"{RUN}/{{kind}}/accepted.json", kind=["analytic", "inputs"])

rule validate_joint:
    input:
        code=CODE,
        builds=lambda w: [NATIVE, PARTON] if w.kind == "inputs" else []
    output:
        protected(f"{RUN}/{{kind}}/accepted.json")
    log:
        job=f"{RUN}/{{kind}}/job.log",
        report=f"{RUN}/{{kind}}/validation.json",
        audit=f"{RUN}/{{kind}}/audit.json",
        controls=f"{RUN}/{{kind}}/report-controls.log"
    wildcard_constraints:
        kind="analytic|inputs"
    params:
        root=str(ROOT), julia=JULIA,
        directory=lambda w: f"{RUN}/{w.kind}"
    threads: 1
    resources:
        mem_mb=4096
    shell:
        "cd {params.root:q} && test ! -e {log.report:q} && test ! -e {log.audit:q} && "
        "mkdir -p {params.directory:q} && "
        "export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 && "
        "{params.julia:q} --startup-file=no --project=numerics numerics/bin/check_joint_likelihood.jl "
        "{wildcards.kind:q} {input.builds:q} {log.report:q} >{log.job:q} 2>&1 && "
        "{params.julia:q} --startup-file=no --project=numerics numerics/bin/audit_joint_likelihood.jl "
        "{log.report:q} {log.audit:q} >>{log.job:q} 2>&1 && "
        "{params.julia:q} --startup-file=no --project=numerics numerics/test/joint_report_controls.jl "
        "{log.report:q} >{log.controls:q} 2>&1 && cp {log.audit:q} {output:q}"
