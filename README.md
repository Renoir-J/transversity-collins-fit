# Transversity / Collins fitting workflow

This project uses a **Julia numerical core with Python/Snakemake job orchestration** (ADR-001). It targets joint transversity and Collins fits using SIDIS, SIA, STAR hadron-in-jet measurements and optional lattice tensor-charge constraints.

The current implementation supports **synthetic closure runs; production physics fits remain blocked**. The 66-point fixture exercises predictions, covariance and nuisance handling, derivatives, multistart optimization, replicas, holdouts and artifact validation. It contains no experimental measurements. The original KPSY15 fitting program and modified HOPPET implementation have not been located.

Experimental intake is now available separately: it acquires 22 pinned sources and normalizes 773 amplitudes across 82 overlapping candidate views. Candidates preserve source locators, units, statistical/systematic components, bin corrections and unresolved admission requirements. They are rejected by the fitting path. See the [experimental data review](docs/experimental-data.md).

## Running the workflow

Use Linux or WSL2 with Python 3.12 and Julia 1.12.6. Keep environments and run directories on the Linux filesystem. Each numerical job uses one process and one Julia, BLAS and OpenMP thread.

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.lock
python -m pip install --no-deps -e .
julia --project=numerics -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
snakemake --snakefile workflow/Snakefile --cores 2 --config mode=closure run_dir=runs/closure
```

Outputs include stage receipts, validation reports, fit and replica arrays, per-bin CSV predictions, a closure figure and `products/summary.json`. Failed or blocked stages stop their descendants, with failure records retained in `failures/`. Use a new run directory for each version. `collinsflow run --mode closure --run-dir runs/local` executes the same DAG sequentially.

```bash
python -m pytest -q
python scripts/plan_workflow.py validate
julia --project=numerics numerics/test/runtests.jl
make -C fortran
julia --project=numerics numerics/bin/check_bridge.jl "$PWD/fortran/libcollins_operator.so"
collinsflow readiness  # Expected exit code 2 while production prerequisites are missing
```

CI runs the tests, Fortran ABI checks and complete Snakemake closure workflow on Linux. Run artifacts are retained for 30 days. Changes to `main` use pull requests and squash merges; see the [Git policy](docs/git-policy.md).

Acquire and validate experimental candidates with `snakemake --snakefile workflow/experimental-intake.smk --cores 2`. Raw and normalized outputs remain outside Git. A separate source-intake workflow verifies public downloads and hashes when the source manifest or parser changes.

## Project status and interfaces

Track progress through [Issues](https://github.com/Renoir-J/transversity-collins-fit/issues), [Milestones](https://github.com/Renoir-J/transversity-collins-fit/milestones) and the [implementation status](docs/implementation-status.md). Production admission is governed by the [physics contract](config/physics-contract.json). Synthetic receipts do not satisfy production validation requirements.

- [Numerical architecture and validation](docs/numerical-workflow.md)
- [Numerical artifact and data interfaces](docs/kernel-contract.md)
- [Twelve-stage DAG and research scenarios](config/workflow.json)
- [Dataset selection and overlap policy](config/data-registry.json)
- [Pinned source inventory](references/source-inventory.json)

Different theory approximations require separate benchmarks. Keep generated reports, downloaded sources, local instructions, internal notes, credentials and runtime environments out of Git. Core source code, environment lock files, reviewed constraints and synthetic fixtures remain versioned. See [third-party notices](THIRD_PARTY_NOTICES.md) for source attribution.
