# Reconstructed parton inputs

These interfaces supply fixed collinear inputs for a future KPSY reconstruction. They do not implement a physical cross section, reproduce the historical fit or satisfy production admission. The [input contract](../config/parton-contract.json) and [source manifest](../config/parton-sources.json) fix the numerical conventions and raw SHA-256 values. Raw tables, library sources and compiled outputs remain ignored.

## Sources and interface

| Input | Retained implementation | Public quantity |
|---|---|---|
| CT10nlo member 0, DataVersion 4 | Unmodified LHAPDF 6.5.6 `logcubic`, `ForcePositive=1`, through a project C ABI | `parton_pdf` returns f; `ct10_batch` returns x*f |
| DSSV08 central NLO | Julia interpolation of the locked grid, compared with unchanged Fortran | `parton_helicity` returns unweighted u, d, ubar, dbar, gluon |
| DSS14 central pion NLO | Julia interpolation of the locked grid, compared with unchanged Fortran | `pion_ff` returns unweighted D for PDG flavor and charge -1/0/+1 |

The CT10 input comes from the [official LHAPDF set distribution](https://lhapdfsets.web.cern.ch/current/CT10nlo/CT10nlo.info). The library is the [official LHAPDF source release](https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.5.6.tar.gz). DSSV08 and DSS14 are secondary distributions of publication-tagged author routines, retained at the exact commits linked in the manifest. They are not recovered KPSY fitting code. The original KPSY modified HOPPET remains unavailable.

All query scales are Q² in GeV². The CT10 raw grid lists Q, which the comparison squares explicitly. Legacy routines return x times helicity or z times fragmentation; the public flavor interfaces remove that factor once. DSSV valence values require adding the corresponding antiquark before dividing by x. DSS14 contains quark-plus-antiquark totals and valence differences; heavy-quark outputs are each half the total. Pion minus swaps light quarks/antiquarks, and pion zero is the charge average. Public DSSV strange/heavy helicity mapping is withheld pending a separate source-convention check.

The Julia legacy interpolator removes the source endpoint powers, interpolates bilinearly in log(x or z), log(Q²), and restores the powers. It sets the final x/z=1 grid entry to zero, as the source does. DATA knots are parsed from the unchanged fixed-form source, preserving default REAL versus D-exponent precision. In particular DSS14 `**0.3` is Float32 promoted to Float64, equal to 0.30000001192092896. An initial all-Float64 decimal port failed the predeclared tolerance; matching the actual Fortran literal semantics fixed the discrepancy without changing thresholds.

## Range and differentiation policy

Default strict ranges are CT10 x≥1e-8, Q²≥1.3²; DSSV x≥1e-5, 0.8≤Q²≤1e6; DSS14 z≥0.05, 1≤Q²≤1e5. Upper x/z is 1; CT10 Q² upper bound is 1e10. The DSSV executable checks and knots supersede its older header's narrower scale description. DSS14's first grid knot is 0.01, but its documented validity begins at 0.05.

An explicitly selected `:reference` policy enables LHAPDF continuation or DSS14 source log-linear scale extrapolation. DSSV still rejects scales outside its executable checks; DSS14 still requires z≥0.05. `parton_domain` labels outside-range queries and always reports production_ready=false. No policy silently freezes scales. CSS μ_b can be below CT10/DSS14 strict limits; a historical/production prescription and sensitivity study are still required.

CT10 calls are serialized within each Julia process. Different source roots require isolated processes because LHAPDF caches metadata globally. Legacy references run as fresh executables, isolating COMMON/SAVE state. Kinematic AD through these fixed input interfaces is rejected; differentiable fitted parameterizations belong in Julia. No Python or model calls enter numerical evaluation. The CSS coupling remains the explicitly contracted coupling rather than LHAPDF's provider alpha_s.

The old CERN-style Fortran interpolators declare some dummy-array extents smaller than the actual backing grids; modern compilation warns about these legacy declarations. The pinned source is preserved and compiled without bounds instrumentation for reference comparison. Agreement on the reported checks does not establish that the legacy program is free of all undefined behavior or make it a production implementation.

## Reproduce the component comparison on Linux

```sh
python scripts/prepare_parton_inputs.py
julia --project=numerics numerics/bin/check_parton_inputs.jl
```

The builder verifies seven locked downloads, builds LHAPDF from a freshly extracted source tree with Python disabled, compiles the C ABI and separate factory executable, and compiles unchanged DSS reference routines with small project drivers. It records source, binary, configuration, compiler and command hashes. Builds must stay inside `.runtime/`. `--offline` uses verified cached or seed bytes. `--fc`, repeated `--fc-flag`, and `--fortran-runtime` support an isolated compiler/runtime; their choices enter the build receipt.

The validation job compares 99 DSSV and 216 DSS14 queries, including off-grid points, boundaries and all pion charges. It checks 42,900 CT10 raw node/flavor entries and 462 separate C++ factory queries, unweighting, PDG ordering, heavy flavors, strict rejection, labeled low-scale reference behavior, closed-handle safety, AD rejection and changed-source rejection. Tolerances are frozen in the contract. CT10's factory and ABI use the same official interpolator: their agreement tests the binding, while raw-node values provide an independent file-level check. Neither is an independent interpolation algorithm or a historical-fit benchmark.

The dedicated GitHub workflow reacquires all bytes and builds references on Linux. Reports and build evidence are retained as artifacts. A component PASS remains reconstructed_theory with production_ready=false; synthetic closure and validated production results are separate scopes.
