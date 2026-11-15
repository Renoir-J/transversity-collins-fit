# KPSY partial tensor-charge benchmark

This is a **reconstructed-theory component benchmark** against [KPSY15 v1, Eqs. 127–131 and Table 1](https://arxiv.org/html/1505.05589v1). It integrates unweighted collinear `h1` over `0.0065 <= x <= 0.35` at **Q² = 10 GeV²**, starting at Q₀² = 2.4 GeV². The historical family sets antiquark transversity to zero. This is a partial moment in the measured-x window, with no full-x extrapolation. The future tensor-charge product at Q = 2 GeV has a separate reporting scale.

The [scientific contract](../config/kpsy-paper-moment-contract.json) was frozen before calculating candidate moments. The primary HTML version, complete byte hash and 25 mathematical anchors lock the definition, scale, printed centers and transversity parameters. A source receipt verifies acquisition only. If the downloaded bytes change, acquisition fails; it does not silently accept a regenerated or different paper.

## Two numerical paths

The native path uses the retained, pinned HOPPET LO transversity solver with the fixed-nf precision controls. Its initial `h1` uses strict CT10/DSSV08 Soffer inputs and the published rounded parameters. The independent dense Julia path uses a different evolution discretization and direct BigFloat powers for the initial formula. These paths share the locked physical input sources; they are not independent experimental determinations or recovered author code.

Each path refines its own grid under a fixed budget. The moment is the exact integral of its piecewise-linear interpolant in x, including partial boundary cells. A separate log-x QuadGK calculation checks the integration Jacobian, and the native adjoint supplies a reusable moment functional. Refinement and the independent evolution comparison must each meet a relative tolerance of `1e-4`. Coarse unresolved differences remain in the report. All numerical predictions and parameter sensitivity calculations run in Julia; Python acquires sources and launches whole jobs.

## Printed-value comparison and current result

The center-value criterion is an absolute difference of at most **0.005**, half a last printed digit, from `u = +0.30` and `d = -0.20`. The paper's statistical intervals and MINUIT errors are not numerical tolerances. Numerical convergence and agreement with the printed centers have separate report fields; both are required for an overall PASS.

The first complete Linux calculation gave:

| Quantity | Native, 1024 intervals | Dense reference, 3201 background nodes | Printed center | Raw-center check |
|---|---:|---:|---:|---|
| Partial up moment | +0.28854591 | +0.28854600 | +0.30 | FAIL |
| Partial down moment | −0.21285772 | −0.21285530 | −0.20 | FAIL |

All 22 other checks passed, including numerical refinement, independent integration, the native adjoint, initial formulas, endpoints and source/build stability. **The scientific report remains FAIL**, with exit code 2 and both failed raw-center checks retained. No native precision error was emitted by this calculation. These facts do not identify the cause of the discrepancy.

The printing sensitivity diagnostic evaluates all 64 vertices of the six-parameter rounding box. `Nd = -1.0` has one printed decimal; its interval is intersected with the physical domain `[-1,1]`. The sampled ranges were `[0.28062927, 0.29652320]` for up and `[-0.21353680, -0.20157143]` for down. They show that printing precision matters, but are not certified continuous extrema, recovered optimizer states or fit uncertainty intervals. They cannot rescue the failed raw-center criterion. Exact historical parameters and input versions remain unresolved.

## Running and evidence

After installing the locked Linux environment:

```bash
python scripts/prepare_kpsy_paper.py
python scripts/prepare_hoppet.py
python scripts/prepare_parton_inputs.py
julia --project=numerics numerics/bin/check_kpsy_paper_moments.jl \
  .runtime/hoppet/build.json .runtime/parton-inputs/build.json \
  .runtime/paper-moments/source/1505.05589v1.html \
  runs/paper-moments/validation.json
```

Use a new run directory for each calculation. Reports include all check operands, coarse and final grids, initial/evolved arrays, moment weights, the native moment functional, rounding vertices, source/build hashes and limitations. Julia matrices serialize as outer arrays of columns. CI preserves failed reports and returns failure when the scientific criterion fails; a red paper-benchmark job must not be described as a validated historical reproduction.

The original KPSY fitting program and modified HOPPET have not been located. This two-moment check does not reproduce the historical χ², SIDIS/SIA/STAR data, covariance, uncertainty coverage or complete twist-3 evolution. It does not populate a production validation gate, admit a new fit, establish G2PERF or change any production-ready flag. The synthetic closure workflow remains a separate software check.
