# Evolved polarized matching rows

[KPSY-MATCHED-ROWS-001](../config/matched-rows-contract.json) combines the separately contracted LO evolution and first-order Wilson matching for polarized legs. It implements SIDIS transversity and SIDIS/SIA homogeneous KPSY Hhat3. The source is [KPSY15 Eqs. 64, 67, 68 and 105](https://arxiv.org/html/1505.05589v1). This is reconstructed theory, with production_ready=false. Original author code and the historical fitted baseline remain unavailable.

For a fixed output scale Q2, the unweighted matched leg is

\[
\left(1+\frac{\alpha_s(Q^2)}{\pi}d_{\rm process}\right)U(\tau)h_0
+\frac{\alpha_s(Q^2)}{\pi}R\otimes U(\tau)h_0.
\]

The delta constant is −2CF for SIDIS and CF(π²−8)/4 for SIA. Transversity has R=0 and is admitted only as a SIDIS leg. The Collins regular coefficient is R(v)=2CF v ln(v)/(1−v), with zero virtual/delta pieces. Multiplication by ln(v) removes the plus subtraction at v=1. There is no polarized gluon matching. Inputs remain dimensionful KPSY Hhat3, with the existing Trento sign/z/mass conversion; no extra mass or z factors enter this collinear matching operation.

Matching uses the **physical** contracted alpha_s at every supplied Q2. The auxiliary constant alpha=1 used by retained LO evolution only parametrizes the evolution time and never becomes a Wilson coefficient. If these scales are mu_b squared, the caller must supply those exact fixed scales. No interpolation or hidden scale prescription is inferred.

## Retained row construction

`HoppetMatchingPlan` remains distinct from evolution-only `HoppetPointPlan`. Each required fraction is the exact boundary of its own single negative-order grid. With D=diag(x), the regular unweighted row is

\[
D\,U_w(\tau)^T R_w^T\frac{e_i}{x_i},
\]

projected onto the constrained h(1)=0 subspace. The adapter applies the retained convolution transpose for R once, then for each retained U. At tau=0 it uses the R cotangent directly. It allocates neither a dense x-by-x matrix nor a nodes-by-scales-by-scales cotangent tensor. The regular kernel includes HOPPET's additional v weight exactly once. Nested-grid matching composition is rejected in this contract.

Preparation compares the combined rows with a full retained R*(U*h) calculation on signed analytic inputs before closing the native context. Repeated fractions share preparation and preserve requested order; fraction one returns zero. Process and physical-alpha metadata are explicit, with fixed-setting mutation checks. Julia generic matrix products evaluate the fitted initial family, gradients and nested Hessians after preparation. There are no native/Python callbacks in numerical iterations. Fitting coupling, scales or grids would require dU and dC and is rejected.

```julia
plan = HoppetMatchingPlan(".runtime/hoppet/build.json", [0.05, 0.3, 0.7, 0.9, 1.0],
                         [0.6, 2.4, 10.0, 100.0], 2.4, AsymptoticTwoLoopCoupling();
                         process=:SIDIS, quantity=:KPSY_Hhat3_homogeneous, intervals=1024)
# Supply unweighted initial Hhat3 on plan.grid; axes are query, initial column, scale.
# All native contexts have closed. This plan supplies one matched leg, not an observable.
```

## Independent checks and limits

The independent Mellin reference uses R_N=−2CF trigamma(N+1), an analytic beta-polynomial initial moment, and a separate inversion contour. Complex trigamma is checked against direct moments, recurrence and contour changes. Q0 matching is compared with independent quadrature. Full constrained-basis rows, transpose identities, physical-coupling selection, generic gradients, Hessians/HVPs and invalid-setting rejection are checked explicitly.

These in-memory rows are not an admitted serialized production kernel. Source-bound artifact preparation and complete coefficient-mutation admission remain required alongside the existing metadata checks.

The pinned-input job uses strict CT10/DSSV08/DSS14 initial inputs and the printed historical parameter family. It refines native grids, checks direct Q0 matching split at source knots, and independently evolves source/probe-aligned dense grids with cell-wise matching quadrature. Coarse failed comparisons remain recorded; only unconverged cases receive the contracted extra dense-grid refinement, with unchanged value/derivative targets. Reference finite differences use 256-bit initial-family/row arithmetic; fitting remains Float64. These checks do not calibrate physical interval coverage or admit a likelihood.

The initial physical run at commit `0b57abd` failed one of 93 checks: SIA Hhat3 at Q2=0.6, z=0.05 in the u channel differed by 1.315e-4 between 512 and 1024 intervals, above the unchanged 1e-4 target. That report and contract snapshot remain preserved. The physical runner therefore has an explicit 256/512/1024/2048/4096 budget, stopping each process/quantity only when all supplied scales pass its final adjacent-grid comparison at or above 1024 intervals. Unconverged cases refine; budget exhaustion is FAIL. Constructor defaults and analytic checks alone never certify an arbitrary physical job.

```sh
python scripts/prepare_parton_inputs.py
python scripts/prepare_hoppet.py
julia --project=numerics numerics/bin/check_matched_rows.jl .runtime/hoppet/build.json runs/matched/rows.json
julia --project=numerics numerics/bin/check_matched_inputs.jl \
  .runtime/hoppet/build.json .runtime/parton-inputs/build.json runs/matched/inputs.json
```

The 64-scale analytic resource job records actual cold preparation, warm batched application, environment and coefficient/output bytes. A separate 2336-Q0-scale matching/storage check exercises the linear row interface. Neither is physical 2336-scale throughput, peak RSS, a full likelihood benchmark or G2PERF. Native quadrature warnings stay in the logs; requested epsilon is not an attained accuracy bound. The independent final comparisons determine the scoped component result.

Connecting both matched legs to CSS/radial/flavor assembly, target/axis/angular and publication-specific observables, source-scale policies, the historical benchmark, data covariance/overlap and production derivative/performance/coverage admission remain open. In particular, this component does not supply Belle thrust-frame factorization, STAR jet matching, a data likelihood or a new fit. The earlier rejected global-grid evolution pilot is not admitted by this separate matching contract.
