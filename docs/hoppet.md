# Retained HOPPET evolution adapter

The project-owned `fortran/hoppet_bridge.f90` builds the LO transversity and explicitly selected homogeneous KPSY Hhat3 kernel using the unchanged, pinned HOPPET convolution, interpolation and evolution routines. The contract is [KPSY-HOPPET-001](../config/hoppet-contract.json), extending the same equation in [COLLINEAR-LO-002](../config/evolution-contract.json). This is reconstructed theory. It is not the recovered author's modified HOPPET, a historical fit, or production admission.

The upstream source is [HOPPET commit feda4f2](https://github.com/hoppet-code/hoppet/tree/feda4f2f64d0131f489ad1cf7728ec6be186cc56). Cite [Salam and Rojo](https://arxiv.org/abs/0804.3755) and [the HOPPET v2 paper](https://arxiv.org/abs/2510.09310) when using it. The kernel and Collins approximation follow [KPSY15 Eqs. 64 and 69–71](https://arxiv.org/html/1505.05589v1). `external/` remains pinned and unchanged. The builder acquires a complete clean checkout under ignored `.runtime`, hashes source files and symlinks, compiles fresh objects and records the native library and runtime dependencies. `BUILT_NOT_VALIDATED` is a build receipt, not a numerical PASS.

## Weights and time

The Julia interface accepts unweighted initial ordinates in increasing x. HOPPET receives `x*h(x)` on its native y grids, with `y=-log(x)`. Non-delta splitting pieces carry HOPPET's additional z Jacobian. Outputs are divided by x once. A Trento Collins first moment cannot be supplied as Hhat3 without the separate sign, z and mass conversion.

For the fixed-nf LO equation, all coupling dependence enters through

\[
\tau=\int_{Q_0^2}^{Q^2}\frac{\alpha_s(q^2)}{2\pi}\,d\log q^2,
\qquad h(Q^2)=\exp(\tau P)h_0.
\]

Julia evaluates this integral using the selected one-loop or explicit asymptotic two-loop contract. The retained solver uses an **auxiliary** constant coupling of one with `Q_start=100` and `Q_end=100*exp(pi*tau)`. This interval has exactly the same evolution time. These scales and coupling are a change of integration variable, not physical inputs. The solver explicitly uses evolution `nloop=1` while the auxiliary coupling uses `nloop=0`. No HOPPET running-coupling definition replaces the contracted physical coupling. Flavor thresholds remain unsupported.

In the pinned source, `ev_evolve` accepts `du` but omits it when invoking `ev_SetModuleConsts`. The adapter therefore sets the public `SetDefaultEvolutionDu` during construction and restores the previous value afterward. An explicit check requires a measurable coarse-step change as well as convergence; identical results from an ineffective setting are not accepted as refinement evidence. The initial failed derivative report and ineffective-step observations are retained.

## Grid, cache and transpose

Single uniform y grids and locked nested grids use negative HOPPET interpolation orders. The nested layout adds small-y resolution for slow endpoint shapes, rather than allocating a dense x-by-x matrix at every scale. Exact native x values form the returned sorted union; no rounding merges distinct floating point knots. Shared nodes use their finest native representation. Upstream grid locking adjusts subgrid steps and extents; its warnings and requested settings remain visible.

The declared domain endpoint is canonical: with `xmin=0.05`, that exact Float64 value is sampled and queryable even when `exp(log(0.05))` rounds one ULP away. This convention only defines the domain endpoint, with `ymax=-log(xmin)`; it does not alter interior source knots or admit extrapolation. A query one ULP below xmin is rejected.

One scalar `grid_conv` is retained for each supplied fixed time, extracted from the nonsinglet operator. At LO there is no flavor mixing. Scales are not implicitly interpolated. The reported coefficient payload counts stored convolution weights; it excludes output arrays, AD batches and compiler/library memory and is not peak process memory.

For each native negative-order subgrid, the convolution is lower-triangular Toeplitz, so `Cᵀ=R*C*R`. The adjoint reuses HOPPET convolution with reversed arguments. For nested grids, fine-to-coarse decanting is transposed in reverse assignment order. Node gather/selection and both x weights are also transposed, then the constrained `h(1)` coordinate is removed. Complete constrained-basis comparisons and dot products check this construction; no forward convolution or solver is rewritten.

## Julia derivatives and state

`HoppetPlan` fixes the grid, scales, coupling, nf and tolerances before optimization. ForwardDiff primal/tangent components are batched through the same Float64 linear ABI, including nested dual components for Hessians. The bridge has no callbacks into Python or Julia. Fitted operator settings require `dU`; attempts to differentiate those settings are rejected.

Each loaded bridge permits one active context, with generation handles and a Julia process lock. Close it before opening another. Whole numerical jobs remain isolated processes. Upstream global/SAVE state is not treated as generally thread-safe. Queries below grid support, nonfinite inputs, nonzero input endpoints, closed contexts and mutated grids are rejected.

## Query-bound row caches

The initial global-grid physical-input pilot at commit `afb5c6d` fails the original `1e-4` target. The largest differences occur at source interpolation knots, including z=0.7 and 0.9, and are already visible at Q0. Its 59-check report retains 11 failures in both local Linux validation and [independent CI](https://github.com/Renoir-J/transversity-collins-fit/actions/runs/43167683823); it does not admit the global-grid path. These are errors in the numerical representation of the full initial functions, not evidence that their standalone fitted power factors describe the complete endpoint behavior.

[KPSY-HOPPET-POINTS-001](../config/hoppet-point-contract.json) instead puts each required query at the exact boundary of its own retained native grid. The checked native adjoint exports that query's fixed rows, with a full native-forward comparison before closing the context. The sorted union of initial coordinates is sampled once. Repeated queries reuse preparation and preserve their requested output order. Julia applies the detached rows with generic matrix products, so h1 and Hhat3 plans can coexist and numerical iterations have no foreign calls. This does not reimplement the native solver or forward convolution.

Each query still needs integration-grid refinement: alignment is not itself an accuracy certificate. The physical pilot retains the same `1e-4` value and `1e-7` derivative targets. Reference finite differences use 256-bit initial-family and fixed-row arithmetic to avoid subtracting rounded foreign outputs; the fit remains Float64. The failed Float64/global-grid pilot remains recorded. Only unconverged independent dense-reference cases receive additional grid refinement within a fixed budget; coarse failures stay visible.

On the Linux development run, 53 point-input checks passed over Q2=0.6, 2.4, 10 and 100 GeV2. The largest native refinement difference was 3.26e-5 and the largest comparison against the independently refined dense route was 2.35e-5. The independent h1 reference at Q2=100 failed its 3201-background-node refinement target (1.03e-4); 6401 nodes reduced that difference to 3.05e-5 with the original target unchanged. The finite-difference JVP/HVP differences were below 2.31e-9. These are scoped component results, not G2/G2AD, data admission or a physical extraction. Independent CI repeats the submitted source.

```julia
plan = HoppetPointPlan(".runtime/hoppet/build.json", [0.05, 0.3, 0.7, 0.9, 1.0],
                      [0.6, 2.4, 10.0, 100.0], 2.4, AsymptoticTwoLoopCoupling();
                      quantity=:KPSY_Hhat3_homogeneous, intervals=1024)
# Supply initial Hhat3 on plan.grid; output axes are query, initial column, scale.
# Every native context has already closed; the retained cache contains fixed rows.
```

```bash
julia --project=numerics numerics/bin/check_hoppet_points.jl \
  .runtime/hoppet/build.json runs/hoppet/points.json
julia --project=numerics numerics/bin/check_hoppet_inputs.jl \
  .runtime/hoppet/build.json .runtime/parton-inputs/build.json \
  runs/hoppet/physical-inputs.json query_rows
```

Omit `query_rows` to reproduce the unadmitted global-grid diagnostic. The point job separately measures a 64-scale analytic workload at checked accuracy. Stored coefficient bytes and warm row application are informational; they are not peak RSS, complete physical matching cost or G2PERF.

```bash
python scripts/prepare_hoppet.py
julia --project=numerics numerics/bin/check_hoppet.jl \
  .runtime/hoppet/build.json runs/hoppet/validation.json
```

```julia
using CollinsFit
plan = HoppetPlan(".runtime/hoppet/build.json", [2.4, 10.0, 100.0], 2.4,
                  AsymptoticTwoLoopCoupling(); quantity=:transversity)
try
    h0 = plan.grid .* (1 .- plan.grid).^3  # analytic input, not an extraction
    hQ = apply_operator(plan, h0)          # grid nodes by supplied scales
    values = hoppet_evaluate(plan, hQ, [0.01, 0.1, 0.5])
finally
    close(plan)
end
```

The validation job compares analytic generators, independent inverse Mellin solutions and moments, solver refinements, complete-basis transposes, nested AD and finite differences. Uniform-grid endpoint failures remain in local records; mesh refinement does not relax acceptance targets. Upstream quadrature can print `TOO HIGH ACCURACY REQUIRED` at very small cells. Requested epsilon is not claimed as an attained error bound: retain these warnings and require independent overall accuracy and epsilon-stability checks. Such component checks cannot certify production quadrature.

Physical source-knot interpolation, complete many-scale matching/observable assembly, historical backend comparison, interval coverage and production gradient/performance admission remain separate requirements. A polynomial comparison alone does not certify the published slow endpoint family or any physical fit.
