# Native running-coupling diagnostic

This **reconstructed-theory comparator** checks the current pinned HOPPET coupling against an independent Julia solution. It does not recover the original KPSY implementation or select a production prescription. [KPSY15 §III.1](https://arxiv.org/html/1505.05589v1) specifies HOPPET alpha_s(Q0)=0.327 at Q0²=2.4 GeV², together with a two-loop nf=5, Lambda=0.225 GeV prescription and alpha_s(MZ)=0.118. These statements do not identify all historical backend choices.

The [prospective contract](../config/native-coupling-contract.json) was committed before any native candidate calculation. It fixes the primary version/byte hash and five math anchors, upstream commit/tree, beta equation, scales, precision, refinement budgets and thresholds. It preserves the distinction in [COLLINEAR-LO-002](../config/evolution-contract.json) between the explicit asymptotic 1/L expansion and the differential-equation solution.

## Reference and independent calculation

An opt-in extension to the existing native builder compiles `fortran/native_coupling_reference.f90` against unchanged pinned HOPPET. The ordinary native bridge and its default build mode are retained. Each comparator batch runs in a fresh process with a local `running_coupling` object: `InitRunningCoupling`, `Value` and `Delete`. It does not initialize a grid or discarded splitting functions. The pinned `qcd_coupling` implementation delegates to `new_as`; its old analytic entry is disabled.

The independent Julia reference solves the fixed-nf two-loop equation

\[
\frac{d\alpha}{d\ln Q^2}=-B_0\alpha^2-B_1\alpha^3,
\quad B_0=\frac{\beta_0}{4\pi},\quad B_1=\frac{\beta_1}{16\pi^2}.
\]

With c=B1/B0, the monotonic implicit function `F(alpha)=1/alpha+c*log(alpha/(1+c*alpha))` is solved by 256-bit bisection. A separate analytic antiderivative gives tau=integral alpha_s/(2*pi) dln(Q²). Native alpha samples feed Julia Gauss quadrature in log(Q²), with no Python/model callbacks. Julia also computes analytic Mellin evolution-factor controls, composition and all differences between prescriptions. Settings remain fixed; no fitted-setting derivatives are admitted.

The comparator uses a strict, complete batch protocol and emits values only after all native evaluations and cleanup succeed. Requests, responses, stderr, return codes and hashes are preserved. Invalid headers, scales, loop orders, nf, dt, nonfinite settings, truncation and trailers must reject. A zero native precision-error count is required.

## Scoped result and limits

The first complete Linux run passes **36 comparator checks**, including sampled alpha/time agreement, adjacent refinements, initial anchors, nf3–nf6 one-loop controls, Mellin factors, composition, nine rejected native requests and source/build stability. At native dt=0.05 with 64-point quadrature, the largest selected alpha difference is 3.53e-14 scaled and the time difference 7.98e-15. The dt=0.2 and 0.1 levels remain unaccepted under the frozen 1e-13 anchor criterion. No threshold changed.

The three definitions give the following diagnostic values:

| Definition | alpha_s at Q²=10 | tau from Q0²=2.4 to 10 | alpha_s at MZ=91.1876 GeV |
|---|---:|---:|---:|
| Existing asymptotic Lambda=0.225 | 0.2458316670 | 0.064089238695 | 0.1179026616 |
| Native RG, literal alpha0=0.327 | 0.2482642732 | 0.064469877780 | 0.1189679196 |
| Native RG, matched asymptotic alpha0 | 0.2482554605 | 0.064467207084 | 0.1189659913 |

These are different definitions; their difference is not a numerical error tolerance or a recovered historical setting. The matched-initial-value diagnostic separates the running-definition difference from initial-value rounding. Deviations from the paper's printed MZ value are recorded without declaring historical compatibility. Flavor thresholds and original backend flags remain unresolved.

No physical partial moment was evaluated here. The failed printed-center checks in [draft #35](https://github.com/Renoir-J/transversity-collins-fit/pull/35) remain unchanged. This result does not establish the cause of that discrepancy, replace its contract, validate full-x extrapolation, certify a fit/coverage result or pass G2PERF. Sampled coupling comparisons are not a uniform certificate over arbitrary scales. The fitting core, numerical predictions and uncertainty work remain governed by ADR-001 and their separate admission gates.

## Running

Use the locked Linux environment:

```bash
python scripts/prepare_native_coupling.py
julia --project=numerics numerics/bin/check_native_coupling.jl \
  .runtime/native-coupling/build.json runs/native-coupling/validation.json
```

Use a new run directory for each job. The builder also supports an offline pinned source, explicit compiler/CMake paths and compiler flags. A build receipt is `BUILT_NOT_VALIDATED`; the separate numerical report validates only this comparator scope. CI archives the raw primary evidence, both compiled objects, upstream test log, all native batches and failed levels. Other runtime dependencies retain receipt linkage unless their bytes are separately archived and checked.
