# Reconstructed collinear evolution

`numerics/src/CollinearLO.jl` implements one component of issue #4. It builds a fixed-nf LO collinear evolution matrix in Julia and applies that matrix in generic arithmetic. Its scientific contract is [COLLINEAR-LO-001](../config/evolution-contract.json). This is reconstructed theory tested on analytic inputs, separate from the synthetic fitting fixture and from validated production physics.

The supported quantities are unweighted spacelike transversity `h1(x)` and explicitly selected `KPSY_Hhat3_homogeneous`: the twist-three fragmentation function defined in [KPSY Eq. 56](https://arxiv.org/html/1505.05589v1#S2.E56), using the homogeneous approximation in Eqs. 64 and 69–71. Those equations justify the same collinear kernel for these two specifically defined quantities. A Trento Collins first moment has additional z/mass conversion factors in Eq. 58 and cannot be supplied as though it were the same quantity. Generic Collins and x-weighted distribution tags are rejected.

## Equation and discretization

The convention is

\[
\frac{d h(x,Q^2)}{d\ln Q^2}=\frac{\alpha_s}{2\pi}(P\otimes h)(x),\qquad
P(z)=C_F\left(2z\left[\frac{1}{1-z}\right]_++\frac32\delta(1-z)\right).
\]

Here the plus prescription acts on the denominator. The truncated convolution includes its lower-limit term:

\[
(P\otimes h)(x)=C_F\left[2\int_x^1\frac{h(x/z)-h(x)}{1-z}\,dz+
\left(\frac32+2\ln(1-x)\right)h(x)\right].
\]

The matrix integrates a piecewise-linear function of x analytically on each cell after setting `y=x/z`. The first cell uses the exact cancellation of `h(y)-h(x)` against `y-x`. The grid must increase from positive xmin to exactly one; `h(1)=0` is imposed, and no extrapolation below xmin is performed. The physical distribution may have either sign.

For the selected one-loop coupling with fixed integer nf, `U=exp(tau*P)` and `tau=-2*log(alpha_s(Q2)/alpha0)/beta0`. This is not the two-loop coupling or the complete CSS evolution used by KPSY. Threshold matching, Sudakov factors, coefficient functions and the historical fit remain open work.

```julia
using CollinsFit
grid = vcat(exp.(range(log(1e-6), log(0.1); length=200)),
            collect(range(0.101, 1; length=400)))
operator = collinear_lo_operator(grid, 100.0, 4.0, 0.3; nf=3, quantity=:transversity)
h0 = grid .* (1 .- grid).^3  # analytic test input, not an extracted distribution
hQ = apply_operator(operator, h0)
```

This example demonstrates the interface; its grid has no automatic accuracy certification. A production grid requires checks for its actual functions and scale range.

## Independent checks

Run the complete Julia validation job with

```bash
julia --project=numerics numerics/bin/check_collinear_lo.jl runs/collinear-lo/validation.json
```

The report records source/contract hashes, settings, all refinement results, checks and its `reconstructed_theory` scope. Changes during a job fail its source-stability check. Any failed check produces a non-PASS report and a nonzero exit code. The report cannot authorize production execution.

The analytic polynomial input `h0=x*(1-x)^3` has known moments `6/((N+1)*(N+2)*(N+3)*(N+4))`. Its inverse Mellin reference uses the complex harmonic number `digamma(N+1)+gamma` and adaptive integration on two deformed contours, independent of the x-space cell weights. A stable reflection formula avoids overflow in complex digamma at the contour tail. The initial-scale inversion must reproduce the polynomial. The comparison uses [the LO moment solution](https://arxiv.org/html/2409.00221v3#S2.E25), not a legacy pickle or a Python numerical callback.

The fixed target is `1e-4` for fine-grid moment relative errors and scaled point errors. The scale floor is 1% of the largest reference value. On the Linux development run:

| Nodes | Maximum moment relative error, N=1,2,3 | Maximum scaled x-space/Mellin point error |
|---:|---:|---:|
| 201 | 8.23e-4 | 4.93e-3 |
| 401 | 2.25e-4 | 1.58e-3 |
| 801 | 6.11e-5 | 4.80e-4 |
| 1601 | 1.64e-5 | 1.41e-4 |
| 3201 | 4.40e-6 | 4.05e-5 |

The 801- and 1601-node point comparisons do not meet the fixed target. The same mesh is refined to 3201 nodes; the threshold is unchanged. The tested scale interval is Q2=4 to 100 GeV2, alpha0=0.3, nf=3, with seven points from x=0.005 to 0.95. These results do not certify other endpoint shapes, scale intervals, schemes or heavy-flavor matching. Build timings in the report are informational and do not pass G2PERF.

Further checks cover an exact affine-function convolution, identity and scale composition, forward/backward composition, parameter Jacobians, two adjacent finite-difference step sizes, the constrained adjoint and an independent parameter HVP. The differentiated coordinates belong to the initial distribution. Grid, coupling, scales and nf remain fixed; attempts to differentiate those settings are rejected. A fit that varies them needs a separately validated derivative of U.
