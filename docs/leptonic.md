# SIDIS and SIA reconstruction components

The [KPSY-LEPTONIC-001 contract](../config/leptonic-contract.json) implements process-dependent matching, Fourier-Bessel transforms and charge-weighted structure-function assembly. It is reconstructed theory, with production_ready=false. It does not reproduce the historical fit or deliver a new extraction.

The source is [KPSY15](https://arxiv.org/html/1505.05589v1), equations 37–51, 58, 67–74 and 101–109. The CSS coupling/Sudakov and pinned input contracts remain separate. Raw inputs, build files and numerical reports stay outside Git.

## Matching and conventions

`matching_rule(process, channel, x, alpha; order, knots)` returns fixed coordinates and quark/gluon weights. `match_leg` applies them to generic Julia values, retaining derivatives of fitted distributions. Processes are `:SIDIS` and `:SIA`; channels are `:pdf`, `:transversity`, `:ff`, `:collins`. SIA admits fragmentation channels only. Inputs are unweighted f, h1, D or dimensionful KPSY Hhat3.

Convolutions use dy/y, integrated at fixed Gauss nodes in log(y). Optional input knots split interpolation cells. The absorbed hard delta term is −2CF per SIDIS leg and CF(π²−8)/4 per SIA leg. The log(v) multiplying a plus distribution makes it regular because log(1)=0; its delta contribution vanishes. The endpoint log(v)/(1−v) is evaluated by its −1 limit. Gluon values are required where the corresponding nonzero matching term exists; transversity/Collins legs reject gluon input.

Two first-order matched legs are multiplied as written in the structure-function formula. This retains their higher-order product; it does not establish how the missing historical fitter handled this choice. Matching alpha_s and query coordinates are fixed settings. Fitted shape values and their derivatives remain in Julia.

## Radial integrals and flavor assembly

`RadialRule(edges; order)` integrates an explicitly supplied finite interval [0,B]. Nodes exclude b=0. No tail is silently certified negligible. `sidis_fourier` and `sia_fourier` accept amplitudes that already contain CSS damping. Their prefactors are:

| Channel | Harmonic and radial power | Prefactor outside db/(2π) |
|---|---|---|
| SIDIS UU | b J0(Pb/z) | 1/z² |
| SIDIS Collins | b² J1(Pb/z) | −1/(2z³) |
| SIA UU | b J0(Pb/z1) | 1/z1² |
| SIA Collins | b³ J2(Pb/z1) | 1/(4z1³z2) |

P is the momentum of hadron 1 in the hadron-2 frame for SIA. The Collins inputs use Hhat3=−2zM_h H1perp(1)_Trento; adding another mass factor to these transforms would duplicate the conversion. Signed results are retained.

`leptonic_plan` caches the fixed kinematics, Fourier weights and CSS exponents. `leptonic_structure` consumes four matrices, each with explicit PDG flavors as rows and b nodes as columns. Flavor lists contain unique quark/antiquark pairs. SIDIS multiplies same-flavor legs; SIA pairs q in the first hadron with anti-q in the second and sums over both quarks and antiquarks. Gluons enter matching, not the electromagnetic flavor sum. The perturbative exponent and channel-specific nonperturbative damping each appear once. Width dependence is evaluated as exp(−fixed−(gh−gc)*shape), avoiding an underflow-times-overflow rearrangement near gc=gh.

Gauss roots/weights are prepared once at 192 bits and rounded to Float64, with symmetric roots and polynomial-derived weights. Input sampling and transforms use Float64. Cache callers receive copies. Preparations are serialized; the tested workers use one thread. The initial eigenvector-weight rule failed a strict Gaussian cancellation check. Its report is retained locally. Accurate roots/weights and panel refinement satisfy the unchanged target; a coarse mesh can still fail and is reported as such.

## Pinned-input unpolarized pilot

`prepare_unpolarized_sidis` and `prepare_unpolarized_sia` sample locked CT10/DSS14 and prepare fixed matching matrices outside the fitted-parameter loop. CT10 queries are batched with bounded temporary arrays. DSS14 queries obtain all flavors from one interpolation at each coordinate. These methods return plans, two matched matrices and source/range metadata. They do not infer polarized distributions or an experimental likelihood.

Default strict preparation rejects the unsupported μ_b range. A caller can explicitly choose `:reference`, consistently with the CT10 context; the result records out-of-range node counts and remains production_ready=false. Both small-b high scales and large-b low scales occur. This is an exploratory reconstruction convention requiring historical comparison and sensitivity checks.

On Linux, after preparing the locked input build:

```sh
julia --project=numerics numerics/bin/check_leptonic.jl
julia --project=numerics numerics/bin/check_leptonic_inputs.jl
```

The first compares coefficient moments with closed harmonic/trigamma expressions, convolutions with 256-bit quadrature in a different variable, and all Fourier prefactors with exact Gaussian transforms and angular Fourier projections. It checks flavor/CSS assembly against independent adaptive integration, radial bounds/nodes, fitted-value gradients and HVPs. The Gaussian refinement table preserves coarse failures. The second exercises pinned-input unpolarized SIDIS/SIA values, strict rejection, reference-policy labeling and joint matching/radial refinement. Reports bind source and build hashes. The two checks run in the ordinary and locked-input Linux CI jobs respectively.

Fitted initial transversity/Collins parameterizations, scalable evolution integration, historical 262-point reproduction, publication-specific angular/depolarization factors and acceptance/bin treatment remain required. In particular these structure functions do not themselves implement the Belle thrust-frame observable, SIA double ratios, STAR jets or a reviewed experimental likelihood. Physical gradient/performance admission and calibrated interval coverage remain open. A component PASS does not authorize production.
