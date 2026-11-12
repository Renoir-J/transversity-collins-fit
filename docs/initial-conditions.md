# Initial transversity and Collins conditions

[KPSY-INITIAL-001](../config/initial-contract.json) reconstructs the 13-parameter initial family from [KPSY15 equations 110–115 and Table 1](https://arxiv.org/html/1505.05589v1). It is reconstructed theory, production_ready=false. The rounded printed parameter vector is a source reference, not an optimizer result from this repository and not a full historical fit state.

## Coordinates and boundary behavior

`KPSYParameters(theta; gh=0.042)` accepts direct physical coordinates in the order recorded in the contract. Nh lies in the closed interval [-1,1]. Shape exponents are nonnegative in this declared reconstruction domain; the original MINUIT bounds were not recovered. Nc is a finite dimensionful amplitude in GeV, with no automatic unit-interval or positivity constraint. gc must be below the supplied fixed gh. An exp/tanh-only decoder cannot represent the printed Nd=-1 and beta=0 values exactly.

Transversity is Nh times a peak-normalized x^a(1−x)^b shape times the fixed Soffer input (f1+g1)/2. Logarithms avoid unstable products. The value limit t log(t)=0 at t=0 permits a or b to be zero for point evaluation; their normalization derivatives there are singular and are rejected. The printed a,b are positive. Negative Soffer input is rejected without clipping. The historical family has zero sea/heavy transversity.

Favored and unfavored Hhat3 use separate pi+/u and pi+/d bases. Antiquark light flavors follow the paper's common favored/unfavored functions, rather than substituting a separate antiquark FF. Strange and antistrange retain separate pi+ FF bases with the unfavored parameters. Zero heavy Collins and the neutral-pion charge average are explicit reconstruction assumptions. They do not certify experimental charm/background treatment. The Trento conversion requires an explicit positive hadron mass: H1perp(1)=−Hhat3/(2zM).

At x/z=1, distributions and their parameter derivatives are exactly zero before evaluating endpoint logs or powers. This includes beta=0 and nested differentiation, avoiding 0 times infinity. Gradients at Nh and beta boundaries are checked against the inward derivative of the formula. Simultaneously reversing the two Nh and two Nc amplitudes leaves leptonic products unchanged before sign-sensitive external constraints. The paper chooses positive u; lattice-constrained fits must retain their own sign information.

## Inputs and evolution

`prepare_initial_conditions` samples locked CT10, DSSV08 and DSS14 once at Q0²=2.4. No provider calls enter parameter evaluation. It requires strict source support: x≥1e-5, z≥0.05, and grids ending at one. Unweighted u,d Soffer inputs and four independent pion FF bases are retained. No endpoint below the available x support is invented for a full tensor charge.

`initial_source_grids` creates logarithmic/sqrt backgrounds augmented by the locked source knots and declared prediction probes. Distinct Fortran literal values remain distinct. A background-only pilot failed the frozen 1e-4 point-refinement target, with differences 1.22e-4 for h1 and 4.80e-4 for Hhat3. Aligning source/probe knots passed that unchanged target; analytic affine-generator checks also verify the resulting nonuniform grids. The failed report remains local.

`initial_distributions` returns explicit flavor-by-grid matrices. `evolve_initial_conditions` requires the matching grids, Q0, separate transversity/Hhat3 quantity tags and the contracted asymptotic two-loop, nf=5, Lambda=0.225 settings. Generic Julia matrix multiplication retains parameter derivatives. Other coupling definitions or mismatched scales are rejected. This connects the family to the validated dense evolution reference; it does not implement an efficient many-b backend.

## Verification and remaining work

```sh
julia --project=numerics numerics/bin/check_initial.jl
julia --project=numerics numerics/bin/check_initial_inputs.jl
```

The analytic check uses independent 256-bit direct-power expressions, dictionary-based flavor relations, closed analytic gradients/Hessians, inward beta=0 differences and independently constructed initial Jacobians acted on by the fixed evolution operator. Unit tests include boundaries, endpoints, sign symmetry, matrix applications and rejected setting mismatches.

The pinned-input pilot refines background sizes 801/1601/3201 with all source/probe knots, checks Soffer inputs, evolved endpoint values, charged/neutral relations and finite boundary derivatives, and compares evolved distributions at Q²=10. Reports bind code, contracts and build hashes. It records dense operator bytes and informational preparation time; the many-scale memory figure is a hypothetical dense-cache estimate, not measured fitting performance or a speedup claim.

Scalable many-scale evolution and matching integration, physical gradient/performance admission, target/axis/angular/acceptance adapters, the historical 262-point reproduction, covariance/overlap and interval coverage remain required. The original author-modified HOPPET/fitter remains unlocated. These component checks do not authorize production fits.
