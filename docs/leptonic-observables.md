# Reconstructed leptonic observables

`config/leptonic-observable-contract.json` defines `KPSY-LEPTONIC-OBS-002` on top of the [polarized structure functions](polarized-processes.md). The implementation produces fixed-kinematics SIDIS ratios and hadron-2-frame SIA linearized double-ratio amplitudes, including finite transverse-momentum integration and explicit angular source types. Every context and report remains `reconstructed_theory`, `production_ready=false` and `physical_extraction=false`. No experimental data are admitted by this interface.

The source is [KPSY, arXiv:1505.05589v1](https://arxiv.org/html/1505.05589v1), Eqs. 3–5 and 76, 81–92. Publication angle conventions, acceptance and normalization still need dataset-specific review.

## Definitions and fixed inputs

For SIDIS, `sidis_observable` requires an explicit `definition`:

- `:structure_function_ratio` returns `FUT/FUU` with the depolarization factor removed.
- `:sin_phi_sum_amplitude` returns `DNN(y)*FUT/FUU`, where `DNN=2(1-y)/(1+(1-y)^2)`, and requires fixed `0<y<=1`.

The angle convention is `sin(phi_h+phi_s)`. No convention is inferred from an experiment's name. At fixed other kinematics the common cross-section prefactor cancels. Integrating over varying x, z, Q² or y generally requires cross-section weights and experimental acceptance; this interface does not supply that calculation. A prediction at a reported mean and an integral over a bin are distinct prescriptions.

The second positional argument `SIDISAsymmetryConvention` makes a source's phase and DNN normalization explicit under the child contract [`KPSY-SIDIS-CONVENTION-001`](../config/sidis-convention-contract.json). With `R=FUT/FUU`, the supported predictions are:

| Required phase | Required normalization | Prediction |
| --- | --- | --- |
| `:sin_phi_sum` | `:DNN_removed` | `R` |
| `:sin_phi_sum` | `:fixed_y_DNN_included` | `DNN(y)*R` |
| `:sin_phi_sum_minus_pi` | `:DNN_removed` | `-R` |
| `:sin_phi_sum_minus_pi` | `:fixed_y_DNN_included` | `-DNN(y)*R` |

```julia
convention = SIDISAsymmetryConvention(
    phase=:sin_phi_sum_minus_pi, normalization=:DNN_removed,
    source="explicitly reviewed same-axis, spin/dilution/DNN-normalized source")
observable = sidis_observable(context, convention)
```

The [COMPASS 2008 primary paper](https://arxiv.org/abs/0802.2160v2) defines `Phi_C=phi_h+phi_S-pi` and `epsilon_C=f*P_T*DNN*A_Coll`. With the same virtual-photon/lepton-plane azimuth axes, `sin(sum-pi)=-sin(sum)` gives `A_Coll=-R`. The [COMPASS 2022 final journal article](https://doi.org/10.1103/PhysRevLett.133.101903) uses that phase and says its extracted amplitudes were divided by `mean(DNN*f*P_T)`. The source selection is explicit: the interface does not decide that a dataset matches this prescription. This final-observable phase conversion does not redefine transversity, repeat the KPSY/Trento Collins conversion or alter the existing J1 minus sign and mass/z factors.

These amplitudes must already be normalized for target spin magnitude and dilution; raw count modulations are unsupported. DNN-removed mode rejects a supplied `y` to prevent a repeated correction. DNN-included mode requires fixed `0<y<=1` and never divides by DNN, including the zero at `y=1`. A reported mean may be supplied only as an explicit fixed-point approximation: `DNN(mean(y))` generally differs from `mean(DNN)`, and a mean of `DNN*f*P_T` need not factor into means. Source strings, phase and normalization are recorded and guarded against mutation. No species mixture, nuclear model, experimental acceptance or likelihood admission follows.

For hadron-2-frame SIA, `sia_observable` requires `ratio=:UL` or `:UC` and an explicit source type. `SIAAngularFactor(value;source)` retains the supplied fixed ratio of angular means `<sin²θ>/<1+cos²θ>`. `SIAMeanOfRatio(value;source,averaging_measure)` records the supplied mean of pointwise ratios `<sin²θ/(1+cos²θ)>` and requires a nonblank description of the averaging measure. Both accept a finite fixed scalar in [0,1]; neither fits the angular input or reconstructs missing component means. The measure description records caller provenance and does not certify an experiment's weights. The chosen type and measure survive in observable provenance and its mutation guard.

The [Belle tables](https://arxiv.org/abs/0805.2975v1) and [BaBar discussion/tables](https://arxiv.org/abs/1309.5278v2) label a mean of the angular ratio. Belle's [2012 erratum](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.86.039905) also changes asymmetry/error tables, so a matching data adapter must retain the version and avoid a repeated correction. The new source type supplies the representation; it does not admit or normalize those data automatically.

The averaging measure matters. Write `s=sin²θ`, `c=cos²θ`. Under factorized fixed other kinematics, an acceptance-weighted ratio `sum(a*s)/sum(a*(1+c))` equals a mean of `s/(1+c)` with event weights proportional to `a*(1+c)`. Under the same arbitrary weights, the two expressions generally differ. The explicitly synthetic two-angle fixture checks both facts. A general experimental bin can correlate angle, momentum, energy fractions and efficiency; its joint numerator/denominator integrals cannot be replaced by a single scalar without a reviewed approximation or factorization.

The charge categories sum structure functions before division: U is (+,−)+(−,+), L is (+,+)+(−,−), and C=U+L. The amplitudes are `angular*(CollinsU/UUU-CollinsL/UUL)` and `angular*(CollinsU/UUU-CollinsC/UUC)`. These are the paper's linearized double-ratio amplitudes; an exact nonlinear ratio has different Fourier coefficients. Thrust-frame A12 is unsupported.

## Finite momentum intervals

`MomentumRule(Pmin,Pmax;order=32)` supplies an explicitly finite nonnegative interval with the `P*dP` Jacobian. It projects the fixed momentum quadrature onto J0/J1/J2 radial weights once. Each numerator and denominator is integrated separately before normalization and before category differences. Averaging pointwise asymmetries would give a different observable.

All other kinematics and the angular factor stay fixed. Uniform P acceptance is an explicit reconstruction assumption, not an experimental correction. The upper P bound is supplied by the caller; neither an infinite integral nor admission to the low-transverse-momentum approximation is inferred. Signed integrands and numerators are retained. Final unpolarized denominators must be finite and positive; there is no absolute-value clamp or imposed asymmetry cap.

Fixed settings, projected weights, input context and provenance are guarded against mutation. At each parameter value the Julia evaluator prepares the transient evolved initial columns once for all SIA charge categories. These responses are not cached across parameter values. No native evolution, provider, Python or language-model callback runs inside parameter evaluation. The original structure-function interface continues to use the same equations.

## Independent numerical checks

Three complete Julia producers are supplied; Python only prepares locked dependencies and schedules jobs:

```sh
julia --project=numerics numerics/bin/check_observables.jl .runtime/hoppet/build.json runs/observables/analytic.json
julia --project=numerics numerics/bin/check_observable_inputs.jl .runtime/hoppet/build.json .runtime/parton-inputs/build.json runs/observables/inputs.json
julia --project=numerics numerics/bin/check_sidis_conventions.jl .runtime/hoppet/build.json .runtime/parton-inputs/build.json runs/observables/sidis-source.json
```

The first uses explicit polynomial initial inputs and analytic unpolarized bases. Independent Mellin inversion, coupling-time integration and adaptive radial/P integrals check both original SIDIS definitions and both SIA ratios at a point and over three finite P intervals. Independent amplitude/gc moments check gradients and Hessians. Gaussian Hankel identities separately check the projected integration measure. Analytic fixtures are never substituted for physical source inputs.

The second uses the pinned CT10, DSSV08 and DSS14 sources, strict initial samples, and explicitly labeled reference extrapolation at unsupported mu_b scales. It compares the projected results with independent adaptive P integration of the original point structure functions followed by manual charge sums and ratios. This is an independent composition check at shared fixed rows; it does not independently certify native evolution or the source extrapolation policy. Parent matching, Mellin and input validation remain necessary.

Both producers retain joint matching/radial/native refinement levels, finite radial-bound extensions, and the momentum ladder. Coarse native settings are not accepted even when their adjacent difference is small: at least the third level is required. Physical values retain the parent `1e-4` scaled target. The P ladder is 8/16/32/64/128, requires at least order32 and an adjacent difference <=`1e-7`; exhaustion is failure. Full-grid physical comparisons use order128. The independent Gaussian integration target is `1e-10`.

Both angular source types are evaluated at every SIA refinement level. The mean-of-ratio path additionally records its point/three-bin reference comparisons, momentum ladder, source metadata and full 13-parameter inward JVP/HVP arrays. The synthetic event weights and paired acceptance weights give the same effective factor 0.6, so both typed prediction paths share the independent scalar/adaptive finite-difference arrays while being evaluated separately. This reuse is explicit in the report; it is not a second independent evolution solution. Numerical targets and original coarse grids are unchanged.

For both definitions in each process, full 13-parameter inward JVP/HVP checks use a separate scalar initial/charge/CSS assembly, adaptive P integration and generic 256-bit row arithmetic. The native rows are shared fixed evidence: this derivative reference is not an independent evolution solution. The analytic coarse radial grid and the full source-knot-aligned physical coarse grid are retained. Both finite-difference steps (1e-4, 1e-5) must satisfy the unchanged `1e-7` derivative target. Reference values and derivative arrays are exported so audits can recompute the reported errors.

Reports bind the producer, implementation, contracts, reference helpers, lock files and native build receipt to source hashes. CI builds locked native dependencies and preserves reports and native test logs, including failures. Unit fixtures are mathematical/synthetic evidence; the scientific producers report their reconstructed scope separately. Timing pilots compare identical point-output pairs but leave `G2PERF_pass=false`; row payload is not peak RSS and no full-fit performance claim follows.

The third producer checks the source convention separately. A 256-bit adaptive Fourier projection starts with explicitly synthetic positive counts, two signed rational asymmetries and three spin azimuths; it projects onto each declared harmonic independently of the implementation's factor. The pinned-input composition check then evaluates all four conventions at a point and over a finite P interval, comparing separate scalar/adaptive-P values, all 13 gradients and inward JVP/HVP references at both existing finite-difference steps. Its native rows are shared fixed evidence, and its coarse composition context does not replace the parent refinement reports. It retains unshifted-phase and repeated-DNN counterfactuals and leaves production/extraction/performance false.

## Remaining scientific admission

The child [fixed-weight bin interface](leptonic-bins.md) now composes prepared nodes across supplied kinematics by summing angular-weighted numerators and unpolarized yields before division. It requires explicit quadrature, cross-section and acceptance weights. It does not recover experimental acceptance or supply an admitted data adapter.

Historical 262-point reproduction, reviewed publication conventions, multidimensional bins/acceptance, covariance and overlap choices, nuclear prescriptions, STAR, production cache/gradient/performance certification and uncertainty coverage remain open. Original author-code and homogeneous-Collins-evolution limitations remain. This component cannot authorize a real-data extraction or close those requirements.
