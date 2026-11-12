# Reconstructed leptonic observables

`config/leptonic-observable-contract.json` defines `KPSY-LEPTONIC-OBS-001` on top of the [polarized structure functions](polarized-processes.md). The implementation produces fixed-kinematics SIDIS ratios and hadron-2-frame SIA linearized double-ratio amplitudes, including finite transverse-momentum integration. Every context and report remains `reconstructed_theory`, `production_ready=false` and `physical_extraction=false`. No experimental data are admitted by this interface.

The source is [KPSY, arXiv:1505.05589v1](https://arxiv.org/html/1505.05589v1), Eqs. 3–5 and 76, 81–92. Publication angle conventions, acceptance and normalization still need dataset-specific review.

## Definitions and fixed inputs

For SIDIS, `sidis_observable` requires an explicit `definition`:

- `:structure_function_ratio` returns `FUT/FUU` with the depolarization factor removed.
- `:sin_phi_sum_amplitude` returns `DNN(y)*FUT/FUU`, where `DNN=2(1-y)/(1+(1-y)^2)`, and requires fixed `0<y<=1`.

The angle convention is `sin(phi_h+phi_s)`. No convention is inferred from an experiment's name. At fixed other kinematics the common cross-section prefactor cancels. Integrating over varying x, z, Q² or y generally requires cross-section weights and experimental acceptance; this interface does not supply that calculation. A prediction at a reported mean and an integral over a bin are distinct prescriptions.

For hadron-2-frame SIA, `sia_observable` requires `ratio=:UL` or `:UC` and a source-labeled `SIAAngularFactor`. The angular input is the fixed ratio of means `<sin²θ>/<1+cos²θ>`, not the mean of their pointwise ratio. It is not a fitted parameter or an automatically reconstructed experimental average.

The charge categories sum structure functions before division: U is (+,−)+(−,+), L is (+,+)+(−,−), and C=U+L. The amplitudes are `angular*(CollinsU/UUU-CollinsL/UUL)` and `angular*(CollinsU/UUU-CollinsC/UUC)`. These are the paper's linearized double-ratio amplitudes; an exact nonlinear ratio has different Fourier coefficients. Thrust-frame A12 is unsupported.

## Finite momentum intervals

`MomentumRule(Pmin,Pmax;order=32)` supplies an explicitly finite nonnegative interval with the `P*dP` Jacobian. It projects the fixed momentum quadrature onto J0/J1/J2 radial weights once. Each numerator and denominator is integrated separately before normalization and before category differences. Averaging pointwise asymmetries would give a different observable.

All other kinematics and the angular factor stay fixed. Uniform P acceptance is an explicit reconstruction assumption, not an experimental correction. The upper P bound is supplied by the caller; neither an infinite integral nor admission to the low-transverse-momentum approximation is inferred. Signed integrands and numerators are retained. Final unpolarized denominators must be finite and positive; there is no absolute-value clamp or imposed asymmetry cap.

Fixed settings, projected weights, input context and provenance are guarded against mutation. At each parameter value the Julia evaluator prepares the transient evolved initial columns once for all SIA charge categories. These responses are not cached across parameter values. No native evolution, provider, Python or language-model callback runs inside parameter evaluation. The original structure-function interface continues to use the same equations.

## Independent numerical checks

Two complete Julia producers are supplied; Python only prepares locked dependencies and schedules jobs:

```sh
julia --project=numerics numerics/bin/check_observables.jl .runtime/hoppet/build.json runs/observables/analytic.json
julia --project=numerics numerics/bin/check_observable_inputs.jl .runtime/hoppet/build.json .runtime/parton-inputs/build.json runs/observables/inputs.json
```

The first uses explicit polynomial initial inputs and analytic unpolarized bases. Independent Mellin inversion, coupling-time integration and adaptive radial/P integrals check both SIDIS definitions and both SIA ratios at a point and over three finite P intervals. Independent amplitude/gc moments check gradients and Hessians. Gaussian Hankel identities separately check the projected integration measure. Analytic fixtures are never substituted for physical source inputs.

The second uses the pinned CT10, DSSV08 and DSS14 sources, strict initial samples, and explicitly labeled reference extrapolation at unsupported mu_b scales. It compares the projected results with independent adaptive P integration of the original point structure functions followed by manual charge sums and ratios. This is an independent composition check at shared fixed rows; it does not independently certify native evolution or the source extrapolation policy. Parent matching, Mellin and input validation remain necessary.

Both producers retain joint matching/radial/native refinement levels, finite radial-bound extensions, and the momentum ladder. Coarse native settings are not accepted even when their adjacent difference is small: at least the third level is required. Physical values retain the parent `1e-4` scaled target. The P ladder is 8/16/32/64/128, requires at least order32 and an adjacent difference <=`1e-7`; exhaustion is failure. Full-grid physical comparisons use order128. The independent Gaussian integration target is `1e-10`.

For both definitions in each process, full 13-parameter inward JVP/HVP checks use a separate scalar initial/charge/CSS assembly, adaptive P integration and generic 256-bit row arithmetic. The native rows are shared fixed evidence: this derivative reference is not an independent evolution solution. The analytic coarse radial grid and the full source-knot-aligned physical coarse grid are retained. Both finite-difference steps (1e-4, 1e-5) must satisfy the unchanged `1e-7` derivative target. Reference values and derivative arrays are exported so audits can recompute the reported errors.

Reports bind the producer, implementation, contracts, reference helpers, lock files and native build receipt to source hashes. CI builds locked native dependencies and preserves reports and native test logs, including failures. Unit fixtures are mathematical/synthetic evidence; the scientific producers report their reconstructed scope separately. Timing pilots compare identical point-output pairs but leave `G2PERF_pass=false`; row payload is not peak RSS and no full-fit performance claim follows.

## Remaining scientific admission

Historical 262-point reproduction, reviewed publication conventions, multidimensional bins/acceptance, covariance and overlap choices, nuclear prescriptions, STAR, production cache/gradient/performance certification and uncertainty coverage remain open. Original author-code and homogeneous-Collins-evolution limitations remain. This component cannot authorize a real-data extraction or close those requirements.
