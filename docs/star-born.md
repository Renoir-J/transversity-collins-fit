# STAR Born channel and measure contract

`STAR-BORN-001` implements the massless SU(3) observed-parton Born hard factors and complete flavor contractions. It is reconstructed theory. The coupled evolved STAR TMD observable, experimental admission and production fit remain incomplete.

The unpolarized normalization is `d sigma/dt = pi*alpha_s^2/s^2*H_U`. The [PDG 2025 formulae](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-cross-section-formulae.pdf), Eqs.51.4-51.12 on pages1-2, therefore give `H_U=4*s/alpha_s^2*d sigma/dOmega`. The outgoing spin difference is `H_T=T_abc*H_U`, where the independently published [Collins-Heppelmann-Ladinsky Eqs.21-26](https://arxiv.org/html/hep-ph/9305309v1) fix all six nonzero transfer channels and the zero channels. The first incoming parton is polarized; `t` is defined relative to the tagged outgoing parton. Include no extra spin factor2, lepton depolarization factor or identical-observed-particle factor1/2.

`JetBornPlan` owns immutable tuples of PDG flavor IDs and fixed hard terms. The caller explicitly chooses `fragmentation=:quark_only` or `:quark_and_gluon`. Quark-only fragmentation keeps incoming gluon channels. It is the declared approximation used by the [2017 numerical study](https://arxiv.org/html/1707.00913v2); it cannot silently become a justified approximation for every updated STAR bin. The full-flavor Born option supplies no gluon TMD model by itself.

```julia
kin = JetBornKinematics(200., 20., 0.4, 0.4)
plan = JetBornPlan([21,2,-2,1,-1,3,-3], kin; fragmentation=:quark_only)
# fA/fB are unweighted PDFs; D and C are supplied transverse fragmentation responses.
value = jet_born_structure(plan, fA, fB, D, hA, C; alpha_s=0.2)
amplitude = jet_born_asymmetry(value; phase=:sin_phiS_minus_phiH)
```

After the on-shell delta function is integrated over `x_b`, `x_b=-x_a*T/(x_a*S+U)`. The common weight is `alpha_s^2/[S*x_a*x_b*(x_a*S+U)]`, with `S=sqrt_s^2`, `T=-sqrt_s*pT*exp(-eta)` and `U=-sqrt_s*pT*exp(eta)`. Quadrature nodes lie inside the open physical interval. `jet_born_integrals` multiplies this factor and the explicitly sourced coordinate/acceptance weights into both yields before taking their ratio. It supplies no quadrature, acceptance reconstruction or bin-edge inference. Renormalization-scale and coordinate-measure factors that cancel in a fixed-pT ratio generally do not cancel across a pT integral.

The supplied polarized fragmentation response already equals `C=jperp/(z*M_h)*H1perp` in the Trento convention. Do not multiply that factor twice. The moment entering the paper's `J1` Fourier transform is `M_h*H1perp(1)=-Hhat3/(2*z)`, which `jet_collins_fourier_moment` converts from the existing dimensionful KPSY family. A Gaussian control makes the normalization explicit: `C=-Hhat3*jperp*exp(-jperp^2/w)/(pi*w^2)`; its complete transverse integral is `-sqrt(pi)*Hhat3/(2*sqrt(w))`. Differential zero momentum and an integrated observable are different requests.

The pinned JAM Gaussian code and FeynFacet plain Born exports are additional source comparisons. No legacy pickle is read. FeynFacet's WSL readiness check returned `NOT_READY` because no Wolfram kernel was found; no fresh FeynCalc, FORM or Tracer calculation is claimed. Inclusive NLO transfer coefficients alone do not define a finite-R hadron-in-jet Collins observable.

Existing leptonic predictions, likelihood, data cuts/overlap and numerical thresholds are unchanged. Synthetic response/derivative controls, reconstructed hard factors and production results must remain separately labeled. The synthetic-only complete-job JSON executor is not extended by this module.

Run the complete source-reference and numerical check on Linux:

```sh
python scripts/prepare_star_born_references.py
julia --project=numerics numerics/bin/check_jet_born.jl \
  .runtime/star-born-references/reference-data.json runs/star-born-validation.json
```

The preparer downloads three pinned text assets and verifies their Git-blob SHA256 values. An existing pinned external checkout can instead be supplied with `--external-root /path/to/external`. It parses syntax without importing upstream programs. The Julia reference evaluator admits only arithmetic and literal array-index expressions. Algebra checks interpret source numeric literals as exact rationals, rather than claiming execution of the legacy floating-point implementation. A fresh output directory is required; no previous reference preparation is overwritten.

The checks compare all 29 released unpolarized channels and both transfer components of seven released Born classes, then relabel them across the complete five-flavor basis. Seven-flavor one-hot contractions also check the legacy density expressions against both explicit fragmentation choices. Full 13-column gradients and directional HVP use an interior synthetic initial-family point and 256-bit fourth-order finite differences at both predeclared steps. The initial-family functions and analytic inputs are shared, so these checks do not provide independent evolution evidence.

The delta-function measure is checked independently in recoil rapidity. Yields are normalized by `S^2/alpha_s^2` before the relative-floor comparison; all declared quadrature orders remain in the report, with the final order required to pass. Gaussian `J1` transforms test the moment, sign and z factors separately. None of these reports grants production admission.
