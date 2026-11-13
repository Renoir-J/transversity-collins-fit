# Coupled STAR reconstruction

`STAR-EVOLVED-C11-001` couples the reviewed Born algebra to collinear incoming PDFs and h1, a single fragmentation TMD, and the existing thirteen-parameter KPSY initial family. It supplies differential points and explicitly weighted bins. This is reconstructed theory under numerical validation. It does not admit experimental data, production fits or publication results.

The [scientific contract](../config/star-evolved-contract.json) declares a canonical Collins-11 reconstruction with `Q=mu=pT`. The original STAR17 program has not been recovered. The printed STAR17 Sudakov measure and half-leg convention do not establish its actual numerical implementation. This reconstruction instead uses the KPSY pair exponent divided by two, consistently with the single-TMD equations in [1705.08443](https://arxiv.org/html/1705.08443v2). It is distinct from finite-R siTMDFJF matching at `mu_J~pT*R`.

Incoming unweighted PDFs and h1 are evaluated at pT. The beam has no incoming TMD or lepton depolarization factor. The outgoing D and canonical Hhat3 are evolved/matched at every actual mu_b node. C11 matching omits SIDIS/SIA hard-factor constants. The Fourier numerator already includes the Trento moment conversion: `-integral b^2 db J1(jT*b/z) Hhat3/(4*pi*z^3)`. Applying another `-1/(2*z)` would double-count the conversion. The inherited homogeneous Collins evolution is an explicit approximation to full twist3 evolution.

Both fragmentation choices must be explicit:

- `quark_only` requires `gluon_model=:none`. Incoming gluons and collinear Dg mixing into quark TMDs remain included.
- `quark_and_gluon` requires `gluon_model=:casimir_common_same_intrinsic_1705`. Gluon matching sums all ten quark/antiquark FFs. Its Sudakov and Casimir-scaled common nonperturbative term follow the stated paper prescription; the same intrinsic gh is an assumption. Gluon Collins remains zero.

The complete massless nf5 PDG basis is retained in the caller's order. Sea/heavy transversity and heavy/gluon Collins remain zero under the inherited initial family. The standard jet axis, radius and scale prescription are recorded. A radius value alone does not supply finite-R corrections.

## Preparation and evaluation

Julia prepares native HOPPET rows, fixed input samples, hard contractions and transforms before parameter evaluation. The beam uses one common log grid, native evolution/interpolation of an identity basis, and detached query rows. Native precision-error messages reject preparation even when the ABI status is zero. This avoids the original per-query preparation's retained precision-error failure without changing the physical domain or numerical acceptance targets. Julia then evaluates generic initial functions, retained row products, parameter-dependent Collins damping and contractions. No provider, Python or language-model callback enters the parameter loop. Fixed contexts validate live snapshots and reject mutated numerical data.

With verified retained HOPPET and parton build receipts:

```julia
model = KPSYCSS()
radial = RadialRule(parton_radial_edges(model, ff, input_root, 24.); order=32)
point = prepare_star_point(pdf, helicity, ff, hoppet_build,
    model, radial, 200., 12.9, .3, .3, .3, 1;
    fragmentation=:quark_only, gluon_model=:none,
    jet_radius=.6, axis=:standard, scale_prescription=:historical_pT_C11,
    outer_order=256, intervals=1024, matching_order=8,
    range_policy=:reference)
value = star_structure(point, theta)  # separate UU and UT
amplitude = predict_star(point, theta; phase=:sin_phiS_minus_phiH)
```

`range_policy=:reference` explicitly permits and records provider extrapolation at radial scales. It cannot admit a fit. Initial Soffer/FF samples always require strict inherited Q0 support. A default strict query can fail for unsupported radial scales; the interface does not substitute analytic inputs.

`star_bin(points, weights; phase, measure_source)` applies supplied jet/acceptance/Jacobian weights to absolute UU and UT yields separately, then forms the ratio. All points must share beam energy and model choices. `star_momentum_bin(point, MomentumRule(low, high); phase, measure, measure_source)` requires `measure=:jT_djT` or `:d2jT`. The latter includes the azimuthal factor 2*pi in both yields. A differential jT=0 point has zero UT and never selects an integrated expression. Experimental bin edges and acceptance are not inferred from horizontal uncertainties.

## Validation and limits

Run the separately scoped Linux producers:

```bash
julia --project=numerics numerics/bin/check_star_evolved.jl analytic \
  .runtime/hoppet/build.json unused runs/star/analytic.json
julia --project=numerics numerics/bin/check_star_evolved.jl inputs \
  .runtime/hoppet/build.json .runtime/parton-inputs/build.json runs/star/inputs.json
```

The analytic producer compares to Mellin evolution, independent single-TMD matching/Sudakov integration, adaptive Fourier integration and recoil-rapidity measure. Its oracle shares the separately reviewed Born coefficients and the declared analytic family. Pinned-input checks independently integrate source FF matching. Full thirteen-parameter gradient and directional-HVP checks use a scalar reference, 256-bit fourth-order differences and two frozen steps. That derivative reference shares retained native rows and provider samples; it is not an independent evolution solver. Coarse refinement levels and failed attempts remain evidence. No acceptance threshold is adjusted to rescue a failed run.

Unit fixtures explicitly use synthetic rows for flavor/phase, gluon-singlet, fixed-bin, nested-AD and mutation tests. These are distinct from reconstructed analytic and pinned-input reports. Neither set establishes physical fitting throughput or calibrated uncertainty coverage. The finite-kernel JSON executor and Snakemake fitting DAG still accept synthetic closure only.

Remaining work includes independent physical benchmarks, measured STAR acceptance and finite-R treatment, updated covariance/overlap admission, review of the Collins evolution approximation, frozen production performance/coverage gates and a production executor. Native FeynFacet toolchain readiness is a separate requirement for a fresh derivation; released coefficient comparisons do not establish it.

A separate common-grid diagnostic shows that near-endpoint h1 interpolation is not uniformly accurate under the old component's relative-floor criterion (at1024 intervals, the largest scaled point error was about0.014 near x=.999, while the native-row/forward comparison agreed). That diagnostic is retained as a limitation, not promoted to PASS. This contract tests the full xa-integrated STAR observables against independent Mellin/recoil integrals; it does not inherit the old query-aligned component's uniform pointwise accuracy claim. No endpoint cuts or clipping conceal the tail. Production still requires endpoint/evolution error qualification across the intended physical workload.
