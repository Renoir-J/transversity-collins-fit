# Reconstructed polarized leptonic structure functions

[KPSY-POLARIZED-PROCESS-001](../config/polarized-process-contract.json) joins the initial parameter family, retained LO evolution, polarized Wilson matching and CSS/Fourier components. It returns signed `UU` and `Collins` structure functions at fixed kinematics for SIDIS and hadron-2-frame SIA. Conventions follow [KPSY15 Eqs. 67, 68, 73, 74, 105, 108 and 109](https://arxiv.org/html/1505.05589v1). Results remain `reconstructed_theory`, `production_ready=false` and `physical_extraction=false`.

The context prepares strict Q0 samples, matched evolution rows at every actual radial mu_b squared scale, unpolarized matching matrices, flavor maps and fixed CSS arrays once. Parameter evaluation constructs two valence h1 columns and four pi+ Hhat3 columns, applies detached Julia rows and evaluates the signed flavor/Fourier sum. It makes no provider, native, Python or language-model calls. Float64-backed ForwardDiff primal/tangent columns share Julia BLAS products, including nested derivatives; BigFloat-backed derivatives retain generic row arithmetic for comparisons.

Q0 squared=2.4 GeV squared, nf=5 and Lambda=0.225 GeV remain fixed. The 13 physical coordinates and closed boundaries are unchanged. Damping uses gh-gc directly and requires gh-gc>0. Singular normalization derivatives at a=0 or b=0 remain errors. Full twist-three Collins mixing, sea transversity and heavy Collins functions are not added.

## Charges and targets

Unpolarized fragmentation matrices always represent pi+. Pi- conjugates that basis and pi0 averages both charges. SIA pairs the second leg with the opposite PDG flavor once. The polarized four-basis family preserves its parent's charge-independent strange prescription, including unequal s/sbar bases in validation.

`NucleonTarget(:proton)` and `NucleonTarget(:neutron)` use exact u/d isospin interchange, including antiquarks in the unpolarized sum. Other targets require a source and separate fixed proton/neutron weights for the denominator and numerator. These implement impulse-model algebra; test weights are not a reviewed nuclear prescription. No dilution or effective-polarization factor is inferred.

`with_momentum(context, P)` reuses providers, rows, initial samples and CSS arrays, changing Fourier weights and cached UU. `with_pion_charges(context, c1, c2)` reassembles from pi+ bases. SIDIS requires a zero second charge slot. Mutation checks bind the context's fixed unpolarized/CSS arrays, metadata and initial samples. The [fixed-cache snapshot contract](numerical-cache.md) covers all accessible mutable settings and source arrays. Retained row matrices have owned read-only storage and a complete construction-time snapshot; block indices are checked on application. Production serialization and cache admission remain separate work.

## Preparing a context

The default source policy is strict. Actual mu_b scales exceed locked provider support at both ends; strict preparation rejects those queries. Explicit `:reference` policy records extrapolation without production admission. Initial samples remain strict at Q0 under either policy.

```julia
using CollinsFit, JSON
input_build = ".runtime/parton-inputs/build.json"
native_build = ".runtime/hoppet/build.json"
# The constructor verifies the input build; keep its receipt available as well.
pdf = CT10Input(input_build; policy=:reference)
try
    input_root = JSON.parsefile(input_build)["input_root"]
    ff, hel = load_dss14(input_root), load_dssv08(input_root)
    model = KPSYCSS()
    radial = RadialRule(parton_radial_edges(model, ff, input_root, 24.); order=32)
    context = prepare_polarized_sidis(pdf, hel, ff, native_build, model, radial,
                                    2.4, .3, .1, .3, 1;
                                    matching_order=16, intervals=2048,
                                    range_policy=:reference)
    theta = collect(kpsy_published_parameters().values) # rounded source values
    result = polarized_structure(context, theta)       # structure functions
finally
    close(pdf)
end
```

`prepare_polarized_sia` requires the pion FF input, two fractions and charges, and only `frame=:hadron2`. It has no helicity/PDF provider dependency. The common radial pilot helper uses the union of locked PDF/FF scale knots; callers supply an explicit radial rule. A thrust-frame asymmetry cannot be obtained by changing a frame label.

## Validation and limits

```bash
julia --project=numerics numerics/bin/check_hoppet_batch.jl \
  .runtime/hoppet/build.json .runtime/parton-inputs/build.json runs/polarized/batch.json
julia --project=numerics numerics/bin/check_polarized.jl \
  .runtime/hoppet/build.json runs/polarized/analytic.json
julia --project=numerics numerics/bin/check_polarized_inputs.jl \
  .runtime/hoppet/build.json .runtime/parton-inputs/build.json runs/polarized/inputs.json
```

The analytic producer uses polynomial initial inputs, independently inverted evolved/matched Mellin moments, adaptive radial integration, independent coupling-time integration and scalar flavor sums. It checks complete values, charges, amplitude/gc derivatives and Hessians, radial/evolution refinement, finite bounds and all-13-parameter inward JVP/HVP finite differences. These mathematical inputs are not experimental measurements or physical input replacements.

The input producer uses the rounded published family and locked CT10/DSSV08/DSS14 sources. It retains coarse levels, checks matching/radial/evolution refinement at five transverse momenta, range-policy labels, charge reuse, separate target spin weights, radial tails and generic 256-bit finite differences. Resource reports identify the point, radial nodes and coefficient payload; they do not certify peak RSS or full-fit throughput. Synthetic row unit fixtures cannot replace these producers.

No function here computes experimental angular/depolarization factors, SIA double ratios, bin acceptance, STAR jets or an admitted data likelihood. Signed structure functions are not automatically converted to ratios or clamped positive. Historical 262-point reproduction, covariance/overlap, production derivative/performance policy and interval coverage remain required before a real-data fit.
