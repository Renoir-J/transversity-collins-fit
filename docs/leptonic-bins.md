# Fixed-weight leptonic bin composition

`KPSY-LEPTONIC-BIN-001` composes prepared [leptonic observable nodes](leptonic-observables.md) across supplied kinematic quadrature. Julia owns evaluation and derivatives. All results remain `reconstructed_theory`, with production, extraction and performance admission false.

For SIDIS the result is `sum(w*a*N)/sum(w*D)`. Here `N,D` are polarized/unpolarized structure-function integrals, `a` is each node's explicit phase/DNN factor, and `w` includes the supplied quadrature, acceptance and cross-section prefactor. In the `dx dy dz d2P` coordinates of [KPSY Eq.3](https://arxiv.org/html/1505.05589v1), the varying prefactor is proportional to `[1+(1-y)^2]/(y*Q2)`. Changing integration coordinates needs the corresponding Jacobian. A common overall constant cancels. Averaging point asymmetries or applying `DNN(mean(y))` generally gives a different quantity.

For SIA, supply an N-by-2 matrix of U/L weights. Each includes the unpolarized angular prefactor `(1+cos(theta)^2)/Q2`, the declared phase-space measure and any conditional category acceptance. The node's angular factor belongs inside the polarized integral. Sum U and L separately, then form `NU/DU-NL/DL` or `NU/DU-(NU+NL)/(DU+DL)`, retaining the parent's linearized double-ratio approximation ([KPSY Eqs.80–84](https://arxiv.org/html/1505.05589v1)). C is the sum of those same integrated yields. It is not a third independent measurement.

```julia
# nodes have already been prepared, with explicit source conventions.
bin = leptonic_bin(nodes, weights;
    source="reviewed reconstruction specification",
    measure="declared quadrature coordinates, intervals and units",
    weight_definition="quadrature times Jacobian, sigma0 and fixed acceptance",
    acceptance="explicit supplied reconstruction weights")
prediction = predict_bin(bin, theta)
```

Weights must be fixed, finite and nonnegative, with nonzero support in each category. The bin owns copies of the node vector and weight array. It checks every fixed weight, denominator, metadata field and descendant node against complete typed snapshots. Changes to caller-owned arrays do not change a prepared bin; changes to its fixed evidence are rejected. Parameter-dependent weights and provider, native, Python or language-model callbacks are unsupported inside evaluation.

Nodes must share the model, process, flavor basis, target/species, observable convention and point/P-integral mode. Their fixed kinematics and conditional angular factors may vary. Each finite-P node already contains its `P*dP` measure; supply that Jacobian only once. The interface does not infer the beam relation `y=Q2/(x*S)`, support, cuts, coordinate units, experimental efficiency or an event-level DNN correction. A source-averaged angular factor still needs an independently reviewed conditional averaging prescription. This composition interface alone does not establish multidimensional quadrature convergence or an experimental adapter.

Two complete Julia jobs produce separate analytic-input and pinned-input reports:

```sh
julia --project=numerics numerics/bin/check_leptonic_bins.jl analytic .runtime/hoppet/build.json - runs/bins/analytic.json
julia --project=numerics numerics/bin/check_leptonic_bins.jl inputs .runtime/hoppet/build.json .runtime/parton-inputs/build.json runs/bins/inputs.json
```

They vary x/z/Q/P across explicitly supplied nodes and compare point and finite-P bin compositions with separate scalar initial/charge/CSS loops and adaptive P integration. Full 13-column gradients and inward JVP/HVP checks use 256-bit references at the parent's unchanged steps/targets. The coarse fixed native rows are shared evidence: these are composition checks, not independent evolution or radial-convergence solutions. Pinned initial inputs remain strict and unsupported Fourier-scale queries retain explicit reference extrapolation.

Separate synthetic joint-acceptance fixtures check y/auxiliary-coordinate and theta/auxiliary-coordinate quadrature against closed integrals. Coarse SIDIS quadrature errors stay in the report; only the declared final-order comparison is required to meet the unchanged analytic target. These acceptances are mathematical fixtures, not measured efficiencies. Existing kernel refinement reports remain required. Unit tests also reject normalized-asymmetry averaging, convention mixtures, invalid weights and mutated evidence.

No historical or updated data set is admitted. Author-code recovery, Collins-evolution approximation, source support, physical binning/acceptance, covariance/overlap, production derivatives/performance and calibrated uncertainty coverage remain open.
