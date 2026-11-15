# Explicit SIA qT windows

`SIAQTWindow` and `sia_qt_observable` convert a caller-declared finite `qT=P_hperp/z1` window into the existing `P*dP` integration at fixed theory `z1`. They follow [KPSY Eqs.82/85/86](https://arxiv.org/html/1505.05589v1) and the prospective [KPSY-SIA-QT-WINDOW-001 contract](../config/sia-qt-window-contract.json). The source states a `qT<3.5 GeV` cut for Belle/BaBar. This interface requires both bounds and their provenance explicitly; it does not infer a publication's lower limit or acceptance.

```julia
window = SIAQTWindow(0.0, 3.5;
    source="explicit reconstructed KPSY-inspired qT window",
    lower_limit_source="caller-declared radial origin; experimental acceptance unresolved")
observable = sia_qt_observable(context, window;
    ratio=:UL,
    angular=SIAMeanOfRatio(0.6;
        source="explicitly synthetic event angular average",
        averaging_measure="fixed synthetic event weights"),
    order=64)
value = predict_observable(observable, theta)
```

`context` must already be a prepared hadron-2-frame SIA plan. Its fixed `z1` gives `Pmin=z1*qTmin` and `Pmax=z1*qTmax`. The existing weights retain `P*dP=z1^2*qT*dqT`. Each signed polarized numerator and unpolarized denominator is integrated separately before charge-category sums and linearized UL/UC ratios. The original angular source kind, scalar and averaging measure remain explicit. Window limits, kinematics and quadrature are preparation settings; derivatives with respect to those settings reject.

The standalone complete Julia check uses synthetic fixed response rows and uniform momentum acceptance:

```sh
julia --project=numerics numerics/bin/check_sia_qt_windows.jl runs/sia-qt-window/validation.json
```

It retains adjacent quadrature levels, absolute polynomial/Gaussian integrals with 256-bit closed references, independent adaptive qT weights and separate yields, all thirteen gradient columns and Hessian-vector products, and two finite-difference steps. A missing common Jacobian can cancel from a fixed-z ratio; the absolute-yield counterfactual must still reject. First-six transversity derivatives vanish for SIA. These checks share prepared rows/initial functions and do not independently validate evolution. The CLI refuses to overwrite an earlier report.

No data values, cuts or row selection are changed by preparing this interface. It is reconstructed theory with synthetic validation, not a historical262-point fit, experimental admission or production certificate. A fixed theory mean does not identify `mean(P/z)` with `mean(P)/mean(z)` or recover a full z-bin/event acceptance. Finite-mass energy/momentum fractions, corrected Belle primary bytes, joint angular acceptance, covariance/overlap, the absent author fitter and the retained paper-center failures remain unresolved. Production, physical extraction, historical-fit readiness and G2PERF stay false. PR #39's fitting evidence remains at its original commit.
