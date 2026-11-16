# Conditional historical BaBar integrals

`KPSY-HISTORICAL-SIA-INTEGRAL-001` was frozen in `738bbb6` before this implementation. The adapter connects 72 BaBar RF0 z-bin amplitudes, 36 UC and 36 UL, to separate finite-momentum integrals at the unchanged rounded KPSY parameters. Complete real-source numerical/audit/cache validation remains pending. This is reconstructed theory with historical candidate rows, not an accepted likelihood, fitted optimum or production extraction.

The existing intake checks 792 fields against the pinned [BaBar v2 source](https://arxiv.org/abs/1309.5278v2), `paper.tex`, `tab:totZuc` and `tab:totZul`. These include reported P means. The means remain source metadata: an integrated amplitude cannot be replaced by a prediction at mean P. Both charge-pair yields must be integrated separately before constructing the normalized UL/UC ratios, as in [KPSY Eq.82](https://arxiv.org/html/1505.05589v1#S2.E82).

The declared conditional reconstruction uses `0 <= qT <= 3.5 GeV`, mapped to P with the fixed reported z1 mean and retaining `P*dP=z1^2*qT*dqT` in absolute yields. The upper cut is printed in the primary selection. Zero is an explicitly declared radial origin, not a recovered lower experimental acceptance. Uniform momentum acceptance, fixed z means and factorization of the primary mean-of-angular-ratio scalar remain model assumptions. The source-literal Q2=112.36 remains conditional. No source value, error component, angular definition or selection changes.

The [contract](../config/historical-sia-integral-contract.json) requires separate momentum refinement, native/matching/radial refinement and finite tails, plus independent adaptive-qT scalar assembly. It preserves the existing 1e-4 physical and 1e-7 momentum/derivative targets. All coarse misses remain visible. Representative all13 gradient/HVP checks are declared at three source geometries with both ratios; their scope is explicit. The one-sided fourth-order reference avoids crossing the rounded Table1 boundaries and must first pass independent mathematical controls.

The 36 exact geometries contain 14 distinct z queries. This implementation uses unchanged separate-geometry preparation and adds no query or parameter-response cache. A future fixed-query cache requires exact source/settings identity, unchanged separate-query comparators and mutation checks. The inventory alone demonstrates no runtime speedup or production throughput. Native state remains isolated, with all native/provider work outside parameter evaluation.

## Whole jobs

```sh
snakemake --snakefile workflow/historical-sia-integrals.smk --cores 2 \
  --config run_dir=runs/babar-integral-example timeout_seconds=18000
python scripts/check_cached_historical.py config/historical-sources.json \
  data/raw/historical runs/babar-integral-example/historical
python scripts/check_cached_historical_sia_integrals.py \
  runs/babar-integral-example/request.json runs/babar-integral-example/diagnostic
```

Current native/parton receipts and a verified raw cache can be supplied with the usual workflow settings. Each Julia worker owns a complete prediction or audit job. Python schedules processes and renders completed values; no Python, provider, native or language-model callback enters parameter evaluation.

Each native/radial preparation first refines momentum over orders16/32/64/128/256, with minimum64 and target1e-7 for all four pair UU/Collins integrals, category sums and amplitudes. At least three native/matching/radial levels must satisfy1e-4, followed by finite bounds12/24/32 and the same1e-4 extension target. Exhausted budgets remain FAIL. Independent scalar charge/CSS/Bessel assembly uses adaptive qT integration and shares prepared rows/initial functions, so it is not independent evolution.

The three prospectively selected low/middle/high diagonal geometries require both ratios and all13 gradient/HVP checks. The 256-bit fourth-order inward reference retains five samples at each of two fixed steps. Its closed-function controls are explicitly synthetic. This is representative evidence, not an all72-geometry or production-derivative certificate.

The separate reader verifies current source/request/contract/build bindings and recomputes retained refinement, category, ratio, derivative and stencil arithmetic. It reevaluates closed mathematical controls but does not rerun native evolution/predictions, physical adaptive integrals, AD, Fourier quadrature or raw primary parsing. Object reuse and source stability remain labeled producer attestations; cold timings and row payload sizes remain informational.

## Retention and products

Fresh request/output paths are required. The scheduler declares only a completion marker; its final target always replays source preparation and the numerical audit. Initial FAIL invocations and atomic in-progress reports retain existing evidence through failures. The shared time budget applies to the complete Julia prediction/audit jobs. Interruptions remain nonzero even when a cooperative child exits0. Explicitly synthetic scheduler controls cover ordinary failure, resource timeout and parent SIGTERM, plus forged-marker rejection. They execute no numerical predictions and cannot approve the integral path. CI reserves time to retain available reports/logs after an interrupted numerical step.

Fifteen cache controls reject altered windows, Jacobian/angle/source provenance, errors/covariance, point substitution, retained predictions/derivative samples, table/figure and production scope, then restore original bytes and timestamps. A completed diagnostic and audit produce an ordered CSV and discrete UL/UC panels by the primary z1 bins. Statistical bars and reported systematic fields remain separate. There is no chi2, pull, combined covariance, normalization nuisance or uncertainty interval. A source review, synthetic closure or completed P-point diagnostic cannot approve this integral path.

The original BaBar P-point refinement misses and paper-center/coupling failures remain separate. Missing author code, covariance/overlap, acceptance, corrected Belle primary bytes, production gates and the complete B0-B4/V1/uncertainty/publication scope remain open.
