# Prospective historical BaBar integrals

`KPSY-HISTORICAL-SIA-INTEGRAL-001` specifies the next 72 BaBar RF0 z1/z2 amplitudes, 36 UC and 36 UL, before new numerical candidates. It is a prospective contract: no source adapter, real-source integral producer, audit, figure or numerical PASS is supplied by this document.

The existing intake checks 792 fields against the pinned [BaBar v2 source](https://arxiv.org/abs/1309.5278v2), `paper.tex`, `tab:totZuc` and `tab:totZul`. These include reported P means. The means remain source metadata: an integrated amplitude cannot be replaced by a prediction at mean P. Both charge-pair yields must be integrated separately before constructing the normalized UL/UC ratios, as in [KPSY Eq.82](https://arxiv.org/html/1505.05589v1#S2.E82).

The declared conditional reconstruction uses `0 <= qT <= 3.5 GeV`, mapped to P with the fixed reported z1 mean and retaining `P*dP=z1^2*qT*dqT` in absolute yields. The upper cut is printed in the primary selection. Zero is an explicitly declared radial origin, not a recovered lower experimental acceptance. Uniform momentum acceptance, fixed z means and factorization of the primary mean-of-angular-ratio scalar remain model assumptions. The source-literal Q2=112.36 remains conditional. No source value, error component, angular definition or selection changes.

The [contract](../config/historical-sia-integral-contract.json) requires separate momentum refinement, native/matching/radial refinement and finite tails, plus independent adaptive-qT scalar assembly. It preserves the existing 1e-4 physical and 1e-7 momentum/derivative targets. All coarse misses remain visible. Representative all13 gradient/HVP checks are declared at three source geometries with both ratios; their scope is explicit. The one-sided fourth-order reference avoids crossing the rounded Table1 boundaries and must first pass independent mathematical controls.

The 36 exact geometries contain 14 distinct z queries. A future fixed-query preparation cache requires exact source/settings identity, unchanged separate-query comparators and mutation checks. This inventory alone demonstrates no runtime speedup or production throughput. Native state remains isolated, with all native/provider work outside parameter evaluation.

Failed evidence must survive scheduler cleanup. Cache controls must reject rehashed window/Jacobian/source changes and point substitution. A source review, synthetic closure or a completed P-point diagnostic cannot approve this integral path.

The original BaBar P-point refinement misses and paper-center/coupling failures remain separate. Missing author code, covariance/overlap, acceptance, corrected Belle primary bytes, production gates and the complete B0-B4/V1/uncertainty/publication scope remain open.
