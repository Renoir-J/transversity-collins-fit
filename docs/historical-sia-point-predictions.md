# BaBar fixed-mean SIA point diagnostics

`KPSY-HISTORICAL-SIA-POINT-001` was frozen in `2e93e1d` before new BaBar numerical predictions. The implementation connects 18 RF0 transverse-momentum amplitudes, nine UL and nine UC, to the existing Julia SIA prediction plans at the rounded KPSY parameters. It is reconstructed theory with real historical source rows, not a likelihood, fitted optimum, experimental admission or production extraction. Numerical validation of this new path remains pending until a complete source-bound job and its audit finish.

The existing historical intake checks 162 primary P-table fields for these rows, including both z means, P means/edges, the angular mean and signed central/statistical/systematic values. The [BaBar v2 primary source](https://arxiv.org/abs/1309.5278v2), `paper.tex`, RF0 section of `tab:totPt`, is independently hashed. The source-literal legacy Q²=112.36 GeV² remains unverified by that selected table. It is neither replaced by 10.58² nor relabeled as a primary-verified beam energy.

The report retains the mean of the angular ratio, its publication measure and unresolved joint acceptance. For unlike U, like L and all C=U+L categories, it applies that angular factor once to `C_U/D_U-C_L/D_L` or `C_U/D_U-C_C/D_C`. Each signed Collins numerator and positive unpolarized denominator comes from all four pion charge pairs. The inherited J2, z/mass, matching and CSS factors are unchanged. This linearized, factorized comparison at reported means does not reconstruct the experimental joint acceptance. The integrated SIA source views require their own P integrals; their displayed mean P is not a point replacement.

## Complete jobs

```sh
snakemake --snakefile workflow/historical-sia-points.smk --cores 2 \
  --config run_dir=runs/babar-point-example
python scripts/check_cached_historical.py config/historical-sources.json \
  data/raw/historical runs/babar-point-example/historical
python scripts/check_cached_historical_sia_points.py \
  runs/babar-point-example/request.json runs/babar-point-example/diagnostic
```

Existing current native/parton build receipts and a verified raw cache can be supplied using the corresponding workflow settings. The default final target replays source preparation and the saved numerical audit even when the completion marker is cached. Julia owns every prediction and numerical comparison; Python schedules whole jobs and renders completed values. No provider/native, Python or language-model callback enters a parameter evaluation.

At least three frozen native/matching/radial refinements must pass the unchanged 1e-4 target for every pair yield, category sum and final amplitude. The finite 24-to32 inverse-GeV radial extension has the same target. Budget exhaustion remains FAIL and every coarse result remains visible. Independent scalar charge/CSS/Bessel composition shares prepared native rows and initial functions; it is not independent evolution. The 256-bit synthetic angular fixture projects the derivative of normalized-yield double ratios at zero perturbation and retains all seven counterfactuals. It validates the declared linearization convention, not finite-amplitude experimental extraction or acceptance.

The producer requires 185 checks. The separate Julia reader recomputes 166 retained arithmetic comparisons and distinguishes 19 producer object-reuse/temporal attestations. It rehashes current sources, requests and build dependencies but does not reexecute native evolution, predictions, primary parsing or Fourier quadrature. Sign and first-six-coordinate controls are fixed-value probes, not global-optimality or production-derivative certificates.

Only a completion marker is declared as the numerical scheduler output. Failed reports, invocations, progress and logs survive cleanup. Fresh request/output paths are required. Six cache controls reject altered geometry, angular definition, Q-status label, predictions, figures and production scope, then restore the original bytes and timestamps. A separate explicitly synthetic failed-job control checks FAIL retention and forged-marker rejection without executing numerical predictions.

## Remaining scope

Statistical and reported systematic errors remain separate, with UL/UC overlap and missing covariance preserved. No χ², pulls, normalization nuisance, fitting or uncertainty interval is constructed. The actual comparison table/figure is produced only after the numerical report and audit succeed; it shows discrete fixed-parameter values and statistical bars, with systematic errors in the table.

The historical 262-point fit, updated SIDIS/SIA, STAR200/510, full-x paired lattice fits, robustness and uncertainty products remain required. Missing author code/modified HOPPET/unrounded states, reconstructed homogeneous Collins evolution, input continuation, corrected Belle primary bytes and existing paper-center FAIL results remain open. Production, extraction, fitting, performance and interval-coverage admission stay false.
