# Bounded joint fitting

`fit_joint` optimizes the prepared SIDIS/SIA/STAR thirteen-parameter predictions and correlated Gaussian nuisance-profile likelihood. Julia owns the predictions, solver, derivatives and postchecks; uncertainty computations also belong in Julia. Snakemake launches complete Julia jobs. Python and language models do not enter a numerical callback.

**Current scope is synthetic closure with reconstructed theory.** This is a fitting executor for the physical parameter family, using explicitly synthetic source rows, geometry, observations and covariance. It has no real-data loader/admission, physical extraction, replica coverage, full-x lattice moment or production performance certificate. The original paper-center failures and missing author-code/covariance/acceptance/overlap limitations remain.

## Direct coordinates and stopping

The domain is `Nh_u,Nh_d ∈ [-1,1]`; all eight shape coordinates are nonnegative; `Nc_u,Nc_d` are unbounded; `gc < gh`. The binary64 upper limit for `gc` is `prevfloat(gh)`. There are no finite surrogate bounds, transformed fit coordinates, post-hoc clipping, fabricated boundary gradients or penalty objectives. Finite direct objective values at the zero transversity-shape boundaries remain supported even though their normalization derivatives are singular.

The original default is `LN_SBPLX`. Explicit `algorithm=:LN_BOBYQA` selects the separately frozen local quadratic-model comparison and tighter operational tolerances. Its native name and complete settings are recorded. It is not an automatic fallback or a production selection. BOBYQA's quadratic model has no guarantee at nonsmooth boundaries. Both algorithms are local methods; neither proves a global optimum.

BOBYQA internally rescales unequal steps and prepares near-bound starts. The requested start and actual physical callback coordinates remain visible. The [initial arithmetic contract](../config/bounded-initial-roundoff-contract.json) binds the native source/version and replays that first-coordinate map exactly in binary64. Fifty-one assertions cover twelve native cases; the reader adds no empirical ULP tolerance. The original retained reader failure and incomplete fitting report remain preserved.

Every supplied start is retained. `FTOL_REACHED` and `XTOL_REACHED` are configured tolerance stops; `MAXEVAL_REACHED` and `MAXTIME_REACHED` are **INCOMPLETE**. Negative/unknown stops, including roundoff limitation, and callback/nonfinite/source failures remain **FAIL**. A positive native return code alone cannot pass the job. A finite candidate remains diagnostic when the job is incomplete. Returned objectives are reevaluated and native/callback counts must agree.

The [original contract](../config/bounded-joint-fit-contract.json) was committed before any candidate. Its 4,000- and 20,000-evaluation Subplex jobs both exhausted their budgets and failed parameter recovery; those results are not promoted. The [solver-comparison contract](../config/bounded-solver-comparison-contract.json) was committed before BOBYQA candidates. It preserves every physical and scientific criterion. Comparing distinct stopping tolerances is not a matched-accuracy performance ranking.

## Whole-job execution

Use the locked Linux environment and a fresh directory. Algorithm and resource budget must be explicit. Candidate/production modes reject.

```sh
snakemake --snakefile workflow/bounded-joint-fit.smk --cores 1 \
  --config run_dir=runs/bounded-example algorithm=LN_BOBYQA maxeval=20000
python scripts/check_cached_bounded_workflow.py runs/bounded-example \
  --algorithm LN_BOBYQA --maxeval 20000
```

The job constructs 24 SIDIS, 24 SIA and four STAR observations, predicts a noiseless frozen truth, and fits from two fixed starts. The full covariance and correlated additive/theory-scaled nuisance model are synthetic. Its whitened/scaled prediction Jacobian must have local rank13; this is not global identifiability or coverage. The unchanged recovery criteria are chi2≤1e-8, maximum absolute prediction error≤1e-6 and maximum scaled parameter error≤0.005 up to the simultaneous transversity/Collins sign flip. Fitted signs are never canonicalized.

The producer retains both complete theta/chi2 evaluation traces, counters, postchecks, nuisance profile, local singular values and a separate deliberately incomplete finite boundary test. The reader checks source/contract/runtime bindings and all trace domains, replays returned candidates and the selected predictions, independently recomputes augmented least-squares QR and all eleven check outcomes. It does **not** replay every callback objective or certify global optimality. The reader's status preserves the producer's failure/incomplete status. Altered-report controls must reject for the expected reason.

`accepted.json` is copied only after producer, independent reader and controls succeed. The final Snakemake rule has no output and revalidates cached reports and the requested algorithm/budget on every invocation. The separate cache check changes a numerical result, requires rejection and restores the original bytes. Wrong algorithm/budget and production mode/scope also reject. Failed reports are kept; a failed directory is not silently overwritten or resumed as a passing one. All receipts retain synthetic scope and false production/extraction/G2PERF flags.

## Integrity and remaining work

The immutable Born subtree is fingerprinted by its owning constructor once; callers cannot supply a trusted digest. Every mutable point/operator/sample/fragmentation/provenance guard still executes. Six fixed cases preserve predictions/objectives and supported interior derivatives/HVPs; the retained before/after local comparison is bitwise identical in all 23 comparisons. Its timing is informational, not G2PERF or an experimental fitting-throughput result. The [cache contract](../config/star-immutable-born-cache-contract.json) records the unchanged numerical targets and limitations.

The existing 32-replica toy closure is separate evidence. It does not establish uncertainty coverage for this thirteen-parameter model. Real-data source/phase/bin/acceptance/covariance/overlap review, independent process accuracy, bounded replica fitting and coverage, full-x endpoint/tail validation, paired lattice fits and final publication products remain required. Production readiness must continue to fail explicitly until its actual prerequisites are satisfied.

Primary software references: [NLopt algorithms and comparison limits](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/), [native stopping and resource controls](https://nlopt.readthedocs.io/en/latest/NLopt_Reference/), [pinned Julia wrapper](https://github.com/jump-dev/NLopt.jl/tree/v1.2.1).
