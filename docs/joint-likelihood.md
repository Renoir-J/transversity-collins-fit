# Joint predictions and profiled likelihood

`KPSY-JOINT-PROFILE-001` joins prepared SIDIS, SIA, STAR200 and STAR510 observations under the same thirteen direct physical parameters. It is a prerequisite to the real-data fitting executor. The current likelihood accepts **explicitly synthetic observations only**; historical/experimental candidate and production scopes fail. There is no implemented real-data fit, new extraction, uncertainty coverage or production admission.

## Shared predictions

`JointObservation` assigns a unique point identifier and source to an existing point/bin plan. `JointPredictionModel` preserves point order and rejects inconsistent physical model settings. Numerical quadrature orders may differ between children. Child snapshots and the fixed observation/covariance/nuisance snapshots are rechecked before use; changed arrays or provenance fail.

Leptonic observable/bin phases and STAR bin phases retain their existing definitions. A STAR point requires an explicit `sin_phiS_minus_phiH` or reverse phase. No second correction is applied to a prepared observable. The simultaneous global sign flip of transversity and Collins normalizations remains ambiguous for these processes. Unsupported observations fail rather than being omitted. In particular, a full-x lattice tensor charge requires a separately validated endpoint/extrapolation plan; measured-x partial moments cannot replace it.

## Statistical contract

For fixed observations y, theory t(theta), data covariance C, additive design A, relative design R and nuisance prior covariance S:

```text
r = y - t(theta)
nu ~ N(0, S), epsilon ~ N(0, C), independently
y = t(theta) + [A + diag(t(theta))*R]*nu + epsilon
S = L*L', B(theta) = [A + diag(t(theta))*R]*L
chi2(theta) = min_a (r-B*a)' C^-1 (r-B*a) + a'*a
nu_profile = L*a_profile
```

C and S must be explicitly supplied, finite, symmetric and positive definite. The constructor does not infer missing correlations, symmetrize experimental matrices, clip eigenvalues or resolve overlapping views. Multiplicative errors scale with **predicted theory**. Julia differentiates that dependence through the complete profile. This is a Gaussian nuisance-profile chi2, without a marginal log determinant; it is not a marginalized Bayesian posterior.

`joint_subset` retains the selected marginal covariance, every nuisance column and the full nuisance prior. `joint_holdout` returns the training subset after removing one explicit process. Neither assumes independent process blocks or conditions covariance on held-out measurements. No experimental normalization/polarization uncertainty has been assigned a nuisance model by this component.

## Whole-job validation

Freeze [the statistical contract](../config/joint-likelihood-contract.json) and [the input preparation contract](../config/joint-pinned-input-contract.json) before running candidates. Julia owns prediction, profiling, derivatives and independent references. Python/Snakemake schedules complete jobs.

```bash
# Pinned native/input builds; existing valid build receipts may also be supplied.
python scripts/prepare_parton_inputs.py
python scripts/prepare_hoppet.py
snakemake --snakefile workflow/joint-likelihood.smk --cores 2 \
  --config run_dir=runs/joint-example
```

The analytic job uses explicitly synthetic source rows. The input job uses pinned CT10/DSSV08/DSS14 and native matching/evolution rows, with explicit reference extrapolation. **Both retain synthetic observations, covariance, kinematics, angular factors and bin weights.** The small fixed preparation orders test statistical integration, not physical convergence or throughput.

Each producer records 22 checks: independent augmented QR profiling, correlated nuisance/covariance identities, zero-nuisance reduction, nuisance-coordinate invariance, global-sign ambiguity, ordered subsets and process holdouts, source/mutation rejection, complete thirteen-coordinate gradient/HVP and source/build stability. Derivative references use 256-bit QR and two frozen finite-difference steps with Richardson cancellation. Raw function values and both coarse derivatives are retained. Thresholds remain 1e-11 for algebra and 1e-7 for derivatives, with a 1e-8 relative floor.

The separate report auditor recomputes profile and retained finite-difference algebra, checks exact comparison names/policies, and verifies current source/build dependencies. Altered reports must reject. This report audit does not replay process predictions or the finite-difference function calls; those remain the complete producer's responsibility. Fresh paths are required, failed evidence is retained, and `accepted.json` is created only after producer, audit and altered-report controls succeed. It certifies this synthetic validation scope only.

The default workflow target revalidates both cached report/audit pairs on every invocation, including source and native dependencies. Use `kind=analytic` or `kind=inputs` for a single path. Direct internal job targets do not perform this final cache check. `python scripts/check_cached_joint_workflow.py runs/joint-example` changes each report in turn, requires cached rejection and restores the original bytes before requiring success. Use a fresh run directory when code or scientific inputs change.

The report reader validates stored H*v using the separately frozen [IEEE roundoff model](../config/joint-report-roundoff-contract.json). Exact rational arithmetic on the binary64 operands checks the bound gamma(2n)*sum(abs(Hij*vj)) plus a conservative subnormal allowance, with gamma(k)=k*u/(1-k*u), u=2^-53. This permits different legal BLAS summation orders without requiring identical bits. Nonfinite operands, possible overflow and violations reject. The earlier cross-machine direct-auditor failure is retained; producer comparisons and scientific thresholds are unchanged. Cancellation and corruption tests check the rule independently of a particular report.

## Remaining fitting work

The existing finite synthetic-kernel fitter is not connected to these physical prediction plans. A real-data loader, a bounded optimizer that preserves closed physical boundaries, replica fitting and interval-coverage checks remain to be implemented and validated. The original KPSY fitter/modified HOPPET and unrounded optimum remain unlocated. The two original paper-center failures and six coupling-study center failures remain FAIL under their unchanged thresholds. Acceptance/binning, covariance/overlap, low-scale policies, full-x lattice moments and production performance remain unresolved. Production, physical extraction and G2PERF stay false.
