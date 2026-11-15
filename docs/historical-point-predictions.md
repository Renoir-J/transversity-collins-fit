# Historical fixed-parameter point predictions

This workflow connects eight source-checked Hall A extracted-neutron Collins
amplitudes to the existing Julia SIDIS prediction interface. It is a reconstructed
historical diagnostic at the literal rounded KPSY Table 1 parameters. It performs
no fit, constructs no covariance or likelihood, and cannot approve production or
experimental admission. The full 262-point historical fit and B1–B4/V1 remain open.

The prospective numerical policy is frozen in
[`historical-point-contract.json`](../config/historical-point-contract.json).
The input view identities, source hashes and order are fixed there. Historical
intake must first replay the complete current source set; the request retains
each reported mean, measured amplitude, separate statistical/experimental/model
errors, primary-row mapping and unresolved source-convention fields.

## Observable convention

Hall A [arXiv:1106.0363v4](https://arxiv.org/abs/1106.0363v4), Eq. 1, defines
the Collins coefficient of the polarization-normalized yield asymmetry with the
Trento `sin(phi_h + phi_S)` harmonic. The selected points are already extracted
neutron results; this workflow applies no second helium/effective-polarization
extraction.

Under the inherited **conditional massless, same-axis, fixed-y model**, the
comparison is `A_C = DNN(mean_y) * F_UT / F_UU`, with
`DNN(y) = 2(1-y)/(1+(1-y)^2)`. This is a declared approximation. The historical
intake's existing `depolarization_removed=false` label is preserved; it does not
establish the complete experimental correction chain. Low-Q² finite-mass/leptonic and
spin-axis effects, the experimental acceptance average and the original treatment
of model-systematic errors remain unverified. A function evaluated at reported
means is not a demonstrated bin average.

An explicitly synthetic, positive-count fixture checks this *conditional*
conversion with a 256-bit two-angle Fourier projection, both signs, a nonzero
Sivers harmonic and wrong-phase/missing-DNN/repeated-DNN/factor-two controls. The
uniform mathematical angle measure in that fixture is not Hall A acceptance.

## Jobs and products

Run from the repository root on Linux with the locked Python environment, Julia
project and compiler prerequisites described in the existing setup documents:

```sh
snakemake --snakefile workflow/historical-points.smk --cores 2 \
  --config run_dir=runs/historical-points/attempt-1
```

Use a fresh `run_dir` for a new attempt. Preserve earlier failures. Source intake
accepts the existing `seed_root`, `cache` and `offline` settings. Build preparation
accepts `hoppet_source`; complete unchanged pinned HOPPET sources belong under
`.runtime/`. Existing build receipts may be supplied as `native_build` and
`parton_build`; the workflow then treats them as read-only inputs. A missing
supplied receipt fails explicitly. Both the retained Fortran bridge and LHAPDF
bridge must match the current project sources, in addition to validating every
recorded build source and binary hash.

The `Historical point diagnostics` GitHub Actions workflow runs the same complete
pipeline and source/cache controls on relevant pull requests. It archives scoped
evidence on both success and failure; a failed source download or numerical check
does not select a synthetic replacement.

Python schedules complete jobs and renders completed values. Julia owns the
prepared prediction, parameter evaluation and comparison arithmetic. The job
builds no `JointGaussianLikelihood`, runs no optimizer and admits no data to one.
No Python, provider, native-evolution or language-model callback runs within a
prepared parameter evaluation. Exact full geometries share prepared native rows
between the pion charges; no rounding-based grouping or parameter-response cache
is used.

| Product | Meaning |
|---|---|
| `historical/candidates.json`, `historical/report.json` | Current historical source replay; candidates only |
| `request.json` | Bound source rows, fixed parameters, conditional convention and hashes |
| `diagnostic/predictions.json` | Complete Julia operands, refinement/tail history, controls and scoped outcome |
| `diagnostic/audit.json` | Current hashes plus recomputed stored comparison arithmetic |
| `diagnostic/invocation.json` and logs | Whole-job commands, exit codes, source/commit identity and product hashes |
| `diagnostic/point-comparison.csv` | Measured values, separate errors, fixed theory and unweighted differences |
| `diagnostic/point-comparison.png` | Discrete fixed-mean theory/data comparison with statistical bars |

The frozen joint native/matching/radial ladder requires at least three levels;
all eight points' `F_UU`, signed Collins numerator and final amplitude must meet
the unchanged `1e-4` adjacent-level target. The same target applies to the
24-to-32 inverse-GeV radial extension. Exhausting the ladder fails. Every coarse
value and missed target remains in the report; the final values use the 32
inverse-GeV bound.

The independent scalar charge/CSS/Bessel composition shares native rows and
initial-function definitions. It is not an independent evolution solution.
The audit reruns source/build hash checks and stored arithmetic; it does not
reexecute the predictions, evolution, two-angle quadrature or raw-table parsing
in Julia. Object reuse and stability during the original job are identified as
producer attestations. Python separately replays the raw source preparation,
including on cached use. These are explicit trust boundaries.

Cached products are revalidated by the `all` rule. Additional mutation controls:

```sh
python scripts/check_cached_historical_points.py \
  runs/historical-points/attempt-1/request.json \
  runs/historical-points/attempt-1/diagnostic
```

They restore original bytes and timestamps and retain rejected numerical-audit reports. Changing
a figure, rehashing altered geometry or a source-status field, rehashing a changed
prediction, or promoting the scope must be rejected.
The CI workflow then repeats the terminal Snakemake target to exercise cached
source/operand revalidation without rerunning the prediction job.

## Interpretation and remaining work

A `PASS` certifies only the declared fixed-parameter diagnostic checks at its
source fingerprint. It is not a new extraction, historical fit, uncertainty
interval, global optimum, gradient/performance certificate or coverage result.
All corresponding flags remain false. Unweighted observed-minus-predicted
differences are not pulls or a statistical test. The plot shows no connecting
curve because the reported means vary together; experimental/model errors remain
separate in the table.

Existing printed-center partial-moment failures remain failures. Original author
fitting code, modified HOPPET, unrounded states, covariance and full acceptance
information remain unresolved. The existing real-data likelihood gate remains
closed; no synthetic receipt can replace it. A real-data fitting contract and
matching numerical/source checks are still required before any optimization.
