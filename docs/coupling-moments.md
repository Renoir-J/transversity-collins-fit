# Coupling-only partial-moment sensitivity

This reconstructed-theory diagnostic asks how three predefined fixed-nf
couplings change the two measured-x transversity moments at the literal KPSY
Table 1 parameters. The input parameters, CT10/DSSV08 sources, LO splitting
kernel, measured-x window and Q2=10 GeV2 are shared. It does not fit data or
select a coupling. The original code and unrounded optimizer state remain
unlocated, and historical flavor thresholds and backend flags remain unknown.
The study's baseline name `asymptotic_two_loop_existing` maps only to the
coupling receipt's `asymptotic_Lambda_0p225`; both identities are retained.

The prospective [contract](../config/coupling-moment-contract.json) was committed
before candidate physical calculations. It preserves both the failed
[paper benchmark](kpsy-paper-moments.md) and the
[native coupling diagnostic](native-coupling.md), including all thresholds.
Statistical fit intervals and printing sensitivity cannot rescue failed gates.

## Fixed-time interface

`PreparedEvolutionTimes` owns immutable Q2/tau tuples and records exact Q02,
fixed nf, alpha0 metadata, model/source identifiers and a SHA256 provenance
link. The supplied quantity is integral alpha_s/(2*pi) d ln Q2. Queries use
exact labels, permit repeats and reordering, and reject interpolation,
extrapolation, mismatched scales/nf and setting derivatives. This descriptor
does not certify its source or reconstruct alpha_s.

The opt-in native path supplies these times to the unchanged retained LO
bridge and auxiliary-time mapping. The dense path uses the unchanged
exp(tau*P) generator. Both bind profile identity in snapshots; prepared dense
operators also verify all matrix and grid bytes before use. Existing scalar
and asymptotic constructors keep their original formulas and signatures.
Point-row and Wilson-matching adapters reject this type until a separate
physical-alpha and provenance contract supports it. The historical initial
operator constructor still requires the original asymptotic prescription.

## Reproducible jobs

On Linux, install the locked Python and Julia dependencies, then run:

```sh
python scripts/prepare_native_coupling.py
python scripts/prepare_parton_inputs.py
snakemake --snakefile workflow/coupling-moments.smk --cores 1 \
  --printshellcmds --config run_dir=runs/coupling-moments/attempt-1
```

Use a fresh run directory for each changed candidate. Accepted stage outputs
are protected, and every job rejects an existing report before redirecting
its log. Scientific validation JSON is retained as a Snakemake log even on
failure. An `accepted.json` copy exists only after that diagnostic job exits
successfully; it is not a production admission artifact.

The scheduler runs three complete Julia jobs: fresh native coupling checks,
analytic/native fixed-time controls, then physical moment sensitivity. It
never calls Python or a language model inside the numerical loop. A consumer
checks current sources, build hashes, complete named checks and operands,
fixed tolerances, native protocol evidence and independent beta references.
Failed, incomplete, absent or stale dependencies reject physical execution.

The physical ladder shares all initial arrays across the three models at each
resolution. It retains coarse trials, independent dense evolution and log-x
quadrature, native adjoint operands and all 13 parameter derivatives. One-sided
256-bit derivative references respect the Nh_d=-1 boundary; the seven
Collins-only derivatives vanish for these transversity moments. Coupling,
grid and evolution-time derivatives remain unsupported.

Overall status is FAIL, with exit 2, if any numerical control or any of the six
model/paper-center comparisons fails. A passed numerical discretization does
not reproduce a paper center. A variant closer to a center cannot replace the
failed baseline or establish the historical prescription.

## Scope of results

The local candidate uses 1,024 native intervals and the 3,201-node dense
background after all three frozen initial refinement levels. At Q2=10 GeV2:

| Fixed coupling | Native up partial moment | Native down partial moment |
| --- | ---: | ---: |
| Existing asymptotic | 0.2885459107 | -0.2128577231 |
| Native RG, literal alpha0=0.327 | 0.2885626656 | -0.2128093975 |
| Native RG, matched asymptotic alpha0 | 0.2885625485 | -0.2128097366 |

The maximum coupling-only changes are 1.676e-5 and 4.833e-5, respectively.
All three models still fail both printed-center comparisons against
+0.30/-0.20 under the unchanged absolute 0.005 criterion. This bounded
comparison does not establish a global cause or recover historical settings.
The baseline matches the failed original benchmark exactly.

The complete physical report remains **FAIL**: 76 numerical checks pass and
six paper-center checks fail. The Julia producer exits 2 and Snakemake exits
1; no sensitivity `accepted.json` is created. All failed JSON, job logs,
coarse trials and nonterminal progress remain available as ignored artifacts.
Independent Julia review recomputes all 82 checks, source/input bindings,
log-x quadrature and 256-bit one-sided derivative operands, and rejects
12 false-admission/report counterfactuals. It confirms this scientific failure.

The preceding current-source coupling comparator passes 36 checks; the
prepared-time controls pass 49. The first development run reused an existing
directory and overwrote its summary when native batch replacement was
rejected; raw batches and invocation hashes were retained, and that run is
not qualifying evidence. The next run passed the numerical controls but its
baseline namespace mismatch was rejected by the consumer before physical
execution. The corrected mapping retains both identities and a fresh full
run passes the dependency checks. Neither fix changed the frozen policy.

Linux regression checks pass 158 Python tests, workflow policy, 1,892 Julia
assertions in 44 testsets including 66 new prepared-time assertions, and the
existing 47-check retained-native regression. The full synthetic closure and
cache rejection/restoration checks pass: all 12 receipts, 33 outputs and
15 input references are verified against the current content fingerprint.
These were run on the hashed working tree before the implementation commit;
the review separately binds those file bytes to the final Git commit.
Remote CI remains pending and cannot be inferred from local checks.
Synthetic results cannot change the six scientific failures. Source-bound
reports and raw inputs are not committed fixtures.

The project still needs complete historical observable benchmarks,
dataset-specific acceptance and covariance/overlap decisions, a real-data
fitting executor, uncertainty coverage and production/performance admission.
These partial moments are neither a new extraction nor a full-x tensor charge.
