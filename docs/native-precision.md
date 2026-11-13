# Native precision qualification

The project adapter previously called `InitDglapHolder` to build ordinary LO
DGLAP splitting functions, deleted that matrix, and replaced it with the
contracted transversity kernel. Phase-labelled diagnostic builds of the parent
adapter locate 14 precision-error messages in this unused initialization across
eight single/nested-grid and independent/batched construction cases. Those
diagnostics do not establish the origin of every warning in earlier runs.

The pinned integrator prints `TOO HIGH ACCURACY REQUIRED` and returns zero on its
failed adaptive-integration branch. ABI status zero therefore cannot by itself
qualify a native calculation. Earlier comparison results remain archived with
their warning counts; the historical runs are not relabelled as passing native
precision qualification.

`KPSY-NATIVE-PRECISION-001`, nested in `config/hoppet-contract.json`, specifies a
project-owned initialization of a fixed-nf LO holder. It uses the pinned public
derived types, `SetNfDglapHolder`, zero `InitSplitMat`, and the existing project
transversity kernel. Three unused coefficient-function headers are explicitly
initialized for safe native deletion. The upstream source, grids, interpolation,
evolution, matching kernels, alpha_s, eps/du and existing acceptance thresholds
are unchanged.

Julia captures native stderr under the existing process lock, flushes C stdio,
replays diagnostics, and rejects precision-error messages. A returned context is
registered before checking those messages, allowing a rejected constructor to
close it. The explicitly synthetic C fixture returns ABI status zero and a live
context while printing the error; tests require both rejection and context
release. The fixture provides no physics implementation or predictions.

The reference producer runs unchanged parent package code with the same driver
and frozen policy. Its status is
`REFERENCE_COMPLETED_WITH_DIAGNOSTICS_NOT_ADMITTED`; it is a comparison reference,
not an admitted native-precision gate. The candidate must retain every exported
Float64 bit of the fixed-grid forward, selected-row, generator, interpolation,
adjoint and applicable Collins matching payloads. Independent Mellin and
row/forward comparisons retain their declared tolerances, and every candidate
native diagnostic count must be zero. Fixed nf3 through nf6 cases also use the
independent Mellin reference. New kinematics still require their own checks.

The dedicated Linux workflow builds both adapters against the same unchanged
upstream pin and archives the failed or successful reports and diagnostics.
Standard CI and the existing physical-input and complete-process workflows
remain separate required evidence. A completed build, synthetic fault-injection
test, or parent comparison cannot replace them.

The STAR common-grid near-endpoint limitation, missing historical author code,
independent physical benchmarks, measured acceptance, data covariance/overlap,
and production performance/coverage remain open. Production readiness, physical
extraction and G2PERF remain false.

Pinned source evidence:

- [Holder initialization and deletion](https://github.com/hoppet-code/hoppet/blob/feda4f2f64d0131f489ad1cf7728ec6be186cc56/src/dglap_holders.f90)
- [Zero splitting-matrix initialization](https://github.com/hoppet-code/hoppet/blob/feda4f2f64d0131f489ad1cf7728ec6be186cc56/src/dglap_objects.f90)
- [Native adaptive integrators](https://github.com/hoppet-code/hoppet/blob/feda4f2f64d0131f489ad1cf7728ec6be186cc56/src/integrator.f90)
