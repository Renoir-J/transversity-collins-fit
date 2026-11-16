# Source-bound fixed context storage jobs

The Julia storage API and four-case workflow implement the separately frozen [prospective contract](../config/fixed-leptonic-storage-contract.json). They save parameter-independent prepared SIDIS/SIA contexts, verify complete data and row aliases, and compare fixed-center values, all thirteen gradients and directional HVPs before and after reload. A fresh Julia auditor loads each context and recomputes its observables and derivatives. This is storage identity evidence; all production, extraction, fitting and coverage flags remain false.

The selected cases are HERMES proton `KPSY15_SIDIS_3000:r2`, COMPASS proton `KPSY15_SIDIS_3025:r2`, free-isoscalar deuteron `KPSY15_SIDIS_4000:r2` and BaBar `KPSY15_SIA_1000:r2`. The BaBar case also includes the source UC row at the same exact geometry, with its own reported angular average. Every source request retains its full original 132-point or 18-point coverage. Request/source verification does not independently reparse the raw primary tables in Julia or admit a likelihood.

Run on Linux with the locked Python/Julia environments, a pinned native build, verified parton inputs and the existing checksum-bound raw caches. Set `PYTHONPATH` to this checkout's `src`; numerical predictions and differentiation remain inside complete Julia jobs.

```bash
snakemake --snakefile workflow/fixed-context-storage.smk --cores 1 \
  --config run_dir=runs/storage-attempt-1 \
  cache=/absolute/path/to/historical-raw-cache \
  primary_cache=/absolute/path/to/sidis-primary-cache \
  seed_root=/absolute/path/to/compass-source-seed offline=True \
  native_build=/absolute/path/to/native/build.json \
  parton_build=/absolute/path/to/parton/build.json
```

The workflow creates current source-bound requests and serializes four complete producer/auditor pairs. Each pair has the unchanged prospective one-hour total budget, a3GiB observed direct-child RSS limit and an8GiB available-memory requirement before each Julia command. A caller may tighten its time limit; it cannot enlarge it. Failed preconditions, nonzero exits, timeout or signal interruption, incomplete reports and missing resource records remain failures, including cooperative children returning zero. Use a fresh attempt path to preserve failures.

Only the scheduler writes `job.completed` after a successful pair. Every workflow invocation freshly checks all case inputs, outputs, raw resources, selected Julia executable and recorded commands, then runs new loaded-context auditors before accepting the storage ledger. A forged completion marker cannot replace missing or changed records. Cache auditors recompute loaded values/derivatives; they do not repeat cold native preparation or certify independent evolution. The selected executable and complete numerical sources must remain available and unchanged.

The manifest and binary payload remain in each case's generated `context` directory. Their receipt must supply the exact external manifest digest and seven expected binding digests. Loaded rows own their memory and retain public mutation guards; later disk changes cannot alter an already loaded context. Hashing, writing, loading and warm evaluations are reported separately. No full-fit runtime or production throughput is inferred.

The original interrupted full132/72 attempts and their replacement attempts use their frozen versions. This workflow does not insert a new cache into those live jobs. Original coarse misses, paper-center/coupling failures, missing author code/unrounded states, source continuation, error/covariance/overlap and acceptance limitations remain. Storage validation cannot authorize B0-B4/V1 fits, STAR predictions, physical intervals or publication claims. All broader scientific deliverables remain required.
