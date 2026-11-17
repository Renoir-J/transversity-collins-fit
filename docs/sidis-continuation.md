# Immutable SIDIS continuation

The interrupted complete132 attempt contains thirteen finished native123 jobs
and an incomplete fourteenth job. Its original141 files remain immutable.
This operational addendum introduces a fresh run namespace and an explicitly
bound complete-job adapter. It changes no numerical source, scientific policy,
parameters, measurements or acceptance threshold.

The prospective policy is `config/sidis-continuation-contract.json`. It must be
committed before implementing or executing the adapter. The supported input is
the independently reviewed attempt2 prefix at frozen numerical head `bec12db`;
attempt1 records cannot replace it. The earlier user-authorized thirteen-job
review establishes saved-record integrity only. It is not a full132 gate,
permission to overwrite original files, fitting authorization or a PR merge.

A typed Julia path argument resolves the thirteen original completed native123
folders to their actual original locations. Every other folder resolves to the
new namespace. Additional methods apply only to this new argument type; no
existing method is replaced, no global is monkey-patched, and no filesystem
link or modified receipt is created. The adapter verifies old jobs at their
actual original paths and rejects dispatch to any original folder. Point14 is
recomputed into a fresh real folder. Every new
job and cached review verifies that the exact original141-file inventory and
hashes are unchanged.

Julia owns all numerical work. New point producers and auditors use the
unchanged frozen package, request and producer/auditor functions. The separate, versioned
global adapter resolves each old/new point at its actual location and uses the
original producer decisions and independent scalar/decimal readers. The adapter
code has its own committed bindings; reporting the numerical fingerprint never
conceals a changed transport or reader. No frozen module is monkey-patched.

The complete original phase policy remains mandatory: all132 first native
pairs, global selection and all132 fourth native pairs if required, common-level
finite tails, all132/all13 candidate and independent256-bit center/inward
references, fourth references if required, source-math controls, and the complete
independent ledger. Only then may comparison CSV and three experiment figures
be generated. Cached use requires fresh saved-arithmetic readers. The original
native3 target miss is retained; it cannot be excused by successful storage or
synthetic checks.

All original time/RSS/availability limits remain in force. New failed or
interrupted folders remain occupied and reject reuse. A monitor timeout alone
never permits restart. No code in this workflow starts fitting or changes PR
readiness. Numerical continuation execution is separate from software checks
and from the previously authorized read-only review.

Missing author code/unrounded states, paper-center/coupling failures, source
transcription, conditional target/acceptance, covariance and projection overlap
remain limitations. The full historical262 fit and B1–B4/V1, uncertainty and
publication deliverables remain unfinished.

## Execution

Use Linux and the existing locked Python/Julia environments. Keep two clean
checkouts: this adapter branch and the exact frozen numerical head. All335
original committed files must match in both checkouts; the adapter only adds
files. The explicit frozen-checkout argument selects trusted source code.
Cached manifests cannot select an executable or an alternative checkout.

Each fresh raw-source replay writes its request and receipt under the new
namespace's `source-replays/` directory. The original request remains an input.
The adapter does not call the legacy Python request verifier, which creates a
temporary directory beside its input request. Original job metadata receives
the same source/input/output/resource checks without those temporary writes;
the unchanged complete Julia reader still audits every original packet.

A per-job `job.completed` marker is created only after both complete commands
and all metadata checks succeed. The original Julia global reader requires it.
Cached verification requires the marker plus the full source, resource, output
and arithmetic evidence. A forged marker never repairs a failed invocation.

Prepare a fresh directory under the frozen checkout's ignored `runs/` tree:

```sh
python scripts/run_sidis_continuation.py prepare \
  --frozen-checkout "$FROZEN" \
  --original "$FROZEN/runs/historical-sidis-jobs-real-data/attempt-2" \
  --review "$REVIEW" --integrity "$INTEGRITY" --loss "$LOSS" \
  --native-build "$NATIVE" --parton-build "$PARTON" \
  --output "$FROZEN/runs/sidis-immutable-continuation/attempt-1"
```

The review, integrity and loss files must have the exact prospectively pinned
digests. Preparation checks original process absence, source bytes, native/input
receipt identities, complete source replay and each original job's metadata.
It dispatches no numerical job and grants no scientific acceptance.

Inspect the plan, then run the fresh typed-path Julia prefix reader:

```sh
python scripts/run_sidis_continuation.py plan \
  --frozen-checkout "$FROZEN" --manifest "$MANIFEST"
python scripts/run_sidis_continuation.py prefix \
  --frozen-checkout "$FROZEN" --manifest "$MANIFEST"
```

Each numerical job writes only into its canonical new folder. For example,
the first remaining source point can be run and freshly reviewed with:

```sh
python scripts/run_sidis_continuation.py run --phase native123 --point-index 14 \
  --frozen-checkout "$FROZEN" --manifest "$MANIFEST"
python scripts/run_sidis_continuation.py verify --phase native123 --point-index 14 \
  --frozen-checkout "$FROZEN" --manifest "$MANIFEST"
```

Run the complete serial phase policy through Snakemake:

```sh
snakemake --snakefile workflow/sidis-continuation.smk --cores 1 \
  --config frozen_checkout="$FROZEN" manifest="$MANIFEST"
```

The standalone `all` action schedules the same complete jobs. Completion markers
are bookkeeping only; cached workflow use invokes the fresh complete Julia
ledger audit again. The existing shared numerical-checkout lock is opened
read-only, and additional lock/log/report files live in the new directory.
An incomplete job remains occupied and explicitly stops subsequent use; this
version does not automatically migrate a later interrupted continuation into a
third run. Every such recovery needs its own provenance and retention review.
