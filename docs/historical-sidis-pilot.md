# Historical SIDIS resource pilot

This prospective operational policy measures complete Julia jobs before scheduling the full 132-point diagnostic. It inherits the unchanged scientific contract at `c2d2281` and source/model interface at `452fef3`. The first retained HERMES, COMPASS proton and COMPASS deuteron rows are selected by source order before any pilot predictions. Their distinct kinematics, errors, source labels and conditional phase/target mappings remain unchanged.

Each case prepares one geometry at the already declared finest native/matching/radial setting, at finite bound 24. It compares signed UU/FUT/source-amplitude values with separate generic scalar assembly sharing the same native rows, then records cold and repeated value/gradient/HVP evaluation times. Every derivative evaluation remains Julia arithmetic on prepared rows. These AD samples have no independent physical derivative certificate. One native level and three rows cannot satisfy the full 132-point ladder, tails, all132/all13 derivative references or coverage requirements.

Python/Snakemake schedules whole Julia producer and audit jobs. GNU time records each complete command's direct-child maximum RSS and wall/user/system seconds, separately from coefficient payload bytes and Julia phase timings. The policy fixes a shared one-hour budget per case, a 3 GiB observed command-RSS limit, at least 4 GiB available before dispatch, and one pilot per checkout. These engineering limits are declared before execution; exhaustion is FAIL and retains evidence. They are not physical acceptance thresholds or a speedup certificate.

Successful measurements use `PILOT_COMPLETED_NOT_DIAGNOSTIC_GATE`, with an independent saved-arithmetic audit and a complete three-case ledger. All production, fitting, full132, derivative-accuracy, performance and coverage admissions remain false. Failure, timeout, parent signal, missing resources and forged completion markers retain or reject evidence explicitly. The policy was frozen at `278cd90` before this implementation. Complete validation and actual resource measurements remain separate from policy publication.

## Complete jobs

From the repository root on Linux, with GNU time and the existing pinned build dependencies available:

```sh
snakemake --snakefile workflow/historical-sidis-pilot.smk --cores 2 \
  --config run_dir=runs/sidis-pilot-example
```

The workflow acquires exact primary bytes, replays historical intake, prepares a full132 source request and runs three complete Julia producer/audit jobs. Cases are chained in their frozen source order even when more cores are available. A nonblocking checkout lock and a memory observation guard dispatch. Existing verified native/parton receipts and raw caches can be supplied with `native_build`, `parton_build`, `cache`, `primary_cache`, `seed_root` and `offline=True`. No earlier failed directory is overwritten; choose a fresh run directory for a new attempt. A corrupt cache fails explicitly, and missing offline inputs never receive a fallback.

The final target always replays source validation and the saved-arithmetic audit, including cached use. GNU time wraps the caller-selected current Julia executables; cached metadata never selects a new executable. Retained resource records include child exit, wrapper outcome, interruption and RSS. A cooperative child exit0 after timeout or parent signal remains a failed command. Scalar reference comparisons retain the inherited precision; AD samples are checked for shape, finiteness and repeatability only. The reader does not reexecute evolution, independent scalar assembly, physical derivative references or resource measurements.

The three-case ledger preserves each source identity, invocation/report/audit hashes, complete-command RSS/time and Julia phase observations. Any all132 storage projection is explicitly conditional on these three samples and is not measured all132 memory. This pilot supplies resource evidence for the complete diagnostic's job layout and budgets. It supplies no likelihood, fit, uncertainty, tensor charge, performance gate or full132 acceptance.
