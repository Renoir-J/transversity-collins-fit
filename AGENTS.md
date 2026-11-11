# Project rules

Follow ADR-001: Julia owns numerical predictions, derivatives, optimization and uncertainty. Python/Snakemake schedules complete jobs. No language-model calls or Python callbacks inside the numerical loop.

Keep `external/` pinned and unchanged. Do not deserialize legacy pickle files. Reference code is evidence, not a substitute for a matching physics contract.

Distinguish synthetic closure, reconstructed theory and validated production results in every report. Never turn a failed, absent or stale gate into PASS. A synthetic receipt cannot authorize a production job. Preserve missing-code and covariance limitations.

Changes to signs, mass/z factors, evolution, likelihood, cuts or dataset overlap need matching numerical checks and an updated scientific contract. Do not tune acceptance thresholds to make a failed run pass.

Use short branches (`feat/`, `fix/`, `data/`, `theory/`, `docs/`) and focused commits. Keep raw downloads, credentials, private chat exports, environments and run outputs out of git. Commit lock files and small, explicitly synthetic fixtures.

Validate with `python -m pytest`, `python scripts/plan_workflow.py validate`, Julia tests and the Snakemake closure pipeline on Linux. Production prerequisites may fail explicitly; synthetic execution must never silently replace them.

