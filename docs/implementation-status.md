# Implementation and research status

ADR-001, 2026-11-11. Issues and milestones track progress. Software capability and scientific admission are separate.

| Component | State | Evidence / remaining work |
|---|---|---|
| Julia prediction/likelihood/AD/optimization | Implemented for finite synthetic kernels | Golden fixture, nuisance identity, derivative and HVP checks |
| Multistart, parametric replicas, process holdouts | Implemented for closure | Numerical reports / CI; replicas are not posterior samples |
| Python data validation and CSV/XLSX normalization | Implemented | Explicit mappings, units/signs, covariance and overlap |
| Snakemake 12-stage DAG / receipts / cache / outputs | Implemented | Linux integration and negative admission tests |
| Fixed Float64 batched Fortran ABI | Implemented | Cross-language application test; no HOPPET operator builder |
| GitHub progress and Git rules | Configured | Milestones, issues, PR/CI, hooks; server enforcement recorded separately |
| Updated experimental data | Pending review and normalization | Candidate/source index exists |
| KPSY CSS historical reproduction | Blocked on source or independent reconstruction | Author's modified HOPPET and fitting program not located |
| Production theory | Pending | Process-accurate bins/acceptance, evolution, matching and independent fixtures |
| Production gradients/performance/interval coverage | Pending | Synthetic checks do not certify physical implementations |
| B0–B4 / V1 real fits and publication grids | Blocked on production gates | No new physical extraction reported |
| NLO jet matching | Separate milestone | Inclusive hard functions alone are insufficient |
| Bayesian sampling / reweighting | Not implemented | Optional after closure and prior/calibration checks |

This repository and its CI artifacts describe implementation progress. Detailed local working notes are excluded from Git. Production does not silently substitute a synthetic or mismatched approximation.

