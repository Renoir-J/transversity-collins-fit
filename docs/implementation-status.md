# Implementation and research status

ADR-001, 2026-11-11. Issues and milestones track progress. Software capability and scientific admission are separate.

| Component | State | Evidence / remaining work |
|---|---|---|
| Julia prediction/likelihood/AD/optimization | Implemented for finite synthetic kernels | Golden fixture, nuisance identity, derivative and HVP checks |
| Multistart, parametric replicas, process holdouts | Implemented for closure | Numerical reports / CI; replicas are not posterior samples |
| Python data validation and CSV/XLSX/HTML/YAML normalization | Implemented for first experimental intake | Twenty-two pinned sources; 773 amplitudes in 82 overlapping views; source corrections and primary/JAM cross-checks |
| Snakemake 12-stage DAG / receipts / cache / outputs | Implemented | Linux integration and negative admission tests |
| Fixed Float64 batched Fortran ABI | Implemented | Cross-language fixed-operator application; separate retained-HOPPET construction contract |
| LO collinear evolution | Reconstructed component; not production admitted | Analytic cell integration, x-space/Mellin refinement and initial-parameter derivatives; [scope and limits](evolution.md) |
| Two-loop coupling / CSS factors | Reconstructed component; not production admitted | Explicit asymptotic coupling, signed Sudakov integral and channel-specific nonperturbative factors; [scope and limits](css.md) |
| CT10 / DSSV08 / DSS14 inputs | Reconstructed interfaces; not production admitted | Locked source builds, native Julia legacy interpolation, official LHAPDF C ABI, weight/flavor checks and explicit range policies; [scope and limits](parton-inputs.md) |
| SIDIS/SIA matching and radial structure functions | Reconstructed components; not production admitted | Process-specific Wilson coefficients, J0/J1/J2 transforms, charge-conjugate flavor sums and CSS assembly; [scope and limits](leptonic.md) |
| Initial h1/Hhat3 parameter family and fixed-scale evolution | Reconstructed components; not production admitted | Direct physical coordinates, exact closed boundaries, strict pinned inputs and generic matrix evolution; [scope and limits](initial-conditions.md) |
| Retained HOPPET kernel/linear ABI and query-row caches | Reconstructed components; not production admitted | Explicit weight/time mapping, checked native transpose and detached Julia rows; global-grid physical failures retained; [scope and limits](hoppet.md) |
| Evolved polarized Wilson matching rows | Reconstructed components; not production admitted | Explicit process/physical-alpha metadata, retained regular convolution and independent Mellin/dense/input/AD checks; [scope and limits](matched-rows.md) |
| Batched retained preparation and polarized leptonic contexts | Reconstructed implementation; not production admitted | Actual radial scales, signed complete structure functions, charge/target assembly, generic parameter derivatives and separate analytic/input validation producers; [scope and limits](polarized-processes.md) |
| Normalized leptonic observables and finite P integration | Reconstructed implementation under validation; not production admitted | Explicit SIDIS depolarization choice, hadron2 SIA U/L/C ratios, projected P*dP quadrature and independent analytic/pinned-input producers; [scope and limits](leptonic-observables.md) |
| GitHub progress and Git rules | Configured | Milestones, issues, PR/CI, hooks; server enforcement recorded separately |
| Updated experimental data | First candidate batch acquired; production review pending | HERMES/Belle/BESIII and official STAR200/510 pion Collins tables; STAR physical binning, COMPASS tables, covariance and admission remain open |
| KPSY CSS historical reproduction | Blocked on source or independent reconstruction | Author's modified HOPPET and fitting program not located |
| Historical 262-amplitude data intake | Candidate implementation; no fit admission | Source-locked historical views, published cut, Hall A primary remap and explicit uncertainty/convention blockers; [scope and limits](historical-data.md) |
| Production theory | Pending | Process-accurate bins/acceptance, evolution, matching and independent fixtures |
| Production gradients/performance/interval coverage | Pending | Synthetic checks do not certify physical implementations |
| B0–B4 / V1 real fits and publication grids | Blocked on production gates | No new physical extraction reported |
| NLO jet matching | Separate milestone | Inclusive hard functions alone are insufficient |
| Bayesian sampling / reweighting | Not implemented | Optional after closure and prior/calibration checks |

This repository and its CI artifacts describe implementation progress. Detailed local working notes are excluded from Git. Production does not silently substitute a synthetic or mismatched approximation.

