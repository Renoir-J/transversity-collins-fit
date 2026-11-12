# Experimental data intake

The intake is an executable part of issue #3. It acquires 22 checksum-locked sources and produces **experimental candidates**, with separate statistical and systematic components. It does not construct an admitted production likelihood or extract physical distributions.

```bash
snakemake --snakefile workflow/experimental-intake.smk --cores 2
```

For an existing pinned research checkout, the acquisition rule can seed the ignored cache without modifying `external/`:

```bash
snakemake --snakefile workflow/experimental-intake.smk --cores 2 --config seed_root=.
```

The equivalent sequential command is `collinsflow intake`. Use `--offline` to require existing locked bytes and `--seed-root .` to copy verified source snapshots. Downloads go to ignored `data/raw/intake/`; candidate arrays, acquisition receipts and the intake report go to ignored `data/normalized/intake/`. Source changes fail the checksum check and need a reviewed manifest update. A missing source is not replaced by a synthetic table.

## Current coverage

| Source | Candidate coverage | Selection and limitations |
|---|---|---|
| HERMES 2020 | 200 amplitudes across 11 charged/neutral-pion views; 128 charged-pion 3D bins | Prefer the two 3D views; overlapping 1D projections are alternatives. Full systematic correlations remain unresolved. |
| Belle 2019 | 133 amplitudes across charged-pion and neutral-pion views in Tables VII–IX | Charged-pion 2D z UL is a primary candidate. UC, transverse-momentum projections and other channels are separate alternatives sharing the event sample. |
| BESIII 3.65 GeV, final revision | 22 amplitudes across z/pT and UL/UC views | Six exchanged-z UL points are a primary candidate. Preserve the union of exchanged bin orientations. |
| STAR 200 GeV, official HEPData v1 | 328 pion Collins amplitudes across 44 tables from Figures 14–18 | Component-year and combined-year views overlap. Other harmonics, non-pion channels and older-energy comparisons are excluded. |
| STAR 510 GeV, official HEPData v1 | 76 pion Collins amplitudes across all 12 tables from Figures 2–4 | Preserve charge, hemisphere and z-slice qualifiers. Horizontal plot errors do not supply physical bin edges. |
| STAR 510 GeV supplement | 14 positive-pion, positive-xF jet-pT points | Partial coverage only: the supplement supplies means and errors, not jet-pT bin edges. |
| Lattice isovector tensor charge | The existing original-paper constraint is cross-checked against the pinned JAM table | Use one isovector constraint; do not also count the correlated flavor charges. |

The 773 amplitudes across 82 views include overlapping projections and sample alternatives; they must not be summed as independent fit points. The official STAR archives are locked to [HEPData 130778 v1](https://doi.org/10.17182/hepdata.130778.v1) and [HEPData 160728 v1](https://doi.org/10.17182/hepdata.160728.v1). Current COMPASS tables remain acquisition work. The production contract remains blocked.

## Reviewed source findings

**HERMES bin-boundary correction.** In pinned JAM commit `3717fe101f6d7d4fb403c7256e3ca10de57b3e5f`, `database/sidis/expdata/3709.xlsx` has `zmin=0.049` in Excel rows 14–17. The corresponding upper-z interval is `(0.49,0.70]` in [HERMES Table V](https://arxiv.org/html/2007.07755v1#S3.T5). The manifest records four explicit old/new cell corrections, bound to the original source hash. Raw bytes stay unchanged. Checks compare all 64 corrected bins with the full 4×4×4 grid in the primary paper and require their means to lie within the bins. The similarly named `depole` column in table 3710 is mapped explicitly, without an inferred convention change.

**HERMES prediction prescription.** The [paper, Sections 3.3–3.4](https://arxiv.org/html/2007.07755v1#S3.SS3) defines structure-function asymmetries with epsilon-dependent prefactors removed. Its acceptance systematic is assessed for interpretation at reported mean kinematics. This must be respected or explicitly rederived before substituting a fully bin-integrated prediction. The 7.3% target-polarization scale error is retained separately as a shared multiplicative uncertainty. An additive fixed-B nuisance model alone does not certify its treatment. Statistical independence of 64 bins for one moment/species does not establish a complete systematic or cross-species covariance.

**BESIII version and units.** [The final v2 paper](https://arxiv.org/html/1507.06824v2#S0.T1) reports fraction-valued asymmetries and revised uncertainties. The older v1 table uses percent and has different errors. Intake uses v2 only and cross-checks all 22 central values, statistical errors, systematic errors and reported mean kinematics against four separately pinned JAM transcriptions. The first final UL z-bin result is `0.0128 ± 0.0085 (stat) ± 0.0114 (sys)`. Its angular factor is a ratio of angular averages, distinct from Belle's average of an angular ratio.

**Belle projection and units.** [The appendix tables](https://arxiv.org/html/1909.01857v1#A2.T8) report percent-valued measurements. Intake converts both values and errors by 0.01. It separates 1D and 2D projections and UL/UC observables. The charged-pion 2D table includes low-z bins absent from some 1D projections; the parser retains the actual published bounds rather than inferring them from a caption or a different projection.

**STAR source intervals and sample overlap.** The HEPData table descriptions identify horizontal boxes as systematic uncertainties in the measured coordinate. Intake stores their `low/high` under `source_intervals`, leaves physical `bins` empty, and does not invent a mean from the interval midpoint. Source qualifiers, table DOIs, licenses, row locators and raw coordinate records are retained. STAR200 separate-year and combined-year measurements share an overlap group, as do alternative projections.

**STAR510 different published binning.** [The supplemental Table I](https://arxiv.org/html/2507.16355v2#as1_S5.T1) has 14 positive-pion jet-pT entries, while the corresponding HEPData Figure 2 view has 13 xT entries. They cannot be treated as a 14-point one-to-one transcription check or concatenated as independent measurements. The supplement's three named systematic components remain separate with unresolved correlations. No bin edges are fabricated from neighboring means. Physical binning, acceptance, polarization scales and matching process kernels still require review.

## Admission and integrity

Every candidate has `scope=experimental_candidate`, source hashes, source-cell locators, an overlap group and explicit production blockers. The normal fit-data validator and Julia numerical loader reject this scope. Intake can validate transcription and representation without passing G1 or any production physics gate.

Successful cached runs verify the locked raw bytes and parser/configuration fingerprint, then reconstruct every view and source-derived decision from those inputs. They compare canonical JSON digests of the complete candidate, including values/errors, means, bin records, harmonic conventions, source locators, overlap groups and admission flags. This preserves scalar types: numeric `0` cannot substitute for `false`, or `1` for `true`. Recomputing an artifact hash after changing a value cannot bypass this check. Report counts, primary-view choices, HERMES/BESIII cross-check decisions, lattice metadata and STAR interval/overlap statements must also match the source replay. The recorded completion time remains event metadata, not a source-derived scientific decision.

Both fresh intake and cached verification recheck the manifest, parser and raw-input snapshot after parsing. Verification also checks that cached artifacts stayed unchanged during replay. This is a complete data-preparation job; it adds no Python callback or language-model call to Julia numerical loops. All existing data adapters, values, uncertainty prescriptions, cuts and overlap policy are unchanged.

Linux CI runs `scripts/check_cached_experimental.py` against the real locked sources. Its nine negative cases cover ordinary artifact changes, rehashed values/harmonics/overlap labels, altered report counts and source decisions, scalar-type aliases, and a false physical-extraction claim. Original candidate/report bytes are restored and reverified before the software-integrity receipt is written. That receipt has scope `candidate_cache_integrity`; it does not admit an experimental candidate to fitting or replace any production gate.

Tests cover unit conversion, ambiguous numbers, MathML duplication, merged cells, source tampering, explicit HERMES corrections, exchanged BESIII bins, safe in-memory archive reading, plot-error semantics, projection overlap and rejection by the fit path. All test arrays are synthetic. Raw downloads, candidate arrays and reports remain local.
