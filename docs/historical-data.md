# Historical KPSY candidate intake

This is the data preparation boundary for the historical benchmark in issues #3/#4. `historical-intake` acquires 26 checksum-locked XLSX transcriptions from JAM commit `331f65606b2e33a511510be361021b7cff8ea410`, an independently locked primary Hall A source archive and six primary COMPASS deuteron HEPData tables. It emits `historical_benchmark_candidate` views. `historical_fit_ready`, `production_ready` and `physical_extraction` remain false. Python and Julia fit loaders reject this candidate scope.

```sh
snakemake --snakefile workflow/historical-intake.smk --cores 1
python scripts/check_cached_historical.py config/historical-sources.json data/raw/historical data/normalized/historical
```

`python -m collinsflow.cli historical-intake` runs the equivalent complete ingestion job. `--offline` requires cached or verified seed bytes; `--seed-root` can seed a research checkout. Raw bytes stay under ignored `data/raw/historical/`, candidates and reports under ignored `data/normalized/historical/`. No legacy database program or compressed pickle `.info` file is executed or deserialized.

## Selection and source identity

The [KPSY PDF](https://arxiv.org/pdf/1505.05589v1), pages 17–20, specifies 20 SIDIS views/140 amplitudes and six SIA views/122 amplitudes. The source contract retains all 26 historical views, including overlapping projections and UL/UC samples. These 262 candidate amplitudes are not 262 independent updated-fit observations.

For SIDIS P-dependent views, the source prescription uses a strict reported-mean `P<0.8 GeV` cut. It removes 14 of the 154 raw SIDIS rows, retaining six P points per pion charge for HERMES and each COMPASS target. Other projections keep their published rows. Every exclusion retains its source row and reason. Thresholds and expected counts are declared before normalization; they are not fitted to a desired result.

Each view binds its source hash, publication, exact headers, experiment, pion species, target, harmonic and projection. Data are evaluated at reported means only after a separate prediction adapter is reviewed. Integrated SIA views retain that source identity even when a reported P mean is also available. A missing mean P is not fabricated from the product of mean z and mean qT. Published edges are retained separately from means and plot spreads.

## Explicit source findings

**HERMES reference mismatch.** Visual inspection of the original KPSY arXiv v1 PDF confirms that Table II cites [95], while its bibliography maps [95] to the Sivers paper 0906.3918. The same PDF's [7] identifies the [2010 Collins paper](https://arxiv.org/abs/1006.4221). The candidate source explicitly selects Collins-labeled HERMES2010 transcriptions and retains the citation mismatch. This does not recover the exact original author input files or certify all transcribed numerical values.

**Hall A legacy column remapping.** The pinned JAM SIDIS files 3001/3002 contain an extra numeric bin ID under `Q2`; later numeric fields are consequently mislabeled. The supplied labels would place Q² values in x and central asymmetries in a statistical-error field. The source contract declares a narrow remapping for these exact file hashes and headers. Every remapped kinematic value, neutron Collins central value, statistical error, experimental systematic and model systematic is checked exactly against the final [Hall A primary paper](https://arxiv.org/pdf/1106.0363v4), Tables I/II. The primary TeX member and its archive are independently hashed; commented obsolete tables are excluded before parsing. Raw XLSX bytes are unchanged.

The selected primary rows are the published **extracted neutron** results. The measured material was ³He; the paper's effective-polarization/dilution extraction must not be applied again. Experimental and model systematics remain separate. The original KPSY likelihood's inclusion of the model uncertainty is unresolved here.

**SIA source columns and units.** The source specifications declare percent-valued central values and uncertainties and convert them by 0.01. BaBar table1003 spells its second uncertainty column `syt_u`; the literal source name is explicitly recorded, with primary numerical/semantic review still pending. The angular column `S2/1+C2` is selected explicitly, while the publication's exact type of angular averaging remains a blocker. A mean of an angular ratio and a ratio of angular means cannot be silently interchanged.

**COMPASS deuteron primary checks.** [HEPData 48532 v1](https://www.hepdata.net/record/ins779473), Tables 1–6, supplies the ALL-pion results from 2003–2004. All five reported means, central values and errors agree exactly with the six pinned JAM deuteron views: 52 raw rows, of which 46 survive the unchanged historical P-projection cut. Intake requires that exact agreement for every row, including excluded rows. It verifies the table DOI, reaction/charge, headers and Collins group separately from the Sivers group. Literal primary projection-bin bounds are retained as `primary_bins`; no edges are inferred from means.

The primary tables explicitly identify these errors as statistical only. Their corresponding `error_u` values become `statistical_error`, with the original source-column label retained. The systematic uncertainty remains unquantified in the tables; an empty `systematic_components` object is accompanied by an explicit missing-systematic status and blocker. It does not mean zero systematic error or an approved covariance.

The [deuteron paper](https://arxiv.org/pdf/0802.2160v2), its Collins-angle definition and Eq.6, specify the source phase `phi_h+phi_S-pi` and divide the measured modulation by the transverse spin-transfer factor. This source convention is recorded separately from the pending mapping to the theory contract. Intake preserves the central-value sign and normalization exactly. The normalized `depolarization_removed` field remains unresolved while that mapping is unapproved; no fit adapter is enabled.

**COMPASS proton version and uncertainty limits.** The proton transcriptions' `error_u` composition remains unresolved. The [arXiv v1 preprint](https://arxiv.org/pdf/1408.4405v1), PDF page 8, gives a 3% scale uncertainty; the [final journal paper](https://mediatum.ub.tum.de/doc/1349896/document.pdf), PDF page 7, gives 5%. The contract binds both versions and records this difference without changing any proton errors or constructing a likelihood. The final paper, pages 5–6, also discusses acceptance folded into projected asymmetries and nonzero correlations between different projections of the same event sample. The historical exclusion of normalization errors and the updated-fit nuisance/covariance policy require separate decisions. No missing systematic covariance is replaced by zero.

## Integrity and remaining work

The report's status is `CANDIDATES_VALIDATED`, which refers to acquisition, representation, source correction and selection. It records all 33 raw-input hashes, parser/configuration fingerprint, candidate hash, eight primary Hall A crosschecks, 52 primary COMPASS deuteron row checks and 14 exclusions. The Snakemake terminal rule revalidates cached outputs and current sources; the cache check deliberately changes one candidate, requires rejection, restores it and checks validity again.

KPSY Eq.116 describes statistical and experimental systematic errors combined in quadrature with normalization uncertainties excluded. This intake preserves components and constructs no likelihood. Primary numerical checks for other transcriptions, error semantics, covariance/overlap choices, SIA acceptance/integration, nuclear treatment and the reconstructed historical prediction contract remain work. The original fit program and modified HOPPET remain unavailable. No χ² reproduction tolerance, historical fit, tensor charge or uncertainty-coverage result is supplied by a candidate report.
