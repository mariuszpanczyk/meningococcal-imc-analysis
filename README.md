# IMD / IMC analytical code, v0.2.2

Code-only companion for the study of immune-mediated complications after invasive meningococcal disease. English documentation; original Polish input column names retained. No raw database, participant-level analysis data, synthetic patient data, manuscript attachments, original HTML or saved research outputs are included.

## Status
A reconstructed and privately tested implementation, not a recovered historical script. The supplied HTML contains results only. Historical random seeds were not available. All Monte Carlo analyses use a newly declared seed. This is a release candidate until the documented manuscript discrepancy, citation metadata and software license are approved.

## Files
- IMD_IMC_analysis_code_only_v0_2_2.ipynb: self-contained Colab/Jupyter notebook.
- imd_imc_analysis.py: matching Python code companion; defaults to no-data mode.
- requirements.txt: pinned reference environment selected for this implementation.
- DATA_DICTIONARY.md: input schema without records.
- VALIDATION_REPORT.md: scope and limitations of performed checks.
- RELEASE_CHECKLIST.md: final scientific, privacy and publication checks.
- CITATION.cff: draft metadata; complete the approved full author list.
- PUBLICATION_METADATA_TEMPLATE.json: manual metadata worksheet, not an automatic .zenodo.json file.
- .gitignore: convenience exclusions only, not a privacy guarantee.

## Colab use
1. Open the .ipynb file in Colab and run all cells. Public defaults execute unit tests only and load no study data.
2. For an authorized private run, work on a private copy. Confirm the approved environment and cohort eligibility, set RUN_PRIVATE_ANALYSIS=True, and set DATA_PATH. For interactive Colab upload, leave DATA_PATH=None and set ALLOW_COLAB_UPLOAD=True.
3. The input is a semicolon-delimited UTF-8 CSV. The source column names and categories are documented in DATA_DICTIONARY.md. No public download link is used.
4. Results are aggregate tables in RESULTS. Optional private export is disabled by default. No participant rows, residual vectors or fitted probabilities are exported.
5. Do not publish the executed private copy. Restore public switches, clear outputs and manually inspect all source text and repository history.

## Analyses
Shared-label permutation tests (19,999), BH in clinical/laboratory-severity families, supplementary single-step WY minP; rank-biserial/HL effects and stratified bootstrap CIs (4,999); exact conditional binary ORs; numeric univariable ML models, multilevel categorical LR tests, descriptive adjusted models; Firth platelet models with penalized-profile inference. The Firth implementation is native Python and does not require R.

## Study populations and remaining source issues
- **Population definition, not a discrepancy:** the full descriptive cohort includes 53 children. Two died within 24 hours of admission and were excluded from IMC analyses because IMC could not be assessed. The IMC analytical population includes 51 children (13 with IMC, 38 without IMC). The supplied analytical database is appropriate for that population; the difference is not missing data or a database error. Early deaths must not be classified as no IMC. The optional full-cohort module requires a separate authorized descriptive input. Narrative variables absent from the analytical input cannot be reproduced.
- Admission-to-antibiotic permutation p is 1.000 in the HTML but .754 in Table 2. Recomputed results are not forced to the table. Reconcile this discrepancy before final release.
- The input has no laboratory measurement/event timestamps; pre-event ordering and the manuscript laboratory window cannot be independently verified.
- Historical random streams are unknown, so bootstrap CIs and permutation p values are not guaranteed to match the last reported decimal. The WY implementation is explicitly reconstructed and its assumptions are documented.

## GitHub / Zenodo
Extract this ZIP and review its contents. Upload only the reviewed code/documentation, not the input CSV or private working directory. Complete citation metadata and add the approved software license. Enable the public repository in Zenodo's GitHub integration and create a versioned GitHub release; verify the resulting archive and DOI. No repository, deposit or DOI has been created by this package.

- [Zenodo release archiving](https://help.zenodo.org/docs/github/archive-software/github-upload/)
- [GitHub archiving and citation](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content)

## Data availability
Participant-level data are not included. Exact numerical reproduction requires the original analysis data. Authors must state actual access conditions and reasons for restriction separately. No on-request access promise is made. Code sharing and data sharing are separate decisions.

## Author-reported environment
The authors report that analyses were performed in Google Colab using Python 3.13.16, pandas 2.2.3, NumPy 2.1.3, SciPy 1.16.3 and statsmodels 0.15.0. These versions are recorded as supplied and have not been independently installed or verified in this delivery.

The optional installation cell and `requirements.txt` retain the previously selected reference environment used for this reconstruction; they do not represent the author-reported study environment. Actual prior validation versions remain in `VALIDATION_REPORT.md`. Do not describe the reconstruction as tested under the author-reported environment until a separate run verifies it. See `AUTHOR_REPORTED_ENVIRONMENT.json`.

## Updated citation metadata
CITATION.cff now includes the supplied title and all six authors in the supplied order, with affiliations and the supplied contact email. The record remains type software; it does not assert that the manuscript is already published. Affiliation wording is reproduced as supplied.
