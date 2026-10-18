# Validation report: v0.2.1

## Performed
- Notebook JSON validated with nbformat; all Python cells and companion script parsed.
- Exact delivered public-mode cells executed successfully; no study data loaded in public mode.
- Exact delivered private orchestration executed on the authorized supplied analysis CSV, with export disabled.
- Full 19,999 shared permutations and 4,999 stratified bootstrap resamples executed successfully.
- All 18 prespecified comparisons and 14 continuous effect modules completed.
- All 14 rank-biserial point estimates agreed with supplied Table 2 at its printed precision.
- Numeric platelet ML OR, Wald CI, p and BH q agreed with the Results text at printed precision.
- Unadjusted and age-adjusted Firth platelet ORs, profile CIs and p values agreed with the Results text at printed precision.
- Sex conditional exact OR/CI agreed with the supplied table at printed precision.
- Duplicate-ID rejection tested. Fixed-array software checks passed.
- Admission-to-antibiotic permutation discrepancy confirmed: recomputation agrees with HTML rather than the conflicting manuscript Table 2 entry.
- Final notebook contains no saved outputs, execution counts, attachments or participant rows. Public archive assembled using an explicit file allowlist.

## Not established
- Recovery of original executable source, historical random seed or identical Monte Carlo stream.
- Exact equality of every Monte Carlo p value or bootstrap CI to historical rounded values.
- Independent R logistf cross-validation; the new Python Firth implementation was compared to manuscript aggregate estimates.
- Full descriptive-cohort Table 1/Figure 1 reproduction was not tested because the supplied input is the correctly defined 51-person IMC analytical population. The full cohort comprises 53 children; two deaths within 24 hours of admission are intentionally excluded from IMC analyses. This is an eligibility distinction, not an incomplete analytical database. The optional full-cohort branch requires a separate authorized descriptive input.
- Verification of laboratory measurement timing or event ordering.
- Fresh Colab installation of the pinned reference environment. Actual delivery environment versions are recorded below; these are not claimed as historical versions.
- GitHub upload, Zenodo deposit, DOI assignment, final author-list approval or software-license selection.

## Delivery test environment
{
  "python": "3.12.15",
  "packages": {
    "numpy": "2.5.4",
    "pandas": "3.0.6",
    "scipy": "1.18.1",
    "statsmodels": "0.15.0",
    "matplotlib": "3.8.4",
    "nbformat": "5.11.1"
  }
}

## Documentation revision 0.2.1
All analytical code cells were compared with 0.2.0 and are identical. The Python companion is unchanged. Notebook schema, code syntax, empty outputs and archive allowlist were checked. The numerical checks above were performed for 0.2.0 and were not rerun for this documentation-only revision.

## Metadata revision 0.2.2
Author-supplied citation and environment metadata added. Analytical code cells are identical to 0.2.1. The author-reported Python/package versions are not the prior delivery validation environment; no new numerical validation or package installation was performed for this metadata revision.
