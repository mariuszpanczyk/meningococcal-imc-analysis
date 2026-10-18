# Input dictionary (no participant records)

The full descriptive cohort comprises 53 children. Two deaths within 24 hours of admission were excluded from IMC analyses because IMC could not be assessed. The analytical database therefore contains 51 eligible children (13 with IMC and 38 without IMC). This is not a missing-data discrepancy.

One row per participant. CSV: UTF-8, semicolon delimiter. Numeric values retain source units; no unit rescaling is applied. Missing cells remain missing. Original column names are retained.

## Required analysis fields
- powiklanie_immuno: tak/nie; mapped to 1/0. Must be observed in the 51-person IMC analytical population. Two early deaths in the 53-person descriptive cohort are excluded from IMC analyses and must not be assigned 0/no IMC.
- plec: K/M; reference K, comparison M.
- wiek_lata: age in years.
- objawy_do_hosp_dni: symptom-onset to admission, days.
- zomr: nie / nie badano / tak; not-assessed retained as an explicit comparison category.
- oitd: nie/tak; PICU admission at any time, treated as a secondary severity marker.
- serotyp: B/C/W/Y; missing serogroup excluded per test.
- crp_max: mg/dL; numeric logistic OR per doubling.
- pct_max: ng/mL; numeric logistic OR per doubling.
- wbc_max, wbc_min, limf_min, plt_min: 10^3/uL; OR per available-case IQR.
- mleczany_max: mmol/L; OR per doubling.
- hco3_min: mmol/L; OR per IQR.
- inr_max: unitless; OR per IQR.
- fibrynogen_min: g/L; OR per IQR.
- d_dimer_max: ng/mL; OR per doubling.

## Antibiotic delay
Provide data_przyjecia and data_atb (DD.MM.YYYY), or explicit przyjecie_do_atb_dni (days). If both are supplied they must agree. The source oczekiwanie_atb_dni is audited against onset-to-antibiotic delay and is NOT substituted for admission-to-antibiotic delay.

## Optional private audit/descriptive fields
- id: private unique identifier, checked but not retained in the analytical frame or exported.
- data_urodzenia, data_objawow, data_przyjecia, data_atb, data_wypisu: private dates in DD.MM.YYYY; chronology audit only.
- czas_hosp_dni, dlugosc_atb_dni: descriptive durations only, excluded as predictors.
- steryd_przed_powiklaniem: not used as a predictor or exported; temporal ambiguity.

## Optional complete-cohort input
The separate full-cohort file follows the same schema and adds death_within_24h (tak/nie). IMC must be missing for early deaths and observed for assessable children. data_objawow is required for annual counts. No missing patients are inferred. The optional module was not validated against an unavailable complete 53-person input.

## Families and scaling
Clinical: age, onset-to-admission, admission-to-antibiotic, sex, meningitis, serogroup (6).
Laboratory/severity: 11 laboratory markers plus PICU (12).
BH is applied separately to each family and separately for permutation and association analyses. Firth models are post hoc nominal inference. No imputation. Zero IQR uses a one-source-unit numeric ML contrast with an explicit label; nonpositive log2 values stop execution.
