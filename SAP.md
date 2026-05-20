# Statistical Analysis Plan (SAP)
**Project**: ML prediction of DFS and OS in Southwestern US pediatric/AYA B-ALL
**Version**: 0.1 (draft)
**Cohort lock**: pending data inventory (Step 2)
**Branch**: `claude/ml-clinical-prediction-models-aQE3m`

## 1. Background and aims
Develop and internally validate clinically-interpretable ML models for
disease-free survival (DFS) and overall survival (OS) from end of induction
(EOI) in a retrospective, ethnically diverse (~50% Hispanic/Latino) cohort
of pediatric/AYA B-ALL patients treated in the Southwestern United States
from 2010 onward. Compare a benchmark roster of standard and
state-of-the-art learners and report the most salient features per outcome.

## 2. Cohort
- **Source**: retrospective records, institution(s) to be specified.
- **Diagnoses**: B-cell acute lymphoblastic leukemia, de novo.
- **Dx window**: 2010-01-01 through 2025-12-31.
- **Age at dx**: 0 to 21 years inclusive.
- **Inclusion**: confirmed B-ALL diagnosis; reached end of induction
  with documented response assessment.
- **Exclusions**:
  - T-cell ALL or mixed-phenotype acute leukemia
  - Down syndrome ALL
  - Relapsed disease at presentation
  - Induction failure or induction death (cohort defined at EOI landmark)
  - Insufficient baseline data to compute the candidate predictor set
- **Time origin (t = 0)**: end of induction (EOI). All survival times
  are measured from EOI.

## 3. Outcomes
| Outcome | Type | Event | Time origin |
|---|---|---|---|
| DFS (primary) | Right-censored survival | relapse OR second malignancy OR death any cause | EOI |
| OS (co-primary) | Right-censored survival | death any cause | EOI |
| Relapse (secondary) | Binary | ever-relapsed by last follow-up | n/a |

- Censoring: last documented contact alive without event.
- Minimum follow-up: patients alive and event-free with < 6 months
  follow-up from EOI are included in the primary analysis (administrative
  censoring) and excluded in a sensitivity analysis.

## 4. Predictors (candidate set, knowable by EOI)
Grouped by acquisition timing. ML models receive the union; the baseline
Cox model uses pre-specified clinically meaningful variables.

- **Demographics**: age at dx, sex, race/ethnicity, global ancestry
  proportions (EUR/AFR/AMR/EAS/SAS), ADI national & state deciles,
  insurance payer, dx_year, treatment era, protocol arm (nested within era).
- **Disease biology (at dx)**: WBC, CNS status, BMI z, cyto subtype,
  binary biology flags (Ph+, Ph-like, KMT2A-r, ETV6-RUNX1, TCF3-PBX1,
  hyperdiploid_high, hypodiploid_low, iAMP21, DUX4-r, MEF2D-r, ZNF384-r,
  PAX5alt, TP53 alt, IKZF1plus, CDKN2A/B del, JAK pathway alt,
  RAS pathway alt), CRLF2 rearrangement, NCI risk.
- **Risk alleles + local ancestry dosing**: ARID5B rs7089424, IKZF1
  rs4132601, GATA3 rs3824662, CDKN2A rs3731217, PIP4K2A rs7088318 —
  each as total dosage and ancestry-attributed dosage.
- **Induction kinetics**: day-8 PB blasts, day-15 marrow blasts.
- **End of induction**: MRD-EOI (continuous + clinical category),
  clinical_trial_enrolled, treatment_type.
- **Omics** (optional, screened inside CV): bulk RNA, tumor methylation,
  metabolomics, germline methylation. Held back for initial models.

Continuous variables modeled with restricted cubic splines (3-4 knots)
in baseline Cox / logistic. Tree learners use raw values.

## 5. Sample size
Target N: 1700-1900 with germline genome ascertainment.
Anticipated DFS event rate (5y): 20-25% (Hispanic-enriched cohort with
higher Ph-like incidence).
Anticipated OS event rate (5y): 8-12%.

Formal calculation in Step 2 via `pmsampsize` with:
- Outcome event rate (above)
- Mean follow-up (TBD from data)
- Anticipated Cox-Snell R^2 = 0.10-0.15
- Target shrinkage 0.9, MAPE <= 0.05

## 6. Missing data
- Multiple imputation by chained equations (`mice`, m = 20) for the
  primary analysis. Outcomes never imputed.
- Auxiliary variables: ancestry proportions, treatment era, dx_year.
- Imputation fit on training fold only inside nested CV (no leakage).
- Sensitivity: complete-case analysis.

## 7. Model development
- **Baseline interpretable**: Cox-RCS for DFS and OS; logistic-RCS for
  relapse. Pre-specified covariates only (no data-driven selection).
  Treatment_era + treatment_era:protocol_arm nested encoding.
- **Penalized linear**: elastic-net Cox / logistic (`glmnet`), alpha
  and lambda tuned by inner CV.
- **Tree ensembles**: ranger RSF, randomForestSRC RSF, AORSF.
- **Gradient boosting**: xgboost (Cox / logistic).
- **Stack**: SuperLearner-style mlr3pipelines stack with glmnet meta.
- **Reference**: Kaplan-Meier (survival) / featureless (classif).

## 8. Resampling and tuning
- Nested CV: outer 10-fold stratified (by event for survival, class for
  binary), inner 5-fold. Repeated outer CV is optional.
- Tuner: random search, 30 evaluations per learner (raise for final run).
- Any feature screening (univariate Cox, omics top-K) executed inside
  inner folds.

## 9. Performance metrics (all reported)
- **Discrimination**: Harrell's C, Uno's C, time-dependent AUC at 1y/3y/5y
  with bootstrap 95% CIs.
- **Calibration**: calibration plot at 3y and 5y, calibration intercept &
  slope, integrated calibration index (ICI).
- **Clinical utility**: decision curve analysis (net benefit) at
  thresholds 5-50%.
- **Survival overall**: Integrated Brier Score.

## 10. Winner-selection rule
Primary criterion: superior decision-curve net benefit at a clinically
relevant 3y DFS threshold (e.g. 15-25% predicted risk), with calibration
slope in [0.85, 1.15] and ICI <= 0.05. Discrimination is a tiebreaker,
not the primary criterion.

## 11. Internal validation
- Optimism-corrected bootstrap (500 reps) of all performance metrics for
  the winning learner.
- Stability: top-K feature selection probability across bootstrap
  resamples.

## 12. External / temporal validation
If a second institution cohort is unavailable, hold out dx 2022-2025 as
temporal validation; develop on 2010-2021.

## 13. Fairness audit (pre-specified subgroups)
Report discrimination, calibration slope, ICI, and DCA net benefit in:
- Hispanic vs non-Hispanic (primary)
- NIH/OMB race categories
- Global AMR ancestry quartile
- ADI national quintile
- MRD-EOI stratum (negative / low-positive / high-positive)
- Sex
- Treatment era
- Age band (<1, 1-9, 10-15, 16-21)

Differential calibration slope outside [0.8, 1.2] flagged as a limitation.

## 14. Interpretability
- Forest plot of Cox HRs (baseline model) with 95% CIs.
- SHAP global beeswarm + per-feature dependence for the winning ML model.
- Permutation importance with correlated-feature blocks.
- Partial dependence / ALE plots for top features.
- Nomogram from baseline Cox.
- KM curves stratified by predicted-risk tertile, with log-rank test.

## 15. Reproducibility
- All code on GitHub (this repo).
- Dependencies pinned via `renv`.
- Random seed = 20260520.
- Trained final model objects archived (e.g. Zenodo) at submission.
- SAP version-controlled (this file); any deviation documented in a
  CHANGELOG section below before unblinding.

## 16. Deviations log
(none yet)
