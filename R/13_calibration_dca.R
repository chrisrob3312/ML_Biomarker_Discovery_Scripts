# 13_calibration_dca.R
# Calibration, time-dependent AUC, and decision-curve analysis.

source(here::here("R", "00_config.R"))

# Calibration intercept + slope from Cox (cv-predicted linear predictors)
cox_calibration <- function(linpred, time, event, horizon) {
  surv_obj <- survival::Surv(time, event)
  fit <- survival::coxph(surv_obj ~ linpred)
  list(slope = unname(stats::coef(fit)),
       intercept = NA_real_,  # intercept defined via calibration plot at horizon
       coxph_fit = fit)
}

# Time-dependent AUC at the configured horizons
td_auc <- function(score, time, event, horizons = CONFIG$horizons_months) {
  timeROC::timeROC(T = time, delta = event,
                   marker = score, cause = 1,
                   times = horizons, iid = TRUE)
}

# Integrated Brier
ibs_score <- function(pred_surv, train_surv, eval_times) {
  riskRegression::Score(list(model = pred_surv),
                        formula = train_surv,
                        times = eval_times,
                        metrics = c("brier"),
                        summary = "ibs",
                        contrasts = FALSE)
}

# Decision curve analysis - binary outcome (relapse) at fixed horizon
dca_binary <- function(dat, outcome, model_probs, thresholds = seq(0.05, 0.5, 0.05)) {
  d <- dat
  d$.pred <- model_probs
  d$.y    <- dat[[outcome]]
  dcurves::dca(stats::as.formula(".y ~ .pred"), data = d,
               thresholds = thresholds)
}

# Decision curve analysis - survival outcome at horizon t
dca_survival <- function(dat, time_col, event_col, model_risk, horizon,
                         thresholds = seq(0.05, 0.5, 0.05)) {
  d <- dat
  d$.risk <- model_risk
  d$.surv <- survival::Surv(d[[time_col]], d[[event_col]])
  dcurves::dca(stats::as.formula(".surv ~ .risk"),
               data = d, time = horizon, thresholds = thresholds)
}
