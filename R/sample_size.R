# sample_size.R
# Riley et al. (2020) sample-size calculation for clinical prediction
# models via the pmsampsize package. Reports minimum N for each outcome,
# given anticipated event rate, follow-up, R^2, and candidate parameter count.
#
# Inputs you (the analyst) supply:
#   anticipated_R2 : Cox-Snell R^2 you expect from a "good" model.
#       Pediatric B-ALL DFS with MRD + cytogenetics + clinical: 0.10-0.18
#       OS is harder: 0.06-0.12
#   prevalence     : event proportion at the time horizon of interest
#   timepoint      : horizon in months matching prevalence
#   mean_followup  : mean follow-up time in months
#   parameters     : number of candidate parameters AFTER one-hot expansion
#                    and RCS basis expansion. Count carefully.

source(here::here("R", "00_config.R"))
suppressPackageStartupMessages({
  library(pmsampsize)
})

# Count candidate parameters for the Cox-RCS baseline model. This is the
# conservative count that drives the formula; ML models with penalization
# tolerate more parameters but the baseline number is what you report.
count_candidate_parameters <- function(dat, feats,
                                       spline_vars = CONFIG$spline_vars,
                                       knots = 4) {
  n <- 0L
  for (v in feats) {
    if (!v %in% names(dat)) next
    if (v %in% spline_vars) {
      n <- n + (knots - 1L)               # RCS basis: knots-1 columns
    } else if (is.factor(dat[[v]]) || is.character(dat[[v]])) {
      lv <- length(unique(stats::na.omit(dat[[v]])))
      n <- n + max(0L, lv - 1L)           # one-hot minus reference
    } else {
      n <- n + 1L
    }
  }
  n
}

sample_size_survival <- function(prevalence,
                                 timepoint,
                                 mean_followup,
                                 parameters,
                                 anticipated_R2 = 0.10) {
  pmsampsize::pmsampsize(
    type = "s",
    rsquared      = anticipated_R2,
    parameters    = parameters,
    rate          = -log(1 - prevalence) / timepoint,  # constant-hazard approx
    timepoint     = timepoint,
    meanfup       = mean_followup
  )
}

sample_size_binary <- function(prevalence,
                               parameters,
                               anticipated_R2 = 0.10) {
  pmsampsize::pmsampsize(
    type = "b",
    rsquared      = anticipated_R2,
    parameters    = parameters,
    prevalence    = prevalence
  )
}

# Convenience: run all three outcomes given realized rates from the data.
run_all_sample_size <- function(dat, feats,
                                horizon_months = 60,
                                anticipated_R2 = list(dfs = 0.12, os = 0.08,
                                                      relapse = 0.10)) {
  params <- count_candidate_parameters(dat, feats)
  log_msg("Candidate parameters (incl. RCS expansion): ", params)

  # Estimate prevalence and mean follow-up empirically when available.
  km_event_rate <- function(time, event, t) {
    fit <- try(survival::survfit(survival::Surv(time, event) ~ 1),
               silent = TRUE)
    if (inherits(fit, "try-error")) return(NA)
    s <- summary(fit, times = t, extend = TRUE)
    1 - s$surv
  }

  out <- list()
  if (all(c("dfs_time_months", "dfs_event") %in% names(dat))) {
    p <- km_event_rate(dat$dfs_time_months, dat$dfs_event, horizon_months)
    mfu <- mean(dat$dfs_time_months, na.rm = TRUE)
    out$dfs <- sample_size_survival(p, horizon_months, mfu, params,
                                    anticipated_R2$dfs)
  }
  if (all(c("os_eoi_time_months", "os_eoi_event") %in% names(dat))) {
    p <- km_event_rate(dat$os_eoi_time_months, dat$os_eoi_event, horizon_months)
    mfu <- mean(dat$os_eoi_time_months, na.rm = TRUE)
    out$os <- sample_size_survival(p, horizon_months, mfu, params,
                                   anticipated_R2$os)
  }
  if ("relapse_event" %in% names(dat)) {
    p <- mean(dat$relapse_event == 1, na.rm = TRUE)
    out$relapse <- sample_size_binary(p, params, anticipated_R2$relapse)
  }
  out
}
