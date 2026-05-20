# 06_cutpoints.R
# (1) Apply pre-specified clinical cut points (already produced in 05).
# (2) Exploratory data-driven cuts via maxstat with bootstrap CI and
#     BH-FDR multiplicity correction. Reported alongside, NEVER used
#     to dichotomize variables that feed the ML models.

source(here::here("R", "00_config.R"))

# Run maxstat on each (variable, outcome) pair for survival outcomes.
# Returns a tibble of cut points, log-rank stats, bootstrap CIs, and BH q-values.
maxstat_scan <- function(dat,
                         vars,
                         time_col,
                         event_col,
                         reps = CONFIG$maxstat$bootstrap_reps,
                         fdr_method = CONFIG$maxstat$fdr_method,
                         seed = SEED) {
  set.seed(seed)
  surv <- survival::Surv(dat[[time_col]], dat[[event_col]])

  results <- purrr::map_dfr(vars, function(v) {
    x <- dat[[v]]
    if (!is.numeric(x)) return(NULL)
    ok <- is.finite(x) & !is.na(surv)
    if (sum(ok) < 30) return(NULL)

    fit <- try(maxstat::maxstat.test(
      surv[ok] ~ x[ok], smethod = "LogRank", pmethod = "Lau94"
    ), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)

    boots <- replicate(reps, {
      idx <- sample(which(ok), replace = TRUE)
      f <- try(maxstat::maxstat.test(
        surv[idx] ~ x[idx], smethod = "LogRank", pmethod = "Lau94"
      ), silent = TRUE)
      if (inherits(f, "try-error")) NA_real_ else f$estimate
    })

    tibble::tibble(
      variable    = v,
      cutpoint    = unname(fit$estimate),
      stat        = unname(fit$statistic),
      p_value     = fit$p.value,
      ci_lo       = stats::quantile(boots, 0.025, na.rm = TRUE),
      ci_hi       = stats::quantile(boots, 0.975, na.rm = TRUE),
      n_used      = sum(ok)
    )
  })

  if (nrow(results) == 0) return(results)
  results$q_value <- stats::p.adjust(results$p_value, method = fdr_method)
  dplyr::arrange(results, q_value)
}

# Apply a cut and return a 2-level factor; for QC / reporting only.
apply_cut <- function(x, cut, labels = c("low", "high")) {
  factor(ifelse(x < cut, labels[1], labels[2]), levels = labels)
}

if (sys.nframe() == 0) {
  dat <- readRDS(file.path(PATHS$data_int, "engineered.rds"))
  cont_vars <- intersect(CONFIG$spline_vars, names(dat))
  for (oc in names(CONFIG$outcomes)) {
    o <- CONFIG$outcomes[[oc]]
    if (o$type != "survival") next
    log_msg("Maxstat scan for outcome: ", oc)
    res <- maxstat_scan(dat, cont_vars, o$time_col, o$event_col)
    readr::write_csv(res,
                     file.path(PATHS$results,
                               sprintf("maxstat_%s.csv", oc)))
  }
}
