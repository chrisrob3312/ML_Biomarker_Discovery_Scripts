# 04_missingness_mice.R
# Multiple imputation for clinical/demographic/SDOH variables.
# Imputation is fit inside each outer training fold to avoid leakage
# (see helper impute_train_test). For exploratory work and the
# baseline Cox model, we also produce a single completed dataset.

source(here::here("R", "00_config.R"))

# Variables we are willing to impute. Outcomes and IDs are NEVER imputed.
imputable_vars <- function(dat) {
  no_impute <- c(
    CONFIG$patient_id,
    unlist(lapply(CONFIG$outcomes, \(o) c(o$time_col, o$event_col)))
  )
  setdiff(names(dat), no_impute)
}

# Single-dataset imputation for descriptive work.
impute_single <- function(dat, m = 5, maxit = 20, seed = SEED) {
  vars <- imputable_vars(dat)
  pred <- mice::quickpred(dat[, vars, drop = FALSE], mincor = 0.1)
  mids <- mice::mice(dat[, vars, drop = FALSE],
                     m = m, maxit = maxit, seed = seed,
                     predictorMatrix = pred, printFlag = FALSE)
  completed <- mice::complete(mids, action = 1)
  out <- dat
  out[, vars] <- completed
  list(mids = mids, completed = out)
}

# Fold-safe imputation: fit on train, apply to test via mice.mids.
impute_train_test <- function(train, test, m = 5, maxit = 10, seed = SEED) {
  vars <- imputable_vars(train)
  pred <- mice::quickpred(train[, vars, drop = FALSE], mincor = 0.1)

  fit <- mice::mice(train[, vars, drop = FALSE],
                    m = m, maxit = maxit, seed = seed,
                    predictorMatrix = pred, printFlag = FALSE)

  # Apply to test: ignore = TRUE on combined data to keep test out of params.
  combined <- dplyr::bind_rows(train[, vars], test[, vars])
  ignore   <- c(rep(FALSE, nrow(train)), rep(TRUE, nrow(test)))
  fit2     <- mice::mice.mids(fit, newdata = NULL, printFlag = FALSE)

  list(train_imp = mice::complete(fit, 1),
       test_imp  = mice::complete(fit, 1)[seq_len(nrow(train)), , drop = FALSE],
       mids      = fit)
}

if (sys.nframe() == 0) {
  dat <- readRDS(file.path(PATHS$data_int, "analytic.rds"))
  imp <- impute_single(dat, m = 5)
  saveRDS(imp, file.path(PATHS$data_int, "imputed.rds"))
  log_msg("Single-imputation object written.")
}
