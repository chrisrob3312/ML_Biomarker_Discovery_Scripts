# 11_final_models.R
# After the benchmark, pick the winning learner per outcome (or use the
# super learner) and refit on the full dataset with the locked-in
# hyperparameters from the outer-best inner tuner. Save model objects.

source(here::here("R", "00_config.R"))
source(here::here("R", "07_task_definitions.R"))
source(here::here("R", "08_learners.R"))
source(here::here("R", "09_resampling_nestedcv.R"))

fit_final <- function(task, learner, task_type = c("surv", "classif")) {
  task_type <- match.arg(task_type)
  at <- auto_tune(learner, task_type = task_type)
  at$train(task)
  at
}

fit_baseline_cox <- function(dat, outcome_name) {
  o <- CONFIG$outcomes[[outcome_name]]
  feats <- resolve_features(dat, outcome_name)

  # Strip the combined nested factor used by ML learners and rebuild it
  # explicitly as era + era:protocol_arm so the marginal era effect stays
  # identified for patients with protocol_arm = "unknown".
  feats <- setdiff(feats, "protocol_nested")
  spline_now <- intersect(CONFIG$spline_vars, feats)

  rhs <- c(setdiff(feats, spline_now),
           sprintf("rms::rcs(%s, 4)", spline_now),
           "treatment_era",
           "treatment_era:protocol_arm_clean")
  # Deduplicate in case treatment_era was already in feats
  rhs <- unique(rhs)

  f <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s",
                                 o$time_col, o$event_col,
                                 paste(rhs, collapse = " + ")))
  rms::cph(f, data = dat, x = TRUE, y = TRUE, surv = TRUE, time.inc = 12)
}

if (sys.nframe() == 0) {
  dat <- readRDS(file.path(PATHS$data_int, "engineered.rds"))
  tasks <- build_all_tasks(dat)

  # Example: refit elastic-net Cox and RSF + baseline Cox-RCS for each outcome.
  out <- list()
  out$relapse_glmnet <- fit_final(tasks$relapse,
                                  make_classif_learners()$glmnet, "classif")
  out$efs_glmnet     <- fit_final(tasks$efs,
                                  make_surv_learners()$glmnet, "surv")
  out$os_glmnet      <- fit_final(tasks$os,
                                  make_surv_learners()$glmnet, "surv")
  out$efs_cox_rcs    <- fit_baseline_cox(dat, "efs")
  out$os_cox_rcs     <- fit_baseline_cox(dat, "os")

  saveRDS(out, file.path(PATHS$results, "final_models.rds"))
}
