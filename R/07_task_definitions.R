# 07_task_definitions.R
# Build mlr3 Tasks for: relapse (binary), EFS (survival), OS (survival),
# in both overall and MRD-stratified forms.

source(here::here("R", "00_config.R"))

# Columns to drop from features (IDs, raw times/events, redundant pre-engineering)
non_feature_cols <- function() {
  out <- c(CONFIG$patient_id, "mrd_eoi_continuous", "mrd_eoc_continuous",
           "wbc_dx")  # drop raw versions, keep log/cat
  outcome_cols <- unlist(lapply(CONFIG$outcomes,
                                \(o) c(o$time_col, o$event_col)))
  c(out, outcome_cols)
}

# Resolve the feature set for a given outcome: core set from features.yaml
# minus exclusions plus extras. Core = clinical + cytogenetics + mrd +
# ancestry_sdoh + ancestry-dosed risk allele columns present in the data.
# All survival models are EOI-anchored, so end-of-induction predictors
# (MRD, induction kinetics) are legitimate baseline covariates.
core_features <- function(dat) {
  cols <- unique(unlist(c(
    CONFIG$clinical,
    CONFIG$cytogenetics,
    CONFIG$mrd,
    CONFIG$ancestry_sdoh
  )))
  # Drop the two raw protocol columns - the nested factor below replaces
  # them in ML learners. Baseline Cox in 11 still reads them directly.
  cols <- setdiff(cols, c("protocol_arm", "treatment_era"))
  # Engineered companions
  cols <- c(cols, "wbc_dx_log10", "mrd_eoi_log10", "mrd_eoc_log10",
            "age_cat_clin", "wbc_cat_clin", "mrd_eoi_cat_clin",
            "protocol_nested")
  # Risk-allele + local-ancestry-dosed columns
  ra_cols <- character()
  for (ra in CONFIG$risk_alleles) {
    ra_cols <- c(ra_cols, paste0(ra$id, "_dosage"))
    for (anc in ra$ancestries)
      ra_cols <- c(ra_cols, sprintf("%s_la_%s_dose", ra$id, anc))
  }
  cols <- c(cols, ra_cols)
  intersect(cols, names(dat))
}

resolve_features <- function(dat, outcome_name) {
  core <- core_features(dat)
  ov   <- CONFIG$per_outcome[[outcome_name]]
  feats <- core
  if (!is.null(ov$exclude))       feats <- setdiff(feats, ov$exclude)
  if (!is.null(ov$include_extra)) feats <- union(feats, intersect(ov$include_extra, names(dat)))
  # Drop the outcome's own time/event columns if they snuck in
  o <- CONFIG$outcomes[[outcome_name]]
  feats <- setdiff(feats, c(o$time_col, o$event_col, CONFIG$patient_id))
  feats
}

build_task_relapse <- function(dat) {
  o <- CONFIG$outcomes$relapse
  d <- dat |>
    dplyr::filter(!is.na(.data[[o$event_col]])) |>
    dplyr::mutate(target = factor(.data[[o$event_col]],
                                  levels = c(0, 1),
                                  labels = c("no", "yes")))
  feats <- resolve_features(d, "relapse")
  mlr3::TaskClassif$new(
    id      = "relapse",
    backend = d[, c("target", feats)],
    target  = "target",
    positive = "yes"
  )
}

build_task_survival <- function(dat, outcome_name) {
  o <- CONFIG$outcomes[[outcome_name]]
  d <- dat |>
    dplyr::filter(!is.na(.data[[o$time_col]]),
                  !is.na(.data[[o$event_col]]),
                  .data[[o$time_col]] > 0)
  feats <- resolve_features(d, outcome_name)
  mlr3proba::TaskSurv$new(
    id      = outcome_name,
    backend = d[, c(o$time_col, o$event_col, feats)],
    time    = o$time_col,
    event   = o$event_col,
    type    = "right"
  )
}

build_all_tasks <- function(dat) {
  list(
    dfs     = build_task_survival(dat, "dfs"),
    os      = build_task_survival(dat, "os"),
    relapse = build_task_relapse(dat)
  )
}

# For MRD-stratified models: returns a list of tasks (one per stratum)
# for each outcome.
build_mrd_stratified <- function(dat) {
  strat_col <- CONFIG$mrd_strata$column
  levels    <- CONFIG$mrd_strata$levels
  out <- list()
  for (lvl in levels) {
    sub <- dat[dat[[strat_col]] %in% lvl, , drop = FALSE]
    if (nrow(sub) < 50) {
      log_msg(sprintf("Stratum %s has only %d rows - skipping.",
                      lvl, nrow(sub)))
      next
    }
    out[[lvl]] <- build_all_tasks(sub)
  }
  out
}

if (sys.nframe() == 0) {
  dat <- readRDS(file.path(PATHS$data_int, "engineered.rds"))
  tasks_overall <- build_all_tasks(dat)
  tasks_mrd     <- build_mrd_stratified(dat)
  saveRDS(list(overall = tasks_overall, mrd_strat = tasks_mrd),
          file.path(PATHS$data_int, "tasks.rds"))
}
