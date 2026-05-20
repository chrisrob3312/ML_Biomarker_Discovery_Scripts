# 10_benchmark.R
# Benchmark the learner roster across outcomes and MRD strata via nested CV.

source(here::here("R", "00_config.R"))
source(here::here("R", "07_task_definitions.R"))
source(here::here("R", "08_learners.R"))
source(here::here("R", "09_resampling_nestedcv.R"))

surv_measures <- function() {
  list(
    mlr3::msr("surv.cindex"),
    mlr3::msr("surv.graf"),     # Integrated Brier
    mlr3::msr("surv.rcll"),
    mlr3::msr("surv.dcalib")
  )
}

classif_measures <- function() {
  list(
    mlr3::msr("classif.auc"),
    mlr3::msr("classif.prauc"),
    mlr3::msr("classif.bbrier"),
    mlr3::msr("classif.logloss")
  )
}

run_benchmark <- function(task, kind = c("surv", "classif"), include_super = TRUE) {
  kind <- match.arg(kind)
  base_learners <- if (kind == "surv") make_surv_learners() else make_classif_learners()
  tuned <- lapply(base_learners, auto_tune, task_type = kind)
  if (include_super) {
    super <- if (kind == "surv") make_super_learner_surv()
             else make_super_learner_classif()
    tuned[["super_learner"]] <- super
  }
  outer <- make_outer_resampling(task)
  design <- mlr3::benchmark_grid(tasks = task, learners = tuned, resamplings = outer)
  bmr <- mlr3::benchmark(design, store_models = FALSE)
  bmr
}

run_all <- function(dat, out_dir = file.path(PATHS$results, "benchmark")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tasks <- build_all_tasks(dat)
  results <- list()

  results$dfs     <- run_benchmark(tasks$dfs,     "surv")
  results$os      <- run_benchmark(tasks$os,      "surv")
  results$relapse <- run_benchmark(tasks$relapse, "classif")

  # MRD-stratified
  strat <- build_mrd_stratified(dat)
  for (lvl in names(strat)) {
    results[[paste0("dfs_",     lvl)]] <- run_benchmark(strat[[lvl]]$dfs,     "surv")
    results[[paste0("os_",      lvl)]] <- run_benchmark(strat[[lvl]]$os,      "surv")
    results[[paste0("relapse_", lvl)]] <- run_benchmark(strat[[lvl]]$relapse, "classif")
  }

  saveRDS(results, file.path(out_dir, "benchmark_results.rds"))

  # Aggregate table
  agg <- purrr::map_dfr(names(results), function(nm) {
    bmr <- results[[nm]]
    msrs <- if (inherits(bmr$tasks$task[[1]], "TaskSurv")) surv_measures() else classif_measures()
    tab <- bmr$aggregate(msrs)
    tab$slice <- nm
    tab
  })
  readr::write_csv(agg, file.path(out_dir, "benchmark_summary.csv"))
  log_msg("Benchmark complete: ", file.path(out_dir, "benchmark_summary.csv"))
  invisible(results)
}

if (sys.nframe() == 0) {
  dat <- readRDS(file.path(PATHS$data_int, "engineered.rds"))
  run_all(dat)
}
