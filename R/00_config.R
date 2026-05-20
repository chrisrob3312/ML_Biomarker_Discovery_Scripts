# 00_config.R
# Central configuration: paths, seeds, options, package loading.
# Source this at the top of every script.

suppressPackageStartupMessages({
  library(yaml)
  library(here)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(forcats)
  library(survival)
  library(rms)
  library(glmnet)
  library(mice)
  library(maxstat)
  library(survminer)
  library(mlr3)
  library(mlr3proba)
  library(mlr3learners)
  library(mlr3extralearners)
  library(mlr3tuning)
  library(mlr3pipelines)
  library(mlr3filters)
  library(mlr3misc)
  library(paradox)
  library(ranger)
  library(xgboost)
  library(randomForestSRC)
  library(riskRegression)
  library(timeROC)
  library(dcurves)
  library(fastshap)
  library(ggplot2)
})

SEED <- 20260520
set.seed(SEED)

PATHS <- list(
  root      = here::here(),
  config    = here::here("config"),
  data_raw  = here::here("data", "raw"),
  data_int  = here::here("data", "interim"),
  data_proc = here::here("data", "processed"),
  results   = here::here("results"),
  reports   = here::here("reports")
)
invisible(lapply(PATHS[c("data_raw","data_int","data_proc","results","reports")],
                 dir.create, recursive = TRUE, showWarnings = FALSE))

CONFIG <- yaml::read_yaml(file.path(PATHS$config, "features.yaml"))

# Resampling defaults
RESAMPLING <- list(
  outer_folds = 10,
  inner_folds = 5,
  repeats     = 1,
  stratify    = TRUE
)

# Tuning budget (kept modest by default; raise for final runs)
TUNING <- list(
  n_evals = 30,
  method  = "random_search"
)

# Logger
log_msg <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}
