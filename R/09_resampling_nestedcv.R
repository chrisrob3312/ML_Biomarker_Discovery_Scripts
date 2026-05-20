# 09_resampling_nestedcv.R
# Nested CV scheme:
#   outer: K-fold (stratified by event for survival, by class for binary)
#         used for honest performance estimation
#   inner: K-fold inside each outer training set, used by mlr3tuning
#         to select hyperparameters
# Omics screening (R/05) and any maxstat cuts must occur inside the inner loop.
# In mlr3 this is achieved by wrapping the screener as a PipeOp upstream of the
# learner so it refits at each tuning iteration.

source(here::here("R", "00_config.R"))

make_outer_resampling <- function(task) {
  rs <- mlr3::rsmp("cv", folds = RESAMPLING$outer_folds)
  if (RESAMPLING$repeats > 1) {
    rs <- mlr3::rsmp("repeated_cv",
                     folds = RESAMPLING$outer_folds,
                     repeats = RESAMPLING$repeats)
  }
  rs$instantiate(task)
  rs
}

make_inner_resampling <- function() {
  mlr3::rsmp("cv", folds = RESAMPLING$inner_folds)
}

# Wrap a learner with an auto-tuner using inner resampling.
auto_tune <- function(learner, task_type = c("surv", "classif")) {
  task_type <- match.arg(task_type)
  measure <- if (task_type == "surv") {
    mlr3::msr("surv.cindex")
  } else {
    mlr3::msr("classif.auc")
  }
  mlr3tuning::auto_tuner(
    tuner       = mlr3tuning::tnr(TUNING$method),
    learner     = learner,
    resampling  = make_inner_resampling(),
    measure     = measure,
    term_evals  = TUNING$n_evals,
    store_models = FALSE
  )
}

# Build the omics-screening pipeline op. Applied as part of the learner graph
# so that the top-K omics features are re-selected on the training fold only.
omics_filter_po <- function(top_k = 50, screen = "cox_univariate") {
  # mlr3filters offers cox.univariate and auc filters; we use cox for surv,
  # auc for binary classification.
  filter_id <- if (screen == "cox_univariate") "surv.univariate"
               else if (screen == "logistic_univariate") "auc"
               else stop("Unsupported screen: ", screen)
  mlr3pipelines::po("filter",
                    filter = mlr3filters::flt(filter_id),
                    filter.nfeat = top_k)
}
