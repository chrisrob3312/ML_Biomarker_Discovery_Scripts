# 08_learners.R
# Benchmark roster. Aim is a defensible "best-in-class" comparison across
# linear, penalized linear, tree-ensemble, gradient-boosted, and stacked
# learners, mirroring the standard clinical-ML benchmark set.
#
# Survival benchmark (mlr3proba):
#   - surv.kaplan          : Kaplan-Meier null reference
#   - surv.coxph           : Cox PH baseline (with RCS in the recipe)
#   - surv.glmnet          : Elastic-net Cox (linear, sparse)
#   - surv.ranger          : Random Survival Forest
#   - surv.rfsrc           : randomForestSRC RSF (alt implementation)
#   - surv.xgboost.cox     : Gradient-boosted Cox
#   - surv.gbm             : Generalized boosted Cox (alt)
#   - surv.aorsf           : Accelerated Oblique Random Survival Forest
#                            (state-of-the-art; Jaeger et al. 2024)
#   - surv.deephit / surv.deepsurv : deep nets (optional; needs torch)
#   - super learner        : weighted stack of the above
#
# Binary relapse benchmark (mlr3):
#   - classif.featureless  : null
#   - classif.log_reg      : logistic baseline
#   - classif.glmnet       : elastic-net logistic
#   - classif.ranger       : Random Forest
#   - classif.xgboost      : XGBoost
#   - classif.lightgbm     : LightGBM (optional)
#   - classif.svm          : RBF SVM (sanity check)
#   - super learner        : stack

source(here::here("R", "00_config.R"))

# ---- Survival learners --------------------------------------------------
make_surv_learners <- function() {
  lrns <- list(
    null   = mlr3::lrn("surv.kaplan"),
    coxph  = mlr3::lrn("surv.coxph"),
    glmnet = mlr3::lrn("surv.glmnet",
                      alpha = paradox::to_tune(0, 1),
                      s     = paradox::to_tune(1e-4, 1, logscale = TRUE)),
    ranger = mlr3::lrn("surv.ranger",
                      num.trees    = 1000,
                      mtry.ratio   = paradox::to_tune(0.1, 0.9),
                      min.node.size= paradox::to_tune(3, 30)),
    rfsrc  = mlr3::lrn("surv.rfsrc",
                      ntree    = 1000,
                      mtry     = paradox::to_tune(2, 30),
                      nodesize = paradox::to_tune(3, 30)),
    xgb_cox= mlr3::lrn("surv.xgboost.cox",
                      nrounds          = paradox::to_tune(100, 1500),
                      eta              = paradox::to_tune(1e-3, 0.3, logscale = TRUE),
                      max_depth        = paradox::to_tune(2, 8),
                      min_child_weight = paradox::to_tune(1, 10),
                      subsample        = paradox::to_tune(0.5, 1),
                      colsample_bytree = paradox::to_tune(0.5, 1))
  )

  # aorsf is a recent (and often top-performing) addition; load if installed.
  if (requireNamespace("aorsf", quietly = TRUE)) {
    lrns$aorsf <- mlr3::lrn("surv.aorsf",
                            n_tree    = 500,
                            mtry      = paradox::to_tune(2, 30),
                            leaf_min_obs = paradox::to_tune(3, 30))
  }
  lrns
}

# ---- Binary classifier learners ----------------------------------------
make_classif_learners <- function() {
  lrns <- list(
    null    = mlr3::lrn("classif.featureless", predict_type = "prob"),
    logreg  = mlr3::lrn("classif.log_reg",     predict_type = "prob"),
    glmnet  = mlr3::lrn("classif.glmnet",      predict_type = "prob",
                        alpha = paradox::to_tune(0, 1),
                        s     = paradox::to_tune(1e-4, 1, logscale = TRUE)),
    ranger  = mlr3::lrn("classif.ranger",      predict_type = "prob",
                        num.trees     = 1000,
                        mtry.ratio    = paradox::to_tune(0.1, 0.9),
                        min.node.size = paradox::to_tune(3, 30)),
    xgb     = mlr3::lrn("classif.xgboost",     predict_type = "prob",
                        nrounds          = paradox::to_tune(100, 1500),
                        eta              = paradox::to_tune(1e-3, 0.3, logscale = TRUE),
                        max_depth        = paradox::to_tune(2, 8),
                        subsample        = paradox::to_tune(0.5, 1),
                        colsample_bytree = paradox::to_tune(0.5, 1)),
    svm     = mlr3::lrn("classif.svm",         predict_type = "prob",
                        kernel = "radial",
                        cost   = paradox::to_tune(1e-2, 1e2, logscale = TRUE),
                        gamma  = paradox::to_tune(1e-3, 1, logscale = TRUE))
  )
  lrns
}

# ---- Super Learner (stack) ----------------------------------------------
# Two options:
#   (a) Native mlr3pipelines stacking with cross-validated base predictions
#       and a glmnet meta-learner (recommended; integrates with the
#       benchmark engine and respects the outer/inner CV).
#   (b) SuperLearner package directly (familiar API; standalone). Use
#       only if you want to mirror a prior SuperLearner workflow exactly.

make_super_learner_surv <- function() {
  base <- make_surv_learners()
  base$null <- NULL  # KM is uninformative as a stack input
  # Each base learner's CV predictions feed a glmnet Cox meta-learner.
  graph <- mlr3pipelines::gunion(lapply(names(base), function(nm) {
    mlr3pipelines::po("learner_cv", learner = base[[nm]],
                      id = paste0("cv_", nm), resampling.folds = 5)
  })) %>>%
    mlr3pipelines::po("featureunion") %>>%
    mlr3::lrn("surv.glmnet", id = "meta_glmnet", alpha = 0.5)
  mlr3::as_learner(graph)
}

make_super_learner_classif <- function() {
  base <- make_classif_learners()
  base$null <- NULL
  graph <- mlr3pipelines::gunion(lapply(names(base), function(nm) {
    mlr3pipelines::po("learner_cv", learner = base[[nm]],
                      id = paste0("cv_", nm), resampling.folds = 5)
  })) %>>%
    mlr3pipelines::po("featureunion") %>>%
    mlr3::lrn("classif.glmnet", id = "meta_glmnet",
              predict_type = "prob", alpha = 0.5)
  mlr3::as_learner(graph)
}
