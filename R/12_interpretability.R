# 12_interpretability.R
# - SHAP for tree models (binary + survival risk score)
# - Partial dependence / ALE
# - Forest plot of adjusted HRs (baseline Cox + RCS)
# - Nomogram from rms::cph
# - KM curves by predicted risk tertile

source(here::here("R", "00_config.R"))

shap_xgb <- function(learner, task, n_bg = 200) {
  X <- task$data(cols = task$feature_names)
  bg <- X[sample(nrow(X), min(n_bg, nrow(X))), ]
  pred_wrap <- function(object, newdata) {
    learner$predict_newdata(newdata)$crank
  }
  fastshap::explain(learner, X = bg, pred_wrapper = pred_wrap, nsim = 50)
}

forest_plot_cox <- function(cph_fit, out_path) {
  s <- summary(cph_fit)
  tab <- data.frame(
    var = rownames(s$coefficients),
    hr  = exp(s$coefficients[, "coef"]),
    lo  = exp(s$coefficients[, "coef"] - 1.96 * s$coefficients[, "se(coef)"]),
    hi  = exp(s$coefficients[, "coef"] + 1.96 * s$coefficients[, "se(coef)"]),
    p   = s$coefficients[, "Pr(>|z|)"]
  )
  p <- ggplot2::ggplot(tab,
                       ggplot2::aes(x = hr, y = forcats::fct_reorder(var, hr))) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0.2) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "HR (95% CI)", y = NULL) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(out_path, p, width = 7, height = max(3, 0.25 * nrow(tab)))
  tab
}

nomogram_from_cph <- function(cph_fit, out_path) {
  dd <- rms::datadist(cph_fit$x); options(datadist = "dd")
  nm <- rms::nomogram(cph_fit, fun = list(
    "1y survival" = function(x) plogis(x),
    "3y survival" = function(x) plogis(x),
    "5y survival" = function(x) plogis(x)
  ))
  grDevices::pdf(out_path, width = 10, height = 8)
  plot(nm)
  grDevices::dev.off()
  invisible(nm)
}

km_by_risk_group <- function(pred_risk, time, event, n_groups = 3, out_path) {
  grp <- cut(pred_risk,
             breaks = stats::quantile(pred_risk, seq(0, 1, length.out = n_groups + 1),
                                      na.rm = TRUE),
             include.lowest = TRUE, labels = paste0("Q", seq_len(n_groups)))
  d <- data.frame(time = time, event = event, grp = grp)
  fit <- survival::survfit(survival::Surv(time, event) ~ grp, data = d)
  p <- survminer::ggsurvplot(fit, data = d, risk.table = TRUE,
                             pval = TRUE, conf.int = TRUE)
  ggplot2::ggsave(out_path, p$plot, width = 7, height = 5)
  invisible(fit)
}
