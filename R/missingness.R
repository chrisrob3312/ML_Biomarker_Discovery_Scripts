# missingness.R
# Visualize missingness pattern and produce a per-variable summary.
# Used to inform: (a) which variables can be imputed, (b) MAR plausibility,
# (c) which auxiliary variables to include in MICE.

source(here::here("R", "00_config.R"))
suppressPackageStartupMessages({
  library(naniar)
  library(VIM)
})

missingness_table <- function(dat) {
  m <- naniar::miss_var_summary(dat)
  m$n_total <- nrow(dat)
  m
}

# Save a missingness heatmap (variables x patients, sorted by missingness)
# and an aggregation plot (patterns of co-missingness).
missingness_plots <- function(dat,
                              out_dir = file.path(PATHS$results, "step2")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  p_heat <- naniar::vis_miss(dat, warn_large_data = FALSE) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 75, size = 6))
  ggplot2::ggsave(file.path(out_dir, "missingness_heatmap.png"),
                  p_heat, width = 10, height = 6, dpi = 150)

  p_upset <- naniar::gg_miss_upset(dat, nsets = 12)
  grDevices::png(file.path(out_dir, "missingness_upset.png"),
                 width = 1000, height = 600, res = 120)
  print(p_upset); grDevices::dev.off()

  log_msg("Missingness plots in ", out_dir)
  invisible(NULL)
}

# Little's MCAR test (formal) - p > 0.05 is consistent with MCAR.
# Even when this rejects, you can still defend MAR if you have rich
# auxiliary variables in the imputation model.
little_mcar <- function(dat) {
  if (!requireNamespace("naniar", quietly = TRUE)) return(NULL)
  numeric_dat <- dat[, sapply(dat, is.numeric), drop = FALSE]
  res <- try(naniar::mcar_test(numeric_dat), silent = TRUE)
  if (inherits(res, "try-error")) return(NULL)
  res
}
