# 03_qc_and_merge.R
# Merge clinical + risk alleles into a single analytic table.
# Omics blocks are kept separate and joined later inside the resampling loop
# (to allow CV-internal feature screening without leakage).

source(here::here("R", "01_load_clinical.R"))
source(here::here("R", "02_load_omics.R"))

build_analytic_table <- function() {
  clin <- load_clinical()
  ra   <- load_risk_alleles()

  dat <- if (nrow(ra) > 0) {
    dplyr::left_join(clin, ra, by = CONFIG$patient_id)
  } else clin

  dat <- apply_cohort_filters(dat)
  qc_report(dat)
  dat
}

apply_cohort_filters <- function(dat) {
  cf <- CONFIG$cohort_filters
  if (is.null(cf)) return(dat)
  n0 <- nrow(dat)
  if (!is.null(cf$min_dx_year)) {
    dat <- dat[!is.na(dat$dx_year) & dat$dx_year >= cf$min_dx_year, , drop = FALSE]
  }
  if (!is.null(cf$max_dx_year)) {
    dat <- dat[!is.na(dat$dx_year) & dat$dx_year <= cf$max_dx_year, , drop = FALSE]
  }
  log_msg(sprintf("Cohort filter dx_year in [%s, %s]: %d -> %d (%d dropped)",
                  cf$min_dx_year %||% "-Inf", cf$max_dx_year %||% "Inf",
                  n0, nrow(dat), n0 - nrow(dat)))
  dat
}

qc_report <- function(dat) {
  log_msg("QC summary -------------------------------------------------")
  log_msg("  N patients: ", nrow(dat))
  if (nrow(dat) == 0) return(invisible(NULL))

  # Outcome event counts
  for (oc in names(CONFIG$outcomes)) {
    o <- CONFIG$outcomes[[oc]]
    ev <- dat[[o$event_col]]
    if (!is.null(ev)) {
      log_msg(sprintf("  %s events: %d / %d (%.1f%%)",
                      oc, sum(ev == 1, na.rm = TRUE), sum(!is.na(ev)),
                      100 * mean(ev == 1, na.rm = TRUE)))
    }
  }

  # Missingness top-10
  miss <- sort(colMeans(is.na(dat)), decreasing = TRUE)
  miss <- miss[miss > 0][1:min(10, length(miss))]
  if (length(miss) > 0) {
    log_msg("  Top missingness:")
    for (i in seq_along(miss)) {
      log_msg(sprintf("    %s: %.1f%%", names(miss)[i], 100 * miss[i]))
    }
  }
  invisible(NULL)
}

save_analytic_table <- function(dat,
                                path = file.path(PATHS$data_int, "analytic.rds")) {
  saveRDS(dat, path)
  log_msg("Saved analytic table to ", path)
}

if (sys.nframe() == 0) {
  dat <- build_analytic_table()
  if (nrow(dat) > 0) save_analytic_table(dat)
}
