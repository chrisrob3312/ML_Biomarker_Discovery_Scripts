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

  consort <- data.frame(step = character(), n = integer(),
                        stringsAsFactors = FALSE)
  add_step <- function(step, n) consort <<- rbind(consort, data.frame(step, n))
  add_step("loaded", nrow(dat))

  if (!is.null(cf$min_dx_year)) {
    dat <- dat[!is.na(dat$dx_year) & dat$dx_year >= cf$min_dx_year, , drop = FALSE]
    add_step(sprintf("dx_year >= %d", cf$min_dx_year), nrow(dat))
  }
  if (!is.null(cf$max_dx_year)) {
    dat <- dat[!is.na(dat$dx_year) & dat$dx_year <= cf$max_dx_year, , drop = FALSE]
    add_step(sprintf("dx_year <= %d", cf$max_dx_year), nrow(dat))
  }
  if (!is.null(cf$min_age_years)) {
    dat <- dat[!is.na(dat$age_at_dx_years) &
                 dat$age_at_dx_years >= cf$min_age_years, , drop = FALSE]
    add_step(sprintf("age >= %g", cf$min_age_years), nrow(dat))
  }
  if (!is.null(cf$max_age_years)) {
    dat <- dat[!is.na(dat$age_at_dx_years) &
                 dat$age_at_dx_years <= cf$max_age_years, , drop = FALSE]
    add_step(sprintf("age <= %g", cf$max_age_years), nrow(dat))
  }
  if (!is.null(cf$lineage_keep)) {
    dat <- dat[!is.na(dat$lineage) & dat$lineage %in% cf$lineage_keep, , drop = FALSE]
    add_step(sprintf("lineage in {%s}",
                     paste(cf$lineage_keep, collapse = ",")), nrow(dat))
  }
  if (isTRUE(cf$exclude_down_syndrome)) {
    dat <- dat[is.na(dat$down_syndrome) | dat$down_syndrome == 0, , drop = FALSE]
    add_step("down_syndrome == 0", nrow(dat))
  }
  if (isTRUE(cf$exclude_mpal)) {
    dat <- dat[is.na(dat$lineage) | dat$lineage != "MPAL", , drop = FALSE]
    add_step("not MPAL", nrow(dat))
  }
  if (isTRUE(cf$exclude_relapsed_at_dx)) {
    dat <- dat[is.na(dat$relapsed_at_dx) | dat$relapsed_at_dx == 0, , drop = FALSE]
    add_step("de novo (not relapsed_at_dx)", nrow(dat))
  }
  if (isTRUE(cf$require_eoi)) {
    dat <- dat[!is.na(dat$reached_eoi) & dat$reached_eoi == 1, , drop = FALSE]
    add_step("reached_eoi == 1", nrow(dat))
  }

  attr(dat, "consort") <- consort
  log_msg("CONSORT --------------------------------------")
  for (i in seq_len(nrow(consort))) {
    log_msg(sprintf("  %-30s n = %d", consort$step[i], consort$n[i]))
  }
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
