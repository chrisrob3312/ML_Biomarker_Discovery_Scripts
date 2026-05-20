# tableone.R
# Table 1 builders using gtsummary. Three views:
#   (a) Overall + stratified by event status
#   (b) Overall + stratified by Hispanic/non-Hispanic ethnicity
#   (c) Overall + stratified by treatment era
#
# These are descriptive only. They do NOT inform variable selection
# and should be inspected BEFORE outcome-aware analyses begin.

source(here::here("R", "00_config.R"))
suppressPackageStartupMessages({
  library(gtsummary)
  library(gt)
})

# Variables to summarize in Table 1. Edit here if you add/remove core vars.
table1_vars <- function(dat) {
  v <- c(
    "age_at_dx_years", "sex", "race_ethnicity",
    "global_ancestry_eur", "global_ancestry_afr",
    "global_ancestry_amr", "global_ancestry_eas",
    "wbc_dx", "cns_status", "bmi_z", "nci_risk",
    "cyto_subtype", "ph_positive", "ph_like", "kmt2a_rearr",
    "etv6_runx1", "tcf3_pbx1", "hyperdiploid_high", "hypodiploid_low",
    "iamp21", "crlf2_rearr", "ikzf1_plus_flag",
    "mrd_eoi_continuous", "mrd_eoi_category",
    "adi_national_decile", "insurance_payer",
    "treatment_era", "clinical_trial_enrolled"
  )
  intersect(v, names(dat))
}

# Helper: standardize Hispanic flag from race_ethnicity if not stored separately.
add_hispanic_flag <- function(dat) {
  if ("hispanic" %in% names(dat)) return(dat)
  dat$hispanic <- ifelse(grepl("hispanic|latin", dat$race_ethnicity,
                               ignore.case = TRUE),
                         "Hispanic/Latino", "Non-Hispanic")
  dat
}

table1_by_event <- function(dat, event_col = "dfs_event") {
  d <- dat[!is.na(dat[[event_col]]), , drop = FALSE]
  d$.event <- factor(d[[event_col]], levels = c(0, 1),
                     labels = c("No event", "Event"))
  gtsummary::tbl_summary(
    d[, c(table1_vars(d), ".event"), drop = FALSE],
    by = ".event",
    missing = "ifany",
    statistic = list(
      gtsummary::all_continuous()  ~ "{median} ({p25}, {p75})",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    )
  ) |>
    gtsummary::add_overall() |>
    gtsummary::add_p() |>
    gtsummary::modify_caption(
      sprintf("Table 1 (stratified by %s).", event_col))
}

table1_by_ethnicity <- function(dat) {
  d <- add_hispanic_flag(dat)
  gtsummary::tbl_summary(
    d[, c(setdiff(table1_vars(d), "race_ethnicity"), "hispanic"),
      drop = FALSE],
    by = "hispanic",
    missing = "ifany"
  ) |>
    gtsummary::add_overall() |>
    gtsummary::add_p() |>
    gtsummary::modify_caption("Table 1B (by Hispanic/Latino ethnicity).")
}

table1_by_era <- function(dat) {
  gtsummary::tbl_summary(
    dat[, c(setdiff(table1_vars(dat), "treatment_era"), "treatment_era"),
        drop = FALSE],
    by = "treatment_era",
    missing = "ifany"
  ) |>
    gtsummary::add_overall() |>
    gtsummary::modify_caption("Table 1C (by treatment era).")
}

save_table1 <- function(tbl, path) {
  gt::gtsave(gtsummary::as_gt(tbl), filename = path)
  log_msg("Wrote ", path)
}
