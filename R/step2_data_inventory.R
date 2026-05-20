# step2_data_inventory.R
# Step 2 orchestrator. Runs end-to-end the data-inventory phase:
#   1. Load + filter cohort, write CONSORT
#   2. Median follow-up via reverse Kaplan-Meier
#   3. Table 1 (overall, by DFS event, by ethnicity, by era)
#   4. Missingness summary + heatmap + UpSet
#   5. Sample-size adequacy via pmsampsize
#   6. Decision support for dx-window extension
#
# Outputs into results/step2/.

source(here::here("R", "00_config.R"))
source(here::here("R", "01_load_clinical.R"))
source(here::here("R", "02_load_omics.R"))
source(here::here("R", "03_qc_and_merge.R"))
source(here::here("R", "05_feature_engineering.R"))
source(here::here("R", "07_task_definitions.R"))
source(here::here("R", "tableone.R"))
source(here::here("R", "missingness.R"))
source(here::here("R", "sample_size.R"))

OUT <- file.path(PATHS$results, "step2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# 1. Load + filter
log_msg("Step 2.1 - cohort assembly")
dat <- build_analytic_table()
dat <- engineer_features(dat)
saveRDS(dat, file.path(PATHS$data_int, "engineered.rds"))

consort <- attr(dat, "consort")
readr::write_csv(consort, file.path(OUT, "consort.csv"))

# 2. Median follow-up (reverse KM)
log_msg("Step 2.2 - median follow-up")
reverse_km <- function(time, event) {
  fit <- survival::survfit(survival::Surv(time, 1 - event) ~ 1)
  q <- stats::quantile(fit, probs = 0.5)
  list(median = q$quantile, lci = q$lower, uci = q$upper)
}
fup_summary <- list(
  dfs = reverse_km(dat$dfs_time_months,    dat$dfs_event),
  os  = reverse_km(dat$os_eoi_time_months, dat$os_eoi_event)
)
saveRDS(fup_summary, file.path(OUT, "followup_summary.rds"))

# 3. Table 1
log_msg("Step 2.3 - Table 1")
save_table1(table1_by_event(dat, "dfs_event"),
            file.path(OUT, "table1_by_dfs.html"))
save_table1(table1_by_event(dat, "os_eoi_event"),
            file.path(OUT, "table1_by_os.html"))
save_table1(table1_by_ethnicity(dat),
            file.path(OUT, "table1_by_ethnicity.html"))
save_table1(table1_by_era(dat),
            file.path(OUT, "table1_by_era.html"))

# 4. Missingness
log_msg("Step 2.4 - missingness")
miss_tab <- missingness_table(dat)
readr::write_csv(miss_tab, file.path(OUT, "missingness_summary.csv"))
missingness_plots(dat, OUT)
mcar <- little_mcar(dat)
if (!is.null(mcar)) saveRDS(mcar, file.path(OUT, "little_mcar.rds"))

# 5. Sample-size adequacy
log_msg("Step 2.5 - sample-size calc")
feats <- core_features(dat)
ss <- run_all_sample_size(dat, feats, horizon_months = 36)
saveRDS(ss, file.path(OUT, "sample_size.rds"))
sink(file.path(OUT, "sample_size.txt"))
print(ss)
sink()

# 6. Decision: dx-window extension
log_msg("Step 2.6 - extension decision")
n_total      <- nrow(dat)
n_dfs_events <- sum(dat$dfs_event == 1, na.rm = TRUE)
n_os_events  <- sum(dat$os_eoi_event == 1, na.rm = TRUE)

tier <- dplyr::case_when(
  n_total >= 1500 & n_dfs_events >= 200 & n_os_events >= 100 ~ "Tier 1: keep 2010 cutoff",
  n_total >= 1000 & n_dfs_events >= 100 & n_os_events >= 50  ~ "Tier 2: keep 2010, lean on penalization",
  n_total >= 700  & n_dfs_events >= 70  & n_os_events >= 35  ~ "Tier 3: extend to 2005",
  TRUE                                                         ~ "Tier 4: extend to 2000 + prune features / PRS"
)
decision <- list(n_total = n_total,
                 n_dfs_events = n_dfs_events,
                 n_os_events = n_os_events,
                 tier = tier)
saveRDS(decision, file.path(OUT, "extension_decision.rds"))
log_msg("Decision tier: ", tier)
log_msg("Step 2 complete -> ", OUT)
