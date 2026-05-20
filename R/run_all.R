# run_all.R
# Convenience: run the pipeline end-to-end. Each step also runs standalone
# when sourced directly.

source(here::here("R", "00_config.R"))

steps <- c(
  "R/01_load_clinical.R",
  "R/02_load_omics.R",
  "R/03_qc_and_merge.R",
  "R/04_missingness_mice.R",
  "R/05_feature_engineering.R",
  "R/06_cutpoints.R",
  "R/07_task_definitions.R",
  "R/10_benchmark.R",
  "R/11_final_models.R",
  "R/12_interpretability.R",
  "R/13_calibration_dca.R",
  "R/99_report.R"
)

for (s in steps) {
  log_msg("==> ", s)
  source(here::here(s), local = new.env())
}
