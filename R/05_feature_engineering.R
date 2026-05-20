# 05_feature_engineering.R
# Deterministic transforms (no learning from outcome). These are safe to apply
# before splitting. Outcome-aware steps (omics screening, maxstat cuts) live
# in 06_cutpoints.R and the resampling pipeline.

source(here::here("R", "00_config.R"))
source(here::here("R", "era_lookup.R"))

engineer_features <- function(dat) {
  dat <- dat |>
    dplyr::mutate(
      wbc_dx_log10 = ifelse(wbc_dx > 0, log10(wbc_dx), NA_real_),
      treatment_era = ifelse(is.na(treatment_era) | treatment_era == "",
                             as.character(derive_treatment_era(
                               dx_year,
                               min_year = CONFIG$cohort_filters$min_dx_year %||% 2010,
                               max_year = CONFIG$cohort_filters$max_dx_year %||% 2025)),
                             treatment_era),
      # Nested factor: protocol_arm nested within treatment_era. When
      # protocol_arm is unknown we record "<era>/unknown" so the era effect
      # is still identified for that patient. Baseline Cox uses era and
      # era:protocol_arm separately (see 11_final_models.R); ML learners
      # use this single combined factor.
      protocol_arm_clean = dplyr::if_else(
        is.na(protocol_arm) | protocol_arm == "" | tolower(protocol_arm) == "unknown",
        "unknown", protocol_arm),
      protocol_nested = factor(paste(treatment_era, protocol_arm_clean, sep = "/")),

      # MRD: floor at limit of detection then log10. Adjust LoD to your assay.
      mrd_eoi_log10 = log10(pmax(mrd_eoi_continuous, 1e-5)),
      mrd_eoc_log10 = log10(pmax(mrd_eoc_continuous, 1e-5)),

      # Pre-specified clinical cuts (kept as parallel categorical variables;
      # use these in stratified tables and the clinical reference model).
      age_cat_clin = cut(
        age_at_dx_years,
        breaks = c(-Inf, CONFIG$clinical_cutpoints$age_at_dx_years, Inf),
        labels = c("<1", "1-<10", "10-<16", ">=16"),
        right  = FALSE),
      wbc_cat_clin = factor(
        ifelse(wbc_dx >= CONFIG$clinical_cutpoints$wbc_dx[1], "ge50", "lt50"),
        levels = c("lt50", "ge50")),
      mrd_eoi_cat_clin = cut(
        mrd_eoi_continuous * 100,  # convert to percent
        breaks = c(-Inf, CONFIG$clinical_cutpoints$mrd_eoi_percent, Inf),
        labels = c("neg", "low_pos", "mid_pos", "high_pos"),
        right  = FALSE)
    ) |>
    dplyr::mutate(dplyr::across(
      c(sex, race_ethnicity, cns_status, cyto_subtype, ikzf1_status,
        treatment_type, nci_risk, insurance_payer),
      \(x) forcats::fct_explicit_na(as.factor(x), na_level = "unknown")
    ))

  dat
}

# Restricted-cubic-spline expansion for baseline Cox / GLM. Uses rms::rcs,
# 4 knots at Harrell's default quantiles. Returned as a model matrix piece
# that you cbind to the main matrix.
rcs_expand <- function(x, knots = 4) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  out <- matrix(NA_real_, nrow = length(x), ncol = knots - 1)
  if (sum(ok) >= knots + 2) {
    out[ok, ] <- rms::rcs(x[ok], knots)[, -1, drop = FALSE]
  }
  colnames(out) <- paste0("rcs", seq_len(ncol(out)))
  out
}

# Outcome-aware top-K omics screen. Call ONLY on training-fold data.
# Returns selected feature names. Use Cox univariate score for survival,
# univariate logistic z for binary relapse.
screen_omics <- function(x, y, top_k, screen = c("cox_univariate",
                                                 "logistic_univariate",
                                                 "limma")) {
  screen <- match.arg(screen)
  if (NCOL(x) == 0) return(character())

  scores <- switch(screen,
    cox_univariate = {
      stopifnot(inherits(y, "Surv"))
      apply(x, 2, function(z) {
        z <- as.numeric(z); ok <- is.finite(z)
        if (sum(ok) < 10 || sd(z[ok]) == 0) return(NA_real_)
        fit <- try(survival::coxph(y[ok] ~ z[ok]), silent = TRUE)
        if (inherits(fit, "try-error")) NA_real_
        else abs(summary(fit)$coefficients[1, "z"])
      })
    },
    logistic_univariate = {
      apply(x, 2, function(z) {
        z <- as.numeric(z); ok <- is.finite(z) & !is.na(y)
        if (sum(ok) < 10 || sd(z[ok]) == 0) return(NA_real_)
        fit <- try(glm(y[ok] ~ z[ok], family = binomial()), silent = TRUE)
        if (inherits(fit, "try-error")) NA_real_
        else abs(summary(fit)$coefficients[2, "z value"])
      })
    },
    limma = stop("limma screen not yet implemented; "
                 , "wire in if you want moderated stats.")
  )

  ord <- order(scores, decreasing = TRUE, na.last = NA)
  head(colnames(x)[ord], top_k)
}

# Ancestry-dosed risk allele features. For each risk allele, produces:
#   <id>_dosage                  (total dosage 0/1/2)
#   <id>_la_<ancestry>_dose      (local-ancestry-attributed dosage 0/1/2)
# These columns are expected to already exist in the analytic table -
# this helper validates and z-standardizes them (within fold).
ancestry_dosed_columns <- function(dat) {
  cols <- character()
  for (ra in CONFIG$risk_alleles) {
    cols <- c(cols, paste0(ra$id, "_dosage"))
    for (anc in ra$ancestries) {
      cols <- c(cols, sprintf("%s_la_%s_dose", ra$id, anc))
    }
  }
  intersect(cols, names(dat))
}

if (sys.nframe() == 0) {
  dat <- readRDS(file.path(PATHS$data_int, "analytic.rds"))
  dat <- engineer_features(dat)
  saveRDS(dat, file.path(PATHS$data_int, "engineered.rds"))
  log_msg("Engineered table saved.")
}
