# 01_load_clinical.R
# Loads the clinical/demographic/outcome table. Replace the placeholder
# loader with your real source (REDCap pull, csv, database) when ready.

source(here::here("R", "00_config.R"))

load_clinical <- function(path = file.path(PATHS$data_raw, "clinical.csv")) {
  if (!file.exists(path)) {
    log_msg("Clinical file not found at ", path, " - returning empty template.")
    return(clinical_template())
  }
  dat <- readr::read_csv(path, show_col_types = FALSE)
  validate_clinical(dat)
  dat
}

# Expected schema; used both as a template when no data is present and
# as the contract that the QC step enforces.
clinical_template <- function() {
  tibble::tibble(
    patient_id              = character(),
    age_at_dx_years         = numeric(),
    sex                     = character(),         # M / F
    race_ethnicity          = character(),         # NHW, NHB, Hispanic, Asian, Other
    wbc_dx                  = numeric(),           # x10^9/L
    cns_status              = character(),         # CNS1 / CNS2 / CNS3
    bmi_z                   = numeric(),
    down_syndrome           = integer(),           # 0/1
    clinical_trial_enrolled = integer(),           # 0/1
    treatment_type          = character(),
    protocol_arm            = character(),         # e.g. AALL0232_HD-MTX_dex
    treatment_era           = character(),         # derived if arm missing
    nci_risk                = character(),         # SR / HR
    dx_year                 = integer(),
    cyto_subtype            = character(),
    ikzf1_status            = character(),         # wt / del / ikzf1plus
    crlf2_rearr             = integer(),           # 0/1
    ph_positive             = integer(),
    ph_like                 = integer(),
    kmt2a_rearr             = integer(),
    etv6_runx1              = integer(),
    tcf3_pbx1               = integer(),
    hyperdiploid_high       = integer(),
    hypodiploid_low         = integer(),
    iamp21                  = integer(),
    dux4_rearr              = integer(),
    mef2d_rearr             = integer(),
    zfn384_rearr            = integer(),
    paxbike_alt             = integer(),
    tp53_alt                = integer(),
    ikzf1_plus_flag         = integer(),
    cdkn2ab_del             = integer(),
    jak_pathway_alt         = integer(),
    ras_pathway_alt         = integer(),
    mrd_eoi_continuous      = numeric(),           # fraction (e.g. 0.0001)
    mrd_eoi_category        = character(),         # neg / low_pos / high_pos
    mrd_eoc_continuous      = numeric(),
    day8_pb_blasts          = numeric(),
    day15_marrow_blasts     = numeric(),
    global_ancestry_eur     = numeric(),
    global_ancestry_afr     = numeric(),
    global_ancestry_eas     = numeric(),
    global_ancestry_amr     = numeric(),
    global_ancestry_sas     = numeric(),
    adi_national_decile     = integer(),
    adi_state_decile        = integer(),
    insurance_payer         = character(),
    # End-of-induction landmark
    reached_eoi             = integer(),           # 0/1: completed induction
    induction_outcome       = character(),         # CR / failure / death / other
    time_dx_to_eoi_months   = numeric(),           # used as offset / for landmark
    # EOI-anchored outcomes
    dfs_time_months         = numeric(),           # from EOI to event/censor
    dfs_event               = integer(),           # 1 = relapse|SMN|death
    dfs_cause               = character(),         # relapse / smn / death / censored
    os_eoi_time_months      = numeric(),
    os_eoi_event            = integer(),
    # Ever-relapsed (for binary relapse model)
    relapse_event           = integer()
  )
}

validate_clinical <- function(dat) {
  required <- names(clinical_template())
  missing  <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop("Clinical table missing required columns: ",
         paste(missing, collapse = ", "))
  }
  if (anyDuplicated(dat$patient_id)) {
    stop("Duplicate patient_id values in clinical table.")
  }
  invisible(TRUE)
}

# Risk allele + local ancestry dosing.
# Expected columns per allele in CONFIG$risk_alleles:
#   <id>_dosage                       (0/1/2)
#   <id>_la_<ancestry>_dose           (0/1/2) for each ancestry
load_risk_alleles <- function(path = file.path(PATHS$data_raw, "risk_alleles.csv")) {
  if (!file.exists(path)) {
    log_msg("Risk allele file not found - returning empty tibble.")
    return(tibble::tibble(patient_id = character()))
  }
  readr::read_csv(path, show_col_types = FALSE)
}

if (sys.nframe() == 0) {
  clin <- load_clinical()
  ra   <- load_risk_alleles()
  log_msg("Clinical rows: ", nrow(clin), " | risk-allele rows: ", nrow(ra))
}
