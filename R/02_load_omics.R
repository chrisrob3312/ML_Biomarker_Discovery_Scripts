# 02_load_omics.R
# Loaders for each omics block. Each loader returns a matrix with patient_id
# rownames OR a tibble with a patient_id column. Files are expected to be
# pre-normalized (RNA VST, methylation M-values, metabolomics on log scale).

source(here::here("R", "00_config.R"))

load_omics_block <- function(block_name) {
  block <- CONFIG$omics[[block_name]]
  if (is.null(block)) stop("Unknown omics block: ", block_name)
  path <- here::here(block$file)
  if (!file.exists(path)) {
    if (isTRUE(block$optional)) {
      log_msg("Optional omics block '", block_name, "' not found - skipping.")
      return(NULL)
    }
    log_msg("Omics file missing for '", block_name, "' at ", path,
            " - returning empty matrix.")
    return(matrix(numeric(0), nrow = 0, ncol = 0,
                  dimnames = list(character(), character())))
  }
  obj <- readRDS(path)
  validate_omics_matrix(obj, block_name)
  obj
}

validate_omics_matrix <- function(obj, block_name) {
  if (!is.matrix(obj)) stop(block_name, ": expected a matrix.")
  if (is.null(rownames(obj))) {
    stop(block_name, ": matrix must have patient_id rownames.")
  }
  if (anyNA(obj)) {
    log_msg("Warning: ", sum(is.na(obj)), " NAs in ", block_name,
            " - will be handled by feature engineering / imputation.")
  }
  invisible(TRUE)
}

load_all_omics <- function() {
  blocks <- names(CONFIG$omics)
  setNames(lapply(blocks, load_omics_block), blocks)
}

if (sys.nframe() == 0) {
  omics <- load_all_omics()
  for (b in names(omics)) {
    m <- omics[[b]]
    log_msg(sprintf("  %s: %d patients x %d features",
                    b, NROW(m), NCOL(m)))
  }
}
