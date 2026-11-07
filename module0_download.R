# ============================================================================
# module0_download.R — Download cMD3 paired training cohorts (BritoIL_2016,
# HMP_2012) via ExperimentHub.
# Input:  none (downloads from ExperimentHub)
# Output:
#   <data_dir>/brito_raw.rds        = list(mat, meta, study)
#   <data_dir>/hmp_raw.rds          = list(mat, meta, study)
#   <data_dir>/sample_metadata.rds  = full cMD3 sample metadata
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
  library(dplyr)
})

## ---- Load one study: relative-abundance species matrix + metadata ----
load_study <- function(study) {
  cat("\n==========", study, "==========\n")
  res <- curatedMetagenomicData(paste0(study, ".*relative_abundance"),
                                dryrun = FALSE)
  es  <- res[[1]]
  m   <- as.matrix(assay(es))
  pd  <- as.data.frame(colData(es))
  
  # colData lacks a sample_id column; restore it from rownames
  pd$sample_id <- rownames(pd)
  stopifnot(identical(pd$sample_id, colnames(m)))
  
  # Harmonize gender -> sex
  pd$sex <- if ("gender" %in% colnames(pd)) pd$gender else NA_character_
  
  cat("dimensions:", paste(dim(m), collapse = " x "), "\n")
  cat("body_site:", paste(unique(pd$body_site), collapse = ","), "\n")
  cat("colSums range:", range(colSums(m, na.rm = TRUE)), "\n")
  
  list(mat = m, meta = pd, study = study)
}

## ---- Download ----
brito_raw <- load_study("BritoIL_2016")
hmp_raw   <- load_study("HMP_2012")

saveRDS(brito_raw, file.path(data_dir, "brito_raw.rds"))
saveRDS(hmp_raw,   file.path(data_dir, "hmp_raw.rds"))

## ---- Full cMD3 sample metadata (reference) ----
sm <- curatedMetagenomicData::sampleMetadata
saveRDS(sm, file.path(data_dir, "sample_metadata.rds"))

cat("\nmodule0_download DONE\n")
