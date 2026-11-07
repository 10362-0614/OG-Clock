# ============================================================================
# module4_prep_gut.R — Prepare RDS files for the three external gut
# validation cohorts
#
# Output (written to <data_dir>, consumed by module4_gut_validation.R):
#   <data_dir>/ext_LifeLinesDeep_2016.rds
#   <data_dir>/ext_AsnicarF_2021.rds
#   <data_dir>/ext_KarlssonFH_2013.rds
#
# Each file: list(mat = species x samples relative-abundance matrix,
#                 meta = data.frame with age / subject_id / body_site, ...)
#
# Requires: Bioconductor package curatedMetagenomicData (installed by
# module0_setup.R)
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(dplyr)
})

studies <- c("LifeLinesDeep_2016", "AsnicarF_2021", "KarlssonFH_2013")

for (st in studies) {
  fp <- file.path(data_dir, paste0("ext_", st, ".rds"))
  if (file.exists(fp)) {
    cat("Already exists, skipping:", fp, "
")
    next
  }
  
  cat("
====", st, "====
")
  # relative_abundance: species x samples (MetaPhlAn3/4 lineage row names)
  res <- curatedMetagenomicData(paste0(st, ".*relative_abundance"),
                                dryrun = FALSE)
  es  <- res[[1]]
  
  m  <- as.matrix(SummarizedExperiment::assay(es))
  pd <- as.data.frame(SummarizedExperiment::colData(es))
  
  # Keep samples with age only
  keep <- !is.na(pd$age) & pd$age > 0
  m  <- m[, keep, drop = FALSE]
  pd <- pd[keep, , drop = FALSE]
  
  cat("Samples with age:", ncol(m), " | features:", nrow(m), "
")
  cat("Age range:", paste(range(pd$age), collapse = "-"), "
")
  
  stopifnot(ncol(m) > 100)          # sanity check
  stopifnot(identical(colnames(m), rownames(pd)))
  
  saveRDS(list(mat = m, meta = pd), fp)
  cat("Saved:", fp, "
")
}

cat("
========== module4_prep_gut DONE ==========
")
cat("Three external gut cohorts ready; run module4_gut_validation.R next
")
