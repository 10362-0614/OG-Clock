# ============================================================================
# module0_setup.R — Environment setup
# ============================================================================

## ---- 0. Paths (must run first) ----
root_dir <- Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace")
lib_dir  <- file.path(root_dir, ".Rlib")
data_dir <- file.path(root_dir, "data")
res_dir  <- file.path(root_dir, "results")

for (d in c(lib_dir, data_dir, res_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  stopifnot(dir.exists(d))
}
.libPaths(c(lib_dir, .libPaths()))

## ---- 1. Reproducibility ----
set.seed(20250101)
options(timeout = 600)

## ---- 2. Dependencies ----
required_pkgs <- c(
  "curatedMetagenomicData", "dplyr", "tidyr", "limma", "glmnet", "mgcv",
  "randomForest", "ggplot2", "svglite", "haven", "SummarizedExperiment",
  "matrixStats", "patchwork", "Matrix"
)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = lib_dir, repos = "https://cloud.r-project.org")

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p, lib = lib_dir, update = FALSE, ask = FALSE)
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf("Package '%s' still unavailable after install.", p))
}

## ---- 3. Safe download helper (used by data-acquisition scripts) ----
# Downloads `url` to `dest` unless the file already exists with at least
# `min_size` bytes; rejects HTML error pages served in place of data.
safe_dl <- function(url, dest, min_size = 1024) {
  if (file.exists(dest) && file.info(dest)$size >= min_size)
    return(invisible(TRUE))
  tmp <- paste0(dest, ".part")
  tryCatch({
    download.file(url, tmp, mode = "wb", quiet = TRUE)
    sz <- file.info(tmp)$size
    if (is.na(sz) || sz < min_size) stop("file too small: ", sz, " bytes")
    first <- readLines(tmp, n = 1, warn = FALSE)
    if (grepl("^<(!DOCTYPE|html)", first, ignore.case = TRUE))
      stop("HTML error page instead of data")
    file.rename(tmp, dest)
    TRUE
  }, error = function(e) {
    unlink(tmp)
    message("FAILED: ", url, " | ", conditionMessage(e))
    FALSE
  })
}

## ---- 4. Logging ----
log_dir  <- file.path(res_dir, "logs")
log_file <- file.path(log_dir, "module0_setup.log")
dir.create(log_dir, showWarnings = FALSE)
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

## ---- 5. Session record ----
writeLines(capture.output(sessionInfo()),
           file.path(res_dir, "sessionInfo_setup.txt"))

cat("setup OK\n")
