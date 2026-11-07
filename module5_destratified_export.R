# ============================================================================
# module5_destratified_export.R — Per-pathway export of de-stratified
# shared pathways
#
# Project: Paired oral-gut microbiome aging clock (OG-Clock)
#          (module 5 / data source for Figure 6d)
#
# Input: <res_dir>/pillarC_oral_age_pathways.csv (significant stratified
#        features, oral)
#        <res_dir>/pillarC_gut_age_pathways.csv  (significant stratified
#        features, gut)
# Output: <res_dir>/destratified_shared_pathways.csv
#         (feature / logFC_oral / logFC_gut / same_dir)
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

EXCL_UNINTEGRATED <- TRUE

o <- read.csv(file.path(res_dir, "pillarC_oral_age_pathways.csv"),
              check.names = FALSE)
g <- read.csv(file.path(res_dir, "pillarC_gut_age_pathways.csv"),
              check.names = FALSE)

strip <- function(x) sub("\\|.*$", "", x)
o$pw <- strip(o$feature); g$pw <- strip(g$feature)
cat("oral sig features:", nrow(o), "-> unique pathways:", length(unique(o$pw)), "
")
cat("gut  sig features:", nrow(g), "-> unique pathways:", length(unique(g$pw)), "
")

ao <- tapply(o$logFC, o$pw, mean)
ag <- tapply(g$logFC, g$pw, mean)
common <- sort(intersect(names(ao), names(ag)))
cat("de-stratified shared:", length(common),
    "| contains UNINTEGRATED:", "UNINTEGRATED" %in% common, "
")

if (EXCL_UNINTEGRATED) {
  common <- setdiff(common, "UNINTEGRATED")
  cat("shared after dropping UNINTEGRATED:", length(common), "
")
}

same <- sign(ao[common]) == sign(ag[common])
rho  <- suppressWarnings(cor(ao[common], ag[common], method = "spearman"))
bp   <- binom.test(sum(same), length(common), p = 0.5)$p.value
cat(sprintf("same direction: %d/%d | Spearman %.3f | binom P = %.2e
",
            sum(same), length(common), rho, bp))

out <- data.frame(
  feature    = common,
  logFC_oral = as.numeric(ao[common]),
  logFC_gut  = as.numeric(ag[common]),
  same_dir   = same,
  stringsAsFactors = FALSE
)
write.csv(out, file.path(res_dir, "destratified_shared_pathways.csv"),
          row.names = FALSE)
cat("wrote destratified_shared_pathways.csv (", nrow(out), "rows)
")
