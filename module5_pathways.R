# ============================================================================
# module5_pathways.R — Pathway-level age associations + reverse validation
#
# Project: Paired oral-gut microbiome aging clock (OG-Clock)
#
# Input: cMD3 pathway abundances (BritoIL_2016 + HMP_2012)
# Output:
#   <data_dir>/pillarC_pathways.rds
#   <res_dir>/pillarC_oral_age_pathways.csv
#   <res_dir>/pillarC_gut_age_pathways.csv
#   <res_dir>/pillarC_reverse_validation.csv
#
# Design:
#   - Shared helpers from module5_utils.R
#   - CLR with within-sample multiplicative replacement
#   - prev >= 0.05 for pathways (the 0.2 species threshold is too strict
#     for sparse pathway abundances)
#   - limma with study_name covariate + sensitivity analysis without it
#   - Bidirectional reverse validation on limma results
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))
source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module5_utils.R"))

set.seed(20250101)

## ============================================================
## 1. Load pathway data
## ============================================================
cat("
========== 1. Loading pathway data ==========
")

studies <- c("BritoIL_2016", "HMP_2012")

# Pathways are sparse; use a lower prevalence threshold
oral_pw <- load_paired_view(studies, "pathway_abundance", "oralcavity",
                            prev_thresh = 0.05)
gut_pw  <- load_paired_view(studies, "pathway_abundance", "stool",
                            prev_thresh = 0.05)

cat("
Oral pathways:", paste(dim(oral_pw$X), collapse = " x "), "
")
cat("Gut pathways:", paste(dim(gut_pw$X), collapse = " x "), "
")

## ============================================================
## 2. Primary analysis: limma age associations (study_name adjusted)
## ============================================================
cat("
========== 2. Primary analysis (limma + study_name) ==========
")

oral_res <- limma_age_assoc(oral_pw$X, oral_pw$meta, adjust_study = TRUE)
gut_res  <- limma_age_assoc(gut_pw$X,  gut_pw$meta,  adjust_study = TRUE)

oral_res$niche <- "oral"
gut_res$niche  <- "gut"

cat("Oral significant (FDR<0.05):",
    sum(oral_res$adj.P.Val < 0.05), "/", nrow(oral_res), "
")
cat("Gut significant (FDR<0.05):",
    sum(gut_res$adj.P.Val < 0.05), "/", nrow(gut_res), "
")

## ============================================================
## 3. Sensitivity analysis: without study_name
## ============================================================
cat("
========== 3. Sensitivity analysis (no study_name) ==========
")

oral_res_ns <- limma_age_assoc(oral_pw$X, oral_pw$meta, adjust_study = FALSE)
gut_res_ns  <- limma_age_assoc(gut_pw$X,  gut_pw$meta,  adjust_study = FALSE)

# Overlap between the two analyses
overlap_oral <- length(intersect(
  oral_res$feature[oral_res$adj.P.Val < 0.05],
  oral_res_ns$feature[oral_res_ns$adj.P.Val < 0.05]
))
overlap_gut <- length(intersect(
  gut_res$feature[gut_res$adj.P.Val < 0.05],
  gut_res_ns$feature[gut_res_ns$adj.P.Val < 0.05]
))

cat("Oral significant-pathway overlap:", overlap_oral, "
")
cat("Gut significant-pathway overlap:", overlap_gut, "
")

## ============================================================
## 4. Reverse validation
## ============================================================
cat("
========== 4. Reverse validation ==========
")

# Oral discoveries -> tested in gut
rv_oral_to_gut <- reverse_validation(oral_res, gut_res,
                                     feature_col = "feature")

cat("
=== Oral discovery -> gut test ===
")
cat("Oral significant:", rv_oral_to_gut$n_sig, "
")
cat("Testable in gut:", rv_oral_to_gut$n_testable, "
")
cat("Same direction:", rv_oral_to_gut$n_same_dir, "/",
    rv_oral_to_gut$n_testable, "
")
cat("Binomial p:", format(rv_oral_to_gut$binom_p, digits = 3), "
")
cat("Gut FDR<0.05:", rv_oral_to_gut$n_dst_sig, "
")
cat("logFC Spearman:", round(rv_oral_to_gut$spearman_lfc, 3), "
")

# Gut discoveries -> tested in oral
rv_gut_to_oral <- reverse_validation(gut_res, oral_res,
                                     feature_col = "feature")

cat("
=== Gut discovery -> oral test ===
")
cat("Gut significant:", rv_gut_to_oral$n_sig, "
")
cat("Testable in oral:", rv_gut_to_oral$n_testable, "
")
cat("Same direction:", rv_gut_to_oral$n_same_dir, "/",
    rv_gut_to_oral$n_testable, "
")
cat("Binomial p:", format(rv_gut_to_oral$binom_p, digits = 3), "
")
cat("Oral FDR<0.05:", rv_gut_to_oral$n_dst_sig, "
")
cat("logFC Spearman:", round(rv_gut_to_oral$spearman_lfc, 3), "
")

## ============================================================
## 5. logFC concordance across all shared pathways
## ============================================================
cat("
========== 5. Shared-pathway concordance ==========
")

shared_feat <- intersect(oral_res$feature, gut_res$feature)
cat("Shared pathways:", length(shared_feat), "
")

sh <- data.frame(
  feature = shared_feat,
  logFC_oral = oral_res$logFC[match(shared_feat, oral_res$feature)],
  logFC_gut  = gut_res$logFC[match(shared_feat, gut_res$feature)],
  fdr_oral   = oral_res$adj.P.Val[match(shared_feat, oral_res$feature)],
  fdr_gut    = gut_res$adj.P.Val[match(shared_feat, gut_res$feature)]
)

rho_all <- suppressWarnings(cor(sh$logFC_oral, sh$logFC_gut,
                                method = "spearman"))
cat("All shared pathways, logFC Spearman:", round(rho_all, 3), "
")

# Significant in both niches
both_sig <- sh %>% filter(fdr_oral < 0.05 & fdr_gut < 0.05)
cat("Significant in both:", nrow(both_sig), "
")

if (nrow(both_sig) > 2) {
  rho_both <- suppressWarnings(cor(both_sig$logFC_oral, both_sig$logFC_gut,
                                   method = "spearman"))
  cat("Significant in both, logFC Spearman:", round(rho_both, 3), "
")
}

## ============================================================
## 6. Hypergeometric overlap test
## ============================================================
cat("
========== 6. Hypergeometric overlap test ==========
")

o_sig <- oral_res$feature[oral_res$adj.P.Val < 0.05]
g_sig <- gut_res$feature[gut_res$adj.P.Val < 0.05]
universe <- shared_feat

shared_sig <- intersect(o_sig, g_sig)

p_hyper <- phyper(
  length(shared_sig) - 1,
  length(intersect(o_sig, universe)),
  length(universe) - length(intersect(o_sig, universe)),
  length(intersect(g_sig, universe)),
  lower.tail = FALSE
)

cat("Oral sig intersect gut sig:", length(shared_sig), "
")
cat("Universe (shared pathways):", length(universe), "
")
cat("Hypergeometric p:", format(p_hyper, digits = 3), "
")

## ============================================================
## 7. Save
## ============================================================
out <- list(
  oral           = oral_res,
  gut            = gut_res,
  oral_ns        = oral_res_ns,
  gut_ns         = gut_res_ns,
  rv_oral_to_gut = rv_oral_to_gut,
  rv_gut_to_oral = rv_gut_to_oral,
  shared         = sh,
  shared_sig     = shared_sig,
  p_hyper        = p_hyper,
  rho_all        = rho_all,
  n_oral_sig     = sum(oral_res$adj.P.Val < 0.05),
  n_gut_sig      = sum(gut_res$adj.P.Val < 0.05),
  overlap_oral   = overlap_oral,
  overlap_gut    = overlap_gut
)

saveRDS(out, file.path(data_dir, "pillarC_pathways.rds"))

# Significant pathway lists
write.csv(
  oral_res %>% filter(adj.P.Val < 0.05) %>% arrange(P.Value),
  file.path(res_dir, "pillarC_oral_age_pathways.csv"),
  row.names = FALSE
)
write.csv(
  gut_res %>% filter(adj.P.Val < 0.05) %>% arrange(P.Value),
  file.path(res_dir, "pillarC_gut_age_pathways.csv"),
  row.names = FALSE
)
write.csv(
  rv_oral_to_gut$table,
  file.path(res_dir, "pillarC_reverse_validation.csv"),
  row.names = FALSE
)

## ============================================================
## 8. Summary
## ============================================================
cat("
========== module5_pathways DONE ==========
")
cat("Oral significant:", out$n_oral_sig, "
")
cat("Gut significant:", out$n_gut_sig, "
")
cat("Oral discovery -> testable in gut:", rv_oral_to_gut$n_testable, "
")
cat("  same direction:", rv_oral_to_gut$n_same_dir, "
")
cat("  binomial p:", format(rv_oral_to_gut$binom_p, digits = 3), "
")
cat("Gut discovery -> testable in oral:", rv_gut_to_oral$n_testable, "
")
cat("  same direction:", rv_gut_to_oral$n_same_dir, "
")
cat("  binomial p:", format(rv_gut_to_oral$binom_p, digits = 3), "
")
