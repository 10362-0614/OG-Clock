# ============================================================================
# module5_marginal_contribution.R — Marginal contribution quantification
#
# Project: Paired oral-gut microbiome aging clock (OG-Clock)
#
# Input: cMD3 oral + gut (species + pathways)
# Output:
#   <data_dir>/marginal_contribution.rds
#   <res_dir>/marginal_contribution.csv
#
# Design:
#   - Four views: oral species, gut species, oral pathways, gut pathways
#     (species prev >= 0.2, pathways prev >= 0.05)
#   - Subject-level alignment across views (oral and gut are different
#     samples; alignment must go through subject_id)
#   - 5x5 repeated CV per view (oral / gut / joint), ENet on top-150
#     |Spearman| features
#   - Subject-level bootstrap CI for delta-MAE (B = 2000)
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))
source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module5_utils.R"))

suppressPackageStartupMessages({
  library(glmnet)
})

set.seed(20250101)

## ============================================================
## 0. run_enet (seeded)
## ============================================================
run_enet <- function(X, y, train_idx, test_idx,
                     nfeat = 150, alpha = 0.5,
                     lambda_choice = "lambda.min",
                     seed = 20250101) {
  set.seed(seed)
  
  X_tr <- X[train_idx, , drop = FALSE]
  y_tr <- y[train_idx]
  X_te <- X[test_idx, , drop = FALSE]
  
  rho <- suppressWarnings(
    apply(X_tr, 2, function(x) cor(x, y_tr, method = "spearman",
                                   use = "complete.obs"))
  )
  rho <- rho[!is.na(rho)]
  n_keep <- min(nfeat, length(rho))
  if (n_keep < 5) return(rep(NA_real_, length(test_idx)))
  keep <- names(sort(abs(rho), decreasing = TRUE))[seq_len(n_keep)]
  
  X_tr_f <- X_tr[, keep, drop = FALSE]
  X_te_f <- X_te[, keep, drop = FALSE]
  X_tr_f[!is.finite(X_tr_f)] <- 0
  X_te_f[!is.finite(X_te_f)] <- 0
  
  fit <- cv.glmnet(X_tr_f, y_tr, alpha = alpha, nfolds = 5)
  as.numeric(predict(fit, X_te_f, s = fit[[lambda_choice]]))
}

## ============================================================
## 1. Load the four views
## ============================================================
cat("
========== 1. Loading views ==========
")

studies <- c("BritoIL_2016", "HMP_2012")

oral_sp <- load_paired_view(studies, "relative_abundance", "oralcavity",
                            prev_thresh = 0.2)
gut_sp  <- load_paired_view(studies, "relative_abundance", "stool",
                            prev_thresh = 0.2)
oral_pw <- load_paired_view(studies, "pathway_abundance", "oralcavity",
                            prev_thresh = 0.05)
gut_pw  <- load_paired_view(studies, "pathway_abundance", "stool",
                            prev_thresh = 0.05)

## ============================================================
## 2. Align subjects
## ============================================================
cat("
========== 2. Aligning subjects ==========
")

# Oral and gut are different samples, so the sample-id intersection is
# empty by construction; alignment must use subject_id (each view already
# keeps one sample per subject)
map_subj <- function(v) {
  i <- match(rownames(v$X), v$meta$sample_id)
  setNames(v$meta$subject_id[i], rownames(v$X))
}

subj_maps <- list(oral_sp = map_subj(oral_sp), gut_sp = map_subj(gut_sp),
                  oral_pw = map_subj(oral_pw), gut_pw  = map_subj(gut_pw))

common_subj <- sort(Reduce(intersect, lapply(subj_maps, unique)))
cat("Common subjects:", length(common_subj), "
")
stopifnot(length(common_subj) >= 50)

# Align each view to the common subjects (one row per subject,
# rownames = subject_id)
align_subj <- function(v, smap) {
  df <- data.frame(sid = names(smap), subj = unname(smap),
                   stringsAsFactors = FALSE)
  df <- df[df$subj %in% common_subj, ]
  df <- df[!duplicated(df$subj), ]          # one sample per subject
  df <- df[match(common_subj, df$subj), ]
  X <- v$X[match(df$sid, rownames(v$X)), , drop = FALSE]
  rownames(X) <- df$subj
  X
}

O_sp <- align_subj(oral_sp, subj_maps$oral_sp)
G_sp <- align_subj(gut_sp,  subj_maps$gut_sp)
O_pw <- align_subj(oral_pw, subj_maps$oral_pw)
G_pw <- align_subj(gut_pw,  subj_maps$gut_pw)

stopifnot(identical(rownames(O_sp), common_subj),
          identical(rownames(G_sp), common_subj),
          identical(rownames(O_pw), common_subj),
          identical(rownames(G_pw), common_subj))

# Age / study (subject level; subject ids carry the study prefix, so they
# cannot collide across studies)
age_of <- setNames(oral_sp$meta$age, oral_sp$meta$subject_id)
y <- unname(age_of[common_subj])
stopifnot(!any(is.na(y)))

study_vec <- sub("_.*$", "", common_subj)

X_oral <- cbind(O_sp, O_pw)
X_gut  <- cbind(G_sp, G_pw)

colnames(X_oral) <- c(paste0("O_",  colnames(O_sp)),
                      paste0("OP_", colnames(O_pw)))
colnames(X_gut)  <- c(paste0("G_",  colnames(G_sp)),
                      paste0("GP_", colnames(G_pw)))

cat("Oral features:", ncol(X_oral), "
")
cat("Gut features:", ncol(X_gut), "
")

## ============================================================
## 3. 5x5 repeated CV
## ============================================================
cat("
========== 3. 5x5 repeated CV ==========
")

views <- list(
  oral  = X_oral,
  gut   = X_gut,
  joint = cbind(X_oral, X_gut)
)

n <- length(y)
nrep <- 5; k <- 5

# Store predictions per repetition and per fold (needed for CIs)
preds_array <- lapply(views, function(v) {
  matrix(NA_real_, nrow = n, ncol = nrep)
})

for (r in 1:nrep) {
  set.seed(20250101 + r)
  folds <- sample(rep(1:k, length.out = n))
  
  for (f in 1:k) {
    te <- which(folds == f); tr <- which(folds != f)
    for (nm in names(views)) {
      preds_array[[nm]][te, r] <- run_enet(
        views[[nm]], y, tr, te,
        seed = 20250101 + r * 100 + f
      )
    }
  }
}

# Average over the 5 repetitions
preds <- lapply(preds_array, function(pp) rowMeans(pp, na.rm = TRUE))

## ============================================================
## 4. Performance metrics
## ============================================================
cat("
========== 4. Performance metrics ==========
")

mae <- function(p) mean(abs(p - y), na.rm = TRUE)
rho <- function(p) suppressWarnings(cor(p, y, method = "spearman",
                                        use = "complete.obs"))
r2  <- function(p) {
  ok <- !is.na(p)
  1 - sum((y[ok] - p[ok])^2) / sum((y[ok] - mean(y[ok]))^2)
}

res <- data.frame(
  model    = names(views),
  MAE      = sapply(preds, mae),
  Spearman = sapply(preds, rho),
  R2       = sapply(preds, r2)
)

cat("CV performance:
")
print(res, row.names = FALSE, digits = 3)

## ============================================================
## 5. Bootstrap CI (subject level)
## ============================================================
cat("
========== 5. Bootstrap ==========
")

B <- 2000
set.seed(20250101)

# Delta-MAE per bootstrap replicate
boot_delta <- function(pA, pB) {
  d <- sapply(1:B, function(i) {
    idx <- sample.int(n, n, replace = TRUE)
    mean(abs(pA[idx] - y[idx]), na.rm = TRUE) -
      mean(abs(pB[idx] - y[idx]), na.rm = TRUE)
  })
  quantile(d, c(0.025, 0.975), na.rm = TRUE)
}

# joint - oral
ci_joint_oral <- boot_delta(preds$joint, preds$oral)
# joint - gut
ci_joint_gut  <- boot_delta(preds$joint, preds$gut)

delta_joint_oral <- mae(preds$joint) - mae(preds$oral)
delta_joint_gut  <- mae(preds$joint) - mae(preds$gut)

cat(sprintf("Delta MAE (joint - oral): %+.2f y [%.2f, %.2f]
",
            delta_joint_oral, ci_joint_oral[1], ci_joint_oral[2]))
cat(sprintf("Delta MAE (joint - gut):  %+.2f y [%.2f, %.2f]
",
            delta_joint_gut, ci_joint_gut[1], ci_joint_gut[2]))

## ============================================================
## 6. Incremental R2 and partial correlations
## ============================================================
cat("
========== 6. Incremental R2 and partial correlations ==========
")

r2_o  <- summary(lm(y ~ preds$oral))$r.squared
r2_g  <- summary(lm(y ~ preds$gut))$r.squared
r2_og <- summary(lm(y ~ preds$oral + preds$gut))$r.squared

cat("R2 oral:", round(r2_o, 3), "
")
cat("R2 gut: ", round(r2_g, 3), "
")
cat("R2 joint (oral + gut):", round(r2_og, 3), "
")
cat("Incremental (oral -> joint):", round(r2_og - r2_o, 3), "
")
cat("Incremental (gut  -> joint):", round(r2_og - r2_g, 3), "
")

# Partial correlations
partial_r <- function(a, b, y) {
  ra <- residuals(lm(a ~ b))
  ry <- residuals(lm(y ~ b))
  cor(ra, ry)
}
pc_oral_given_gut <- partial_r(preds$oral, preds$gut, y)
pc_gut_given_oral <- partial_r(preds$gut, preds$oral, y)

cat("Partial r (oral | gut):", round(pc_oral_given_gut, 3), "
")
cat("Partial r (gut | oral):", round(pc_gut_given_oral, 3), "
")

## ============================================================
## 7. Save
## ============================================================
out <- list(
  res              = res,
  preds            = preds,
  preds_array      = preds_array,
  y                = y,
  study            = study_vec,
  ci_joint_oral    = ci_joint_oral,
  ci_joint_gut     = ci_joint_gut,
  delta_joint_oral = delta_joint_oral,
  delta_joint_gut  = delta_joint_gut,
  r2               = c(oral = r2_o, gut = r2_g, joint = r2_og),
  pc_oral_given_gut = pc_oral_given_gut,
  pc_gut_given_oral = pc_gut_given_oral,
  n_subj = n
)

saveRDS(out, file.path(data_dir, "marginal_contribution.rds"))

# Summary CSV
summary_df <- data.frame(
  metric = c("MAE_oral", "MAE_gut", "MAE_joint",
             "delta_MAE_joint_oral", "delta_MAE_joint_oral_ci_low",
             "delta_MAE_joint_oral_ci_high",
             "delta_MAE_joint_gut", "delta_MAE_joint_gut_ci_low",
             "delta_MAE_joint_gut_ci_high",
             "R2_oral", "R2_gut", "R2_joint",
             "incr_R2_oral_to_joint", "incr_R2_gut_to_joint",
             "partial_cor_oral_given_gut",
             "partial_cor_gut_given_oral"),
  value = c(res$MAE[1], res$MAE[2], res$MAE[3],
            delta_joint_oral, ci_joint_oral[1], ci_joint_oral[2],
            delta_joint_gut,  ci_joint_gut[1],  ci_joint_gut[2],
            r2_o, r2_g, r2_og,
            r2_og - r2_o, r2_og - r2_g,
            pc_oral_given_gut, pc_gut_given_oral)
)

write.csv(summary_df,
          file.path(res_dir, "marginal_contribution.csv"),
          row.names = FALSE)

## ============================================================
## 8. Summary
## ============================================================
cat("
========== module5_marginal_contribution DONE ==========
")
cat("Outputs:
")
cat("  ", file.path(data_dir, "marginal_contribution.rds"), "
")
cat("  ", file.path(res_dir, "marginal_contribution.csv"), "
")
