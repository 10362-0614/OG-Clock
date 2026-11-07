# ============================================================================
# module3_clock_transfer.R — HMP2 longitudinal clock (cMD3 transfer + native
# out-of-fold predictions)
#
# Input:
#   <data_dir>/hmp2_species_matrix.rds
#   <data_dir>/hmp2_sample_meta.rds
#   <data_dir>/BritoIL_2016__relative_abundance.rds
#   <data_dir>/HMP_2012__relative_abundance.rds
#   <data_dir>/sample_metadata.rds
#
# Output:
#   <data_dir>/hmp2_clock_resid.rds
#
# Design:
#   - CLR with within-sample multiplicative replacement (as in modules 1-2)
#   - cMD3 shared-species clock trained on stool, transferred to HMP2
#   - HMP2-native clock evaluated with 5x5 repeated CV (not a single
#     5-fold split)
#   - Age estimate: consent_age + week_num / 52
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(parallel)
})

set.seed(2026)

## ============================================================
## 0. CLR (within-sample multiplicative replacement)
## ============================================================
# Identical to modules 1-2
clr_mat <- function(m, pseudocount_factor = 0.5) {
  # m: samples x features
  m2 <- m
  for (i in seq_len(nrow(m))) {
    row_i <- m[i, ]
    pos <- row_i[row_i > 0]
    if (length(pos) == 0) next
    pmin_i <- min(pos)
    zero_idx <- which(row_i == 0 | is.na(row_i))
    if (length(zero_idx) > 0) {
      m2[i, zero_idx] <- pmin_i * pseudocount_factor
    }
  }
  row_sums <- rowSums(m2)
  m2 <- sweep(m2, 1, row_sums, "/")
  log_m <- log(m2)
  log_m - rowMeans(log_m)
}

## ============================================================
## 1. Load data
## ============================================================
cat("
========== 1. Load data ==========
")

hmp2 <- readRDS(file.path(data_dir, "hmp2_species_matrix.rds"))
md   <- readRDS(file.path(data_dir, "hmp2_sample_meta.rds"))

cat("HMP2:", paste(dim(hmp2), collapse = " x "), "
")

## ============================================================
## 2. cMD3 stool training set
## ============================================================
cat("
========== 2. cMD3 training set ==========
")

get_mat <- function(x) if (is.list(x)) x$mat else x

sp1 <- get_mat(readRDS(file.path(data_dir, "BritoIL_2016__relative_abundance.rds")))
sp2 <- get_mat(readRDS(file.path(data_dir, "HMP_2012__relative_abundance.rds")))

common_sp <- intersect(rownames(sp1), rownames(sp2))
sp <- cbind(sp1[common_sp, , drop = FALSE], sp2[common_sp, , drop = FALSE])

md3 <- readRDS(file.path(data_dir, "sample_metadata.rds"))

stool <- md3 %>%
  filter(body_site == "stool",
         study_name %in% c("BritoIL_2016", "HMP_2012"),
         age > 0) %>%
  group_by(study_name, subject_id) %>%
  slice_head(n = 1) %>%
  ungroup()

sel3 <- intersect(stool$sample_id, colnames(sp))
X3_raw <- t(sp[common_sp, sel3, drop = FALSE])   # samples x species
y3     <- stool$age[match(sel3, stool$sample_id)]

cat("cMD3 stool:", nrow(X3_raw), "samples
")
cat("Age range:", paste(range(y3), collapse = "-"), "
")

## ============================================================
## 3. Shared species (name-level)
## ============================================================
cat("
========== 3. Shared species ==========
")

# Strip the full lineage down to the species name; the leading ^s__
# sub() guards against names that retained the prefix
cmd3_name <- sub(".*\\|s__", "", common_sp)
cmd3_name <- sub("^s__", "", cmd3_name)

hmp2_sp <- rownames(hmp2)
shared_names <- intersect(unique(cmd3_name), hmp2_sp)

cat("cMD3 species:", length(common_sp), "-> species names:", length(unique(cmd3_name)), "
")
cat("HMP2 species:", length(hmp2_sp), "
")
cat("Shared species names:", length(shared_names), "
")

if (length(shared_names) < 50) {
  warning("Fewer than 50 shared species; training may be difficult")
}

## ============================================================
## 4. Align matrices
## ============================================================
cat("
========== 4. Alignment ==========
")

# cMD3: aggregate to species names
X3n <- t(rowsum(t(X3_raw), cmd3_name, reorder = FALSE))
X3s <- X3n[, shared_names, drop = FALSE]

# HMP2: extract shared species
Xh <- t(hmp2[shared_names, , drop = FALSE])   # samples x species

cat("X3s:", paste(dim(X3s), collapse = " x "), "
")
cat("Xh:",  paste(dim(Xh),  collapse = " x "), "
")

# CLR
C3 <- clr_mat(X3s)
NH <- clr_mat(Xh)

cat("After CLR:
")
cat("  C3:", paste(dim(C3), collapse = " x "), "
")
cat("  NH:", paste(dim(NH), collapse = " x "), "
")

## ============================================================
## 5. Training function
## ============================================================
run_enet <- function(X, y, tr, te, nfeat = 150, alpha = 0.5,
                     lambda_choice = "lambda.min", seed = 2026) {
  set.seed(seed)
  
  rho <- suppressWarnings(
    apply(X[tr, , drop = FALSE], 2,
          function(x) cor(x, y[tr], method = "spearman",
                          use = "complete.obs"))
  )
  rho <- rho[!is.na(rho)]
  n_keep <- min(nfeat, length(rho))
  if (n_keep < 5) return(rep(NA_real_, length(te)))
  sel <- names(sort(abs(rho), decreasing = TRUE))[seq_len(n_keep)]
  
  X_tr <- X[tr, sel, drop = FALSE]
  X_te <- X[te, sel, drop = FALSE]
  X_tr[!is.finite(X_tr)] <- 0
  X_te[!is.finite(X_te)] <- 0
  
  fit <- cv.glmnet(X_tr, y[tr], alpha = alpha, nfolds = 5)
  as.numeric(predict(fit, X_te, s = fit[[lambda_choice]]))
}

## ============================================================
## 6. cMD3 shared-species clock: 5x5 CV
## ============================================================
cat("
========== 6. cMD3 5x5 CV ==========
")

n3 <- length(y3)
preds3 <- matrix(NA_real_, nrow = n3, ncol = 5)
for (r in 1:5) {
  folds <- sample(rep(1:5, length.out = n3))
  for (f in 1:5) {
    te <- which(folds == f); tr <- which(folds != f)
    preds3[te, r] <- run_enet(C3, y3, tr, te,
                              seed = 2026 + r * 100 + f)
  }
}
pred3_mean <- rowMeans(preds3, na.rm = TRUE)

mae3 <- mean(abs(pred3_mean - y3))
rho3 <- suppressWarnings(cor(pred3_mean, y3, method = "spearman"))
cat(sprintf("cMD3 5x5 CV: MAE=%.2f rho=%.3f (n=%d)
", mae3, rho3, n3))

## ============================================================
## 7. Transfer: train on full cMD3 -> apply to HMP2
## ============================================================
cat("
========== 7. Transfer ==========
")

# Train on full cMD3
rho_all <- suppressWarnings(
  apply(C3, 2, function(x) cor(x, y3, method = "spearman",
                               use = "complete.obs"))
)
rho_all <- rho_all[!is.na(rho_all)]
sel_final <- names(sort(abs(rho_all), decreasing = TRUE))[
  seq_len(min(150, length(rho_all)))
]

fit_final <- cv.glmnet(C3[, sel_final, drop = FALSE], y3,
                       alpha = 0.5, nfolds = 5)
pred_h <- as.numeric(predict(fit_final,
                             NH[, sel_final, drop = FALSE],
                             s = "lambda.min"))

# Age estimates
md$age_est <- md$consent_age + md$week_num / 52
md$pred    <- pred_h[match(md$sid, rownames(NH))]

# Transfer diagnostic
keep_pred <- !is.na(md$age_est) & !is.na(md$pred)
rho_transfer <- suppressWarnings(
  cor(md$pred[keep_pred], md$age_est[keep_pred], method = "spearman")
)
cat(sprintf("Transfer pred-age rho = %.3f
", rho_transfer))
cat("Note: rho < 0.3 indicates failed transfer
")

md <- md %>% filter(keep_pred)
md$resid <- md$pred - md$age_est

# Recalibrated version
cal <- lm(age_est ~ pred, data = md)
md$resid_cal <- md$age_est - fitted(cal)

cat(sprintf("HMP2 samples: %d | subjects: %d
",
            nrow(md), length(unique(md[["Participant ID"]]))))
cat(sprintf("resid SD: %.2f | resid_cal SD: %.2f
",
            sd(md$resid), sd(md$resid_cal)))

## ============================================================
## 8. HMP2-native clock: 5x5 repeated OOF
## ============================================================
cat("
========== 8. HMP2-native 5x5 OOF ==========
")

Xh_full <- NH
y_h     <- md$age_est[match(rownames(NH), md$sid)]

keep_y <- !is.na(y_h)
Xh_k <- Xh_full[keep_y, , drop = FALSE]
y_hk <- y_h[keep_y]

cat("Training samples:", nrow(Xh_k), "
")

n_h <- nrow(Xh_k)
oof_rep <- matrix(NA_real_, nrow = n_h, ncol = 5)
for (r in 1:5) {
  folds <- sample(rep(1:5, length.out = n_h))
  for (f in 1:5) {
    te <- which(folds == f); tr <- which(folds != f)
    oof_rep[te, r] <- run_enet(Xh_k, y_hk, tr, te,
                               seed = 2026 + r * 100 + f)
  }
}
oof <- rowMeans(oof_rep, na.rm = TRUE)

mae_native <- mean(abs(oof - y_hk))
rho_native <- suppressWarnings(cor(oof, y_hk, method = "spearman"))
cat(sprintf("HMP2-native 5x5 OOF: MAE=%.2f rho=%.3f
",
            mae_native, rho_native))

# Assign back
md$oof_resid <- NA_real_
idx_k <- rownames(Xh_k)
md$oof_resid[match(idx_k, md$sid)] <- oof - y_hk

## ============================================================
## 9. Save
## ============================================================
out <- list(
  md = md,
  cv_perf = c(MAE = mae3, Spearman = rho3),
  native_perf = c(MAE = mae_native, Spearman = rho_native),
  transfer_rho = rho_transfer,
  n_shared = length(shared_names),
  shared_names = shared_names
)

saveRDS(out, file.path(data_dir, "hmp2_clock_resid.rds"))

cat("
========== module3_clock_transfer DONE ==========
")
cat("Output:", file.path(data_dir, "hmp2_clock_resid.rds"), "
")
cat("Key numbers:
")
cat("  Transfer rho:", round(rho_transfer, 3), "(failure threshold 0.3)
")
cat("  Native 5x5 CV:", sprintf("MAE=%.2f rho=%.3f", mae_native, rho_native), "
")
cat("  Shared species:", length(shared_names), "
")
