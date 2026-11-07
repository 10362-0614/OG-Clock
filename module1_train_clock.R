# ============================================================================
# module1_train_clock.R — Train oral, gut, and joint microbiome aging clocks
# on the paired cMD3 training set (BritoIL_2016 + HMP_2012).
#
# Design (matches Methods):
#   - Views = species + pathway abundance concatenated
#     (X_oral = cbind(O: species, OP: pathways); gut analogously G:/GP:)
#   - CLR with multiplicative replacement (zeros -> per-column pos_min/2,
#     column re-normalization, within-sample CLR)
#   - In-fold Spearman |rho| top-150 feature prefilter
#   - cv.glmnet(alpha = 0.5, nfolds = 5) at lambda.min
#   - 5-fold x 5-repetition CV with overwrite-style out-of-fold prediction
#     (all_pred[te] <- pred; not fold-averaged)
#   - Leave-one-study-out (LOSO) cross-cohort extrapolation
#   - Complementarity decomposition (bootstrap B = 2000, partial correlations)
# Execution order: set.seed(42) -> CV (oral, gut, joint) -> LOSO ->
#   complementarity
#
# Output:
#   <data_dir>/03_models.rds, 03_predictions.rds
#   <res_dir>/03_performance.csv, 03_complementarity.csv
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(matrixStats)
  library(glmnet)
})

set.seed(42)

## ============================================================
## 1. Load paired data
## ============================================================
cat("\n========== Loading paired data ==========\n")

load_mat <- function(study, typ) {
  readRDS(file.path(data_dir, paste0(study, "__", typ, ".rds")))
}

mk_subject_table <- function(study, typ) {
  x <- load_mat(study, typ)
  m <- x$mat
  # Use the native row order of sample_metadata.rds rather than the $meta
  # stored in the study RDS (matrix column order) — slice_head() sampling
  # depends on row order
  pd <- readRDS(file.path(data_dir, "sample_metadata.rds")) %>%
    filter(study_name == study)
  pd$sample_id <- as.character(pd$sample_id)
  
  # Prefix subject IDs with the study name (idempotent)
  pd$subject_id <- as.character(pd$subject_id)
  prefix <- paste0(study, "_")
  need_prefix <- !startsWith(pd$subject_id, prefix)
  pd$subject_id[need_prefix] <- paste0(prefix, pd$subject_id[need_prefix])
  
  keep_cols <- c("sample_id", "subject_id", "body_site", "age",
                 "study_name", "sex", "gender")
  pd <- pd[, intersect(keep_cols, colnames(pd)), drop = FALSE]
  
  if (!"sex" %in% colnames(pd) && "gender" %in% colnames(pd))
    pd$sex <- pd$gender
  if (!"study_name" %in% colnames(pd)) pd$study_name <- study
  
  list(mat = m, meta = pd)
}

# One subject per body site (first sample per subject)
build_subject_matrix <- function(study, typ, site) {
  d <- mk_subject_table(study, typ)
  sub <- d$meta %>% filter(body_site == site) %>%
    group_by(subject_id) %>% slice_head(n = 1) %>% ungroup()
  m <- d$mat[, sub$sample_id, drop = FALSE]
  list(mat = m, info = sub)
}

# Features x samples; multiplicative replacement then CLR
clr <- function(m) {
  pos_min <- apply(m, 2, function(x) min(x[x > 0], na.rm = TRUE))
  m2 <- m
  for (j in seq_len(ncol(m))) m2[is.na(m2[, j]) | m2[, j] == 0, j] <- pos_min[j] / 2
  m2 <- sweep(m2, 2, colSums(m2), "/")
  apply(m2, 2, function(x) log(x) - mean(log(x)))
}

prep_view <- function(studies, typ, site) {
  mats <- lapply(studies, function(s) build_subject_matrix(s, typ, site))
  common_feat <- Reduce(intersect, lapply(mats, function(d) rownames(d$mat)))
  common_feat <- sort(common_feat)
  mats <- lapply(mats, function(d) { d$mat <- d$mat[common_feat, , drop = FALSE]; d })
  X <- do.call(cbind, lapply(mats, `[[`, "mat"))
  info <- do.call(rbind, lapply(mats, `[[`, "info"))
  Xc <- clr(X)
  rownames(Xc) <- common_feat
  list(X = t(Xc), y = info$age, subject = info$subject_id, study = info$study_name)
}

studies <- c("BritoIL_2016", "HMP_2012")
oral_sp  <- prep_view(studies, "relative_abundance", "oralcavity")
gut_sp   <- prep_view(studies, "relative_abundance", "stool")
oral_pw  <- prep_view(studies, "pathway_abundance",  "oralcavity")
gut_pw   <- prep_view(studies, "pathway_abundance",  "stool")

cat(sprintf("subjects: oral_sp=%d gut_sp=%d | features: oral_sp=%d gut_sp=%d oral_pw=%d gut_pw=%d\n",
            nrow(oral_sp$X), nrow(gut_sp$X), ncol(oral_sp$X), ncol(gut_sp$X), ncol(oral_pw$X), ncol(gut_pw$X)))

# Align subjects across views (same subject set, same order)
common <- Reduce(intersect, list(oral_sp$subject, gut_sp$subject, oral_pw$subject, gut_pw$subject))
align <- function(v) { idx <- match(common, v$subject); list(X = v$X[idx, , drop = FALSE], y = v$y[idx]) }
o_sp <- align(oral_sp); g_sp <- align(gut_sp); o_pw <- align(oral_pw); g_pw <- align(gut_pw)
y <- o_sp$y
study_vec <- oral_sp$study[match(common, oral_sp$subject)]
cat("aligned subjects:", length(common), "| age range:", min(y), "-", max(y), "\n")

# Joint feature matrices (species + pathways)
X_oral  <- cbind(o_sp$X, o_pw$X); colnames(X_oral)  <- c(paste0("O:", colnames(o_sp$X)), paste0("OP:", colnames(o_pw$X)))
X_gut   <- cbind(g_sp$X, g_pw$X); colnames(X_gut)   <- c(paste0("G:", colnames(g_sp$X)), paste0("GP:", colnames(g_pw$X)))
X_joint <- cbind(X_oral, X_gut)

cat("\n========== View dimensions ==========\n")
cat("oral: ", nrow(X_oral),  " x ", ncol(X_oral),  "\n")
cat("gut:  ", nrow(X_gut),   " x ", ncol(X_gut),   "\n")
cat("joint:", nrow(X_joint), " x ", ncol(X_joint), "\n")

## ============================================================
## 2. Modeling functions
## ============================================================
run_enet <- function(X, y, train_idx, test_idx, nfeat = 150) {
  tr <- X[train_idx, , drop = FALSE]; te <- X[test_idx, , drop = FALSE]
  ytr <- y[train_idx]
  rho <- suppressWarnings(apply(tr, 2, function(x) cor(x, ytr, method = "spearman", use = "complete.obs")))
  keep <- names(sort(abs(rho), decreasing = TRUE))[seq_len(min(nfeat, sum(!is.na(rho))))]
  tr <- tr[, keep, drop = FALSE]; te <- te[, keep, drop = FALSE]
  fit <- glmnet::cv.glmnet(as.matrix(tr), ytr, alpha = 0.5, nfolds = 5)
  pred <- as.numeric(predict(fit, as.matrix(te), s = "lambda.min"))
  pred
}

metrics <- function(pred, obs) c(MAE = mean(abs(pred - obs)),
                                 Spearman = suppressWarnings(cor(pred, obs, method = "spearman")))

## ============================================================
## 3. Repeated 5x5 CV
## ============================================================
cv_eval <- function(X, label, nrep = 5, k = 5) {
  n <- length(y)
  all_pred <- rep(NA_real_, n)
  maes <- c()
  for (r in seq_len(nrep)) {
    folds <- sample(rep(1:k, length.out = n))
    for (f in 1:k) {
      te <- which(folds == f); tr <- which(folds != f)
      pred <- run_enet(X, y, tr, te)
      all_pred[te] <- pred
      maes <- c(maes, mean(abs(pred - y[te])))
    }
  }
  m <- metrics(all_pred, y)
  cat(sprintf("[%s] 5x%dCV  MAE=%.2f  rho=%.3f\n", label, nrep, m["MAE"], m["Spearman"]))
  list(pred = all_pred, mae = m["MAE"], rho = m["Spearman"], label = label)
}

cat("\n=== Repeated CV (elastic net, in-fold prefilter) ===\n")
res_oral  <- cv_eval(X_oral,  "oral-clock")
res_gut   <- cv_eval(X_gut,   "gut-clock")
res_joint <- cv_eval(X_joint, "joint-clock")

## ============================================================
## 4. Leave-one-study-out (LOSO)
## ============================================================
cat("\n=== Leave-one-study-out ===\n")
loso_eval <- function(X, label) {
  preds <- rep(NA_real_, length(y))
  for (s in unique(study_vec)) {
    te <- which(study_vec == s); tr <- which(study_vec != s)
    preds[te] <- run_enet(X, y, tr, te)
    cat(sprintf("  test on %s (n=%d): MAE=%.2f\n", s, length(te), mean(abs(preds[te] - y[te]))))
  }
  m <- metrics(preds, y)
  cat(sprintf("[%s] LOSO  MAE=%.2f  rho=%.3f\n", label, m["MAE"], m["Spearman"]))
  list(pred = preds, mae = m["MAE"], rho = m["Spearman"], label = paste0(label, "-LOSO"))
}
loso_oral  <- loso_eval(X_oral,  "oral-clock")
loso_gut   <- loso_eval(X_gut,   "gut-clock")
loso_joint <- loso_eval(X_joint, "joint-clock")

## ============================================================
## 5. Complementarity decomposition
## ============================================================
cat("\n=== Complementarity ===\n")
res_oral_o <- res_oral$pred - y; res_gut_o <- res_gut$pred - y; res_joint_o <- res_joint$pred - y
r_og <- suppressWarnings(cor(res_oral_o, res_gut_o, method = "spearman"))
cat(sprintf("residual correlation oral vs gut: %.3f\n", r_og))
cat(sprintf("deltaMAE joint vs best single: %.2f years\n",
            min(res_oral$mae, res_gut$mae) - res_joint$mae))

# Full complementarity: bootstrap CIs + partial correlations
set.seed(42)
po <- res_oral$pred; pg <- res_gut$pred; pj <- res_joint$pred
n <- length(y)
bs <- replicate(2000, {
  i <- sample.int(n, n, replace = TRUE)
  c(mean(abs(pj[i]-y[i])) - mean(abs(po[i]-y[i])),
    mean(abs(pj[i]-y[i])) - mean(abs(pg[i]-y[i])))
})
ci_o <- quantile(bs[1, ], c(0.025, 0.975), na.rm = TRUE)
ci_g <- quantile(bs[2, ], c(0.025, 0.975), na.rm = TRUE)

r2 <- function(p) summary(lm(y ~ p))$r.squared
r2_o <- r2(po); r2_g <- r2(pg); r2_j <- r2(pj)

pcor <- function(a, b, c) {
  ra <- residuals(lm(a ~ b)); rc <- residuals(lm(c ~ b))
  cor(ra, rc)
}
pc_oral_given_gut <- pcor(po, pg, y)
pc_gut_given_oral <- pcor(pg, po, y)

comp_table <- data.frame(
  metric = c(
    "residual_cor_spearman", "residual_cor_pearson",
    "deltaMAE_joint_minus_oral", "deltaMAE_joint_minus_oral_ci_low",
    "deltaMAE_joint_minus_oral_ci_high",
    "deltaMAE_joint_minus_gut", "deltaMAE_joint_minus_gut_ci_low",
    "deltaMAE_joint_minus_gut_ci_high",
    "R2_oral", "R2_gut", "R2_joint",
    "incr_R2_oral_to_joint", "incr_R2_gut_to_joint",
    "partial_cor_oral_given_gut", "partial_cor_gut_given_oral"
  ),
  value = c(
    r_og, cor(res_oral_o, res_gut_o, method = "pearson"),
    mean(abs(pj-y)) - mean(abs(po-y)), ci_o[1], ci_o[2],
    mean(abs(pj-y)) - mean(abs(pg-y)), ci_g[1], ci_g[2],
    r2_o, r2_g, r2_j,
    r2_j - r2_o, r2_j - r2_g,
    pc_oral_given_gut, pc_gut_given_oral
  )
)
cat("\nComplementarity metrics:\n"); print(comp_table, digits = 3)

## ============================================================
## 6. Performance summary table
## ============================================================
perf_table <- data.frame(
  clock      = rep(c("oral", "gut", "joint"), 2),
  validation = rep(c("5x5CV", "LOSO"), each = 3),
  n          = length(y),
  MAE        = c(res_oral$mae, res_gut$mae, res_joint$mae,
                 loso_oral$mae, loso_gut$mae, loso_joint$mae),
  rho        = c(res_oral$rho, res_gut$rho, res_joint$rho,
                 loso_oral$rho, loso_gut$rho, loso_joint$rho)
)
cat("\n========== Performance summary ==========\n")
print(perf_table, row.names = FALSE, digits = 3)

## ============================================================
## 7. Save
## ============================================================
results <- list(
  cv    = list(oral = res_oral, gut = res_gut, joint = res_joint),
  loso  = list(oral = loso_oral, gut = loso_gut, joint = loso_joint),
  complementarity = list(table = comp_table, rho_sp = r_og,
                         ci_o = ci_o, ci_g = ci_g,
                         r2 = c(oral = r2_o, gut = r2_g, joint = r2_j),
                         pcor = c(oral_given_gut = pc_oral_given_gut,
                                  gut_given_oral = pc_gut_given_oral)),
  perf  = perf_table,
  y = y, study = study_vec, subject = common
)
saveRDS(results, file.path(data_dir, "03_models.rds"))

predictions <- do.call(rbind, lapply(c("oral", "gut", "joint"), function(cl) {
  cv <- results$cv[[cl]]
  data.frame(
    clock      = cl,
    subject_id = common,
    study      = study_vec,
    age_true   = y,
    age_pred   = cv$pred,
    age_accel  = cv$pred - y
  )
}))
saveRDS(predictions, file.path(data_dir, "03_predictions.rds"))

write.csv(perf_table,  file.path(res_dir, "03_performance.csv"), row.names = FALSE)
write.csv(comp_table,  file.path(res_dir, "03_complementarity.csv"), row.names = FALSE)

## ============================================================
## 8. Summary output
## ============================================================
cat("\n========== Module 1 summary ==========\n")
cat("Pooled training set:", length(y), "subjects\n")
cat("  BritoIL:", sum(study_vec == "BritoIL_2016"),
    "| HMP:", sum(study_vec == "HMP_2012"), "\n")
cat("Age range:", round(range(y), 1), "\n\n")

cat("--- Pooled 5x5 CV ---\n")
cat("Oral clock MAE: ", round(res_oral$mae, 2),  " years\n")
cat("Gut clock MAE:  ", round(res_gut$mae, 2),   " years\n")
cat("Joint clock MAE:", round(res_joint$mae, 2), " years\n\n")

cat("--- Complementarity ---\n")
cat("dMAE (joint - oral):", round(comp_table$value[3], 2), "years",
    sprintf("(95%% CI: %.2f to %.2f)\n", ci_o[1], ci_o[2]))
cat("dMAE (joint - gut): ", round(comp_table$value[6], 2), "years",
    sprintf("(95%% CI: %.2f to %.2f)\n", ci_g[1], ci_g[2]))
cat("Residual correlation (Spearman):", round(r_og, 3), "\n")
cat("Partial correlation (gut | oral):", round(pc_gut_given_oral, 3), "\n")
cat("Partial correlation (oral | gut):", round(pc_oral_given_gut, 3), "\n\n")

cat("--- LOSO ---\n")
cat("Oral LOSO MAE: ", round(loso_oral$mae, 2),  " years\n")
cat("Gut LOSO MAE:  ", round(loso_gut$mae, 2),   " years\n")
cat("Joint LOSO MAE:", round(loso_joint$mae, 2), " years\n")

cat("\nmodule1_train_clock DONE\n")
