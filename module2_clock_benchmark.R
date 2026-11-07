# ============================================================================
# module2_clock_benchmark.R — ENet vs OMAA-RF reproduction and extrapolation
#
# Input:  <data_dir>/nhanes_analysis.rds
# Output: <data_dir>/nhanes_clocks.rds = list(dat, X_all, acc_enet, acc_omaa,
#                                            omaa_genera, perf, preds, ...)
#
# Design:
#   - Training: NHANES 2009-2010 (discovery)
#   - External validation: NHANES 2011-2012
#   - Method comparison: our ENet vs OMAA-RF (GAM screen + random forest)
#   - OMAA genera screened on d09 only; d11 used as an independent test set
#   - CLR with multiplicative replacement (consistent with module 1)
#   - Age acceleration defined uniformly as acc = pred - age
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(mgcv)
  library(glmnet)
  library(randomForest)
})

set.seed(2026)

## ============================================================
## 1. Load preprocessed data
## ============================================================
cat("\n========== 1. Load data ==========\n")

r <- readRDS(file.path(data_dir, "nhanes_analysis.rds"))
meta  <- r$meta
genus <- r$genus

cat("Samples:", nrow(meta), " | genus features:", ncol(genus), "\n")

# Ensure genus rownames align with meta$SEQN
stopifnot(identical(rownames(genus), as.character(meta$SEQN)))

## ============================================================
## 2. CLR transformation (multiplicative replacement, as in module 1)
## ============================================================
cat("\n========== 2. CLR transformation ==========\n")

# Input/output: samples x features matrix; zeros replaced by
# pseudocount_factor x the sample's minimum positive value, rows
# re-normalized, then CLR
clr_mat <- function(m, pseudocount_factor = 0.5) {
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
  
  # Re-normalize (row sums = 1)
  row_sums <- rowSums(m2)
  m2 <- sweep(m2, 1, row_sums, "/")
  
  # CLR: log(x) - mean(log(x))
  log_m <- log(m2)
  log_m - rowMeans(log_m)
}

X_all <- clr_mat(as.matrix(genus))
colnames(X_all) <- colnames(genus)

cat("CLR dimensions:", paste(dim(X_all), collapse = " x "), "\n")
cat("CLR range:", paste(round(range(X_all, na.rm = TRUE), 2), collapse = " - "), "\n")

## ============================================================
## 3. Train/test split (cycle read directly from XPT)
## ============================================================
cat("\n========== 3. Train/test split ==========\n")

if (!"cycle" %in% colnames(meta)) {
  stop("meta lacks the 'cycle' column; re-run module2_preprocess.R")
}

idx09 <- which(meta$cycle == "d09")
idx11 <- which(meta$cycle == "d11")

cat("d09 (training):", length(idx09), "samples\n")
cat("d11 (test):", length(idx11), "samples\n")

y <- meta$RIDAGEYR
cat("Age range:", paste(range(y), collapse = " - "), "\n")

# Factorize covariates
meta$race   <- factor(meta$RIDRETH1)
meta$gender <- factor(meta$RIAGENDR)

## ============================================================
## 4. OMAA-RF reproduction (GAM screen + random forest)
## ============================================================
cat("\n========== 4. OMAA-RF reproduction ==========\n")

# GAM screen: d09 only
gam_screen <- function(idx) {
  cat("  GAM screen:", length(idx), "samples...\n")
  
  df_base <- data.frame(
    age    = y[idx],
    gender = meta$gender[idx],
    race   = meta$race[idx]
  )
  
  res <- lapply(colnames(X_all), function(g) {
    df_base$g <- X_all[idx, g]
    fit <- tryCatch(
      gam(g ~ s(age) + gender + race, data = df_base, method = "REML"),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(data.frame(Genus = g, p = NA_real_, edf = NA_real_))
    }
    s <- summary(fit)$s.table
    data.frame(
      Genus = g,
      p     = s["s(age)", "p-value"],
      edf   = s["s(age)", "edf"]
    )
  })
  
  do.call(rbind, res)
}

t0 <- Sys.time()
g09 <- gam_screen(idx09)
g09$FDR <- p.adjust(g09$p, method = "BH")

cat("  GAM screen time:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
cat("  d09 genera with FDR < 0.05:", sum(g09$FDR < 0.05, na.rm = TRUE), "\n")

# Key design point: screen on d09 only; d11 is never used for screening
omaa_genera <- g09 %>%
  filter(!is.na(FDR), FDR < 0.05) %>%
  arrange(p) %>%
  pull(Genus)

cat("  OMAA genera retained:", length(omaa_genera), "\n")

# RF parameters (as in the original publication): ntree=500, mtry=sqrt(p)
rf_fit <- randomForest(
  x = X_all[idx09, omaa_genera, drop = FALSE],
  y = y[idx09],
  ntree = 500,
  mtry  = max(1, floor(sqrt(length(omaa_genera)))),
  importance = TRUE
)

# In-sample d09 predictions (diagnostic only, not for evaluation)
omaa_pred09_in <- predict(rf_fit, X_all[idx09, omaa_genera, drop = FALSE])
# External d11 predictions
omaa_pred11 <- predict(rf_fit, X_all[idx11, omaa_genera, drop = FALSE])

# 10-fold CV within d09 (for evaluating d09 performance)
set.seed(2026)
folds09 <- sample(rep(1:10, length.out = length(idx09)))
omaa_pred09_cv <- rep(NA_real_, length(idx09))
for (f in 1:10) {
  tr <- which(folds09 != f)
  te <- which(folds09 == f)
  fit_f <- randomForest(
    x = X_all[idx09[tr], omaa_genera, drop = FALSE],
    y = y[idx09[tr]],
    ntree = 500,
    mtry  = max(1, floor(sqrt(length(omaa_genera))))
  )
  omaa_pred09_cv[te] <- predict(fit_f, X_all[idx09[te], omaa_genera, drop = FALSE])
}

cat("  OMAA d09 10-fold CV done\n")

## ============================================================
## 5. Our ENet (in-fold Spearman prefilter + cv.glmnet)
## ============================================================
cat("\n========== 5. Our ENet ==========\n")

enet_pred <- function(idx_tr, idx_te,
                      nfeat = 150, alpha = 0.5,
                      lambda_choice = "lambda.min",
                      seed = 2026) {
  set.seed(seed)
  
  Xtr <- X_all[idx_tr, , drop = FALSE]
  ytr <- y[idx_tr]
  Xte <- X_all[idx_te, , drop = FALSE]
  
  # In-fold Spearman prefilter
  rho <- suppressWarnings(
    apply(Xtr, 2, function(x) cor(x, ytr, method = "spearman",
                                  use = "complete.obs"))
  )
  rho <- rho[!is.na(rho)]
  
  n_keep <- min(nfeat, length(rho))
  if (n_keep < 5) {
    stop("Fewer than 5 features after prefiltering; check idx_tr or data.")
  }
  keep <- names(sort(abs(rho), decreasing = TRUE))[seq_len(n_keep)]
  
  Xtr_f <- Xtr[, keep, drop = FALSE]
  Xte_f <- Xte[, keep, drop = FALSE]
  
  # Guard against NAs (none expected after CLR, but be safe)
  if (anyNA(Xtr_f)) Xtr_f[is.na(Xtr_f)] <- 0
  if (anyNA(Xte_f)) Xte_f[is.na(Xte_f)] <- 0
  
  fit <- cv.glmnet(
    x = Xtr_f, y = ytr,
    family = "gaussian",
    alpha = alpha,
    nfolds = 5,
    standardize = TRUE
  )
  
  pred <- as.numeric(predict(fit, Xte_f, s = fit[[lambda_choice]]))
  list(pred = pred, fit = fit, features = keep)
}

# 10-fold CV within d09 (for evaluating d09 performance)
set.seed(2026)
folds09_enet <- sample(rep(1:10, length.out = length(idx09)))
enet_pred09_cv <- rep(NA_real_, length(idx09))
for (f in 1:10) {
  tr <- which(folds09_enet != f)
  te <- which(folds09_enet == f)
  res <- enet_pred(idx09[tr], idx09[te], seed = 2026 + f)
  enet_pred09_cv[te] <- res$pred
}

cat("  ENet d09 10-fold CV done\n")

# Train on full d09, predict d11
enet_res_09_11 <- enet_pred(idx09, idx11, seed = 2026)
enet_pred11 <- enet_res_09_11$pred

cat("  ENet d11 external prediction done\n")

## ============================================================
## 6. Performance evaluation
## ============================================================
cat("\n========== 6. Performance evaluation ==========\n")

met <- function(pred, obs) {
  ok <- !is.na(pred) & !is.na(obs)
  pred <- pred[ok]; obs <- obs[ok]
  
  MAE  <- mean(abs(pred - obs))
  RMSE <- sqrt(mean((pred - obs)^2))
  # R2: 1 - SSE/SST around the observed mean
  R2   <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
  rho  <- suppressWarnings(cor(pred, obs, method = "spearman"))
  r    <- suppressWarnings(cor(pred, obs, method = "pearson"))
  
  c(n = length(obs), MAE = MAE, RMSE = RMSE, R2 = R2,
    Spearman = rho, Pearson = r)
}

# Pooled 10-fold CV over the full population (native reference; same
# evaluation set as the transfer experiments)
set.seed(2026)
folds_all <- sample(rep(1:10, length.out = nrow(meta)))
enet_pred_all_cv <- rep(NA_real_, nrow(meta))
for (f in 1:10) {
  tr <- which(folds_all != f)
  te <- which(folds_all == f)
  res_f <- enet_pred(tr, te, seed = 2026 + f)
  enet_pred_all_cv[te] <- res_f$pred
}
m_enet_all <- met(enet_pred_all_cv, y)
cat(sprintf("  ENet pooled CV: MAE=%.2f rho=%.3f\n",
            m_enet_all["MAE"], m_enet_all["Spearman"]))

# Within-d09 CV
m_omaa_09 <- met(omaa_pred09_cv, y[idx09])
m_enet_09 <- met(enet_pred09_cv, y[idx09])

# d11 external validation
m_omaa_11 <- met(omaa_pred11, y[idx11])
m_enet_11 <- met(enet_pred11, y[idx11])

cat(sprintf("d09 (10-fold CV): OMAA-RF MAE=%.2f rho=%.3f | ENet MAE=%.2f rho=%.3f\n",
            m_omaa_09["MAE"], m_omaa_09["Spearman"],
            m_enet_09["MAE"], m_enet_09["Spearman"]))

cat(sprintf("d11 (external):   OMAA-RF MAE=%.2f rho=%.3f | ENet MAE=%.2f rho=%.3f\n",
            m_omaa_11["MAE"], m_omaa_11["Spearman"],
            m_enet_11["MAE"], m_enet_11["Spearman"]))

## ============================================================
## 7. Age acceleration (uniformly pred - age)
## ============================================================
cat("\n========== 7. Age acceleration ==========\n")

acc_omaa <- rep(NA_real_, nrow(meta))
acc_omaa[idx09] <- omaa_pred09_cv - y[idx09]
acc_omaa[idx11] <- omaa_pred11 - y[idx11]

acc_enet <- rep(NA_real_, nrow(meta))
acc_enet[idx09] <- enet_pred09_cv - y[idx09]
acc_enet[idx11] <- enet_pred11 - y[idx11]

cat("ENet acceleration SD (years):", round(sd(acc_enet, na.rm = TRUE), 2), "\n")
cat("OMAA acceleration SD (years):", round(sd(acc_omaa, na.rm = TRUE), 2), "\n")

# By cycle
cat("\nAcceleration SD by cycle:\n")
cat("  ENet d09:", round(sd(acc_enet[idx09], na.rm = TRUE), 2),
    "| d11:", round(sd(acc_enet[idx11], na.rm = TRUE), 2), "\n")
cat("  OMAA d09:", round(sd(acc_omaa[idx09], na.rm = TRUE), 2),
    "| d11:", round(sd(acc_omaa[idx11], na.rm = TRUE), 2), "\n")

## ============================================================
## 8. Save
## ============================================================
perf <- data.frame(
  method = c("OMAA-RF", "ENet", "ENet"),
  cohort = c("discovery_CV", "discovery_CV", "validation_ext",
             "pooled_CV")[c(1, 2, 4)],
  MAE    = c(m_omaa_09["MAE"], m_enet_09["MAE"], m_enet_all["MAE"]),
  Spearman = c(m_omaa_09["Spearman"], m_enet_09["Spearman"],
               m_enet_all["Spearman"]),
  Pearson  = c(m_omaa_09["Pearson"], m_enet_09["Pearson"],
               m_enet_all["Pearson"]),
  R2       = c(m_omaa_09["R2"], m_enet_09["R2"], m_enet_all["R2"]),
  stringsAsFactors = FALSE
)
# Append the two d11 external-validation rows (ENet / OMAA)
perf <- rbind(perf, data.frame(
  method = c("ENet", "OMAA-RF"), cohort = c("validation_ext", "validation_ext"),
  MAE = c(m_enet_11["MAE"], m_omaa_11["MAE"]),
  Spearman = c(m_enet_11["Spearman"], m_omaa_11["Spearman"]),
  Pearson = c(m_enet_11["Pearson"], m_omaa_11["Pearson"]),
  R2 = c(m_enet_11["R2"], m_omaa_11["R2"]),
  stringsAsFactors = FALSE
))

out <- list(
  dat        = meta,
  X_all      = X_all,
  acc_omaa   = acc_omaa,
  acc_enet   = acc_enet,
  omaa_genera = omaa_genera,
  n_omaa      = length(omaa_genera),
  perf        = perf,
  preds       = list(
    omaa_d09_cv = omaa_pred09_cv,
    omaa_d11    = omaa_pred11,
    enet_d09_cv = enet_pred09_cv,
    enet_d11    = enet_pred11,
    enet_all_cv = enet_pred_all_cv,
    y_d09       = y[idx09],
    y_d11       = y[idx11]
  ),
  idx09 = idx09,
  idx11 = idx11
)

saveRDS(out, file.path(data_dir, "nhanes_clocks.rds"))

cat("\n========== module2_clock_benchmark DONE ==========\n")
cat("Output:", file.path(data_dir, "nhanes_clocks.rds"), "\n")

cat("\nPerformance summary:\n")
print(perf, row.names = FALSE, digits = 3)
