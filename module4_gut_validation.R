# ============================================================================
# module4_gut_validation.R — External validation of the gut clock in three
# independent shotgun cohorts
#
# Input:
#   <data_dir>/BritoIL_2016__relative_abundance.rds   (bare matrix: species x samples)
#   <data_dir>/HMP_2012__relative_abundance.rds       (bare matrix: species x samples)
#   <data_dir>/sample_metadata.rds                    (cMD3 metadata)
#   <data_dir>/ext_LifeLinesDeep_2016.rds             (list(mat, meta))
#   <data_dir>/ext_AsnicarF_2021.rds                  (list(mat, meta))
#   <data_dir>/ext_KarlssonFH_2013.rds                (list(mat, meta))
#
# Output:
#   <data_dir>/ext_gut_results.rds
#   <res_dir>/external_validation_gut.csv
#   <res_dir>/ext_gut_predictions.csv  (per-sample predictions, used by
#                                       module4_figure5.R panel c)
#   <res_dir>/Figure5a.pdf (external-validation MAE / rho bars)
#
# Design:
#   - Training: Brito + HMP stool species (one sample per subject, same
#     pipeline as module 3)
#   - External: 3 independent shotgun stool cohorts
#   - Feature prefiltering on the training set only (no information leakage)
#   - Missing features imputed with training-set CLR medians
#   - Brito/HMP RDS files are bare matrices (no meta); subject/stool
#     filtering therefore uses sample_metadata.rds (same pipeline as
#     module3_clock_transfer.R)
#   - Fisher-z pooled Spearman uses inverse-variance weighting
#     (pooled_z = sum(w_i * z_i) / sum(w_i), w_i = n_i - 3)
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(openxlsx)
})

set.seed(20250101)

## ============================================================
## 0. Helpers
## ============================================================
# Unified CLR (samples x features, within-sample multiplicative replacement)
clr_samples <- function(m, pseudocount_factor = 0.5) {
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

# Support list(mat, meta) or bare matrices
get_mat_meta <- function(x) {
  if (is.list(x) && !is.null(x$mat)) {
    list(mat = x$mat, meta = x$meta)
  } else {
    list(mat = x, meta = data.frame(sample_id = colnames(x)))
  }
}

# Cell-style theme
theme_cell <- function(base_size = 7, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      text          = element_text(color = "black"),
      plot.title    = element_text(size = base_size + 1, face = "bold",
                                   hjust = 0,
                                   margin = ggplot2::margin(b = 4)),
      axis.title    = element_text(size = base_size),
      axis.text     = element_text(size = base_size - 0.5, color = "black"),
      axis.line     = element_line(color = "black", linewidth = 0.3),
      axis.ticks    = element_line(color = "black", linewidth = 0.3),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text       = element_text(size = base_size, face = "bold"),
      legend.title  = element_text(size = base_size),
      legend.text   = element_text(size = base_size - 0.5),
      legend.key.size = unit(3, "mm"),
      legend.background = element_blank(),
      legend.key    = element_blank(),
      plot.margin   = ggplot2::margin(4, 4, 4, 4)
    )
}

## ============================================================
## 1. Training set (Brito + HMP stool species)
## ============================================================
cat("
========== 1. Training set ==========
")

md3 <- readRDS(file.path(data_dir, "sample_metadata.rds"))
stopifnot(all(c("sample_id", "subject_id", "body_site", "age",
                "study_name") %in% colnames(md3)))

mats_list <- list()
infos_list <- list()

for (st in c("BritoIL_2016", "HMP_2012")) {
  m <- readRDS(file.path(data_dir, paste0(st, "__relative_abundance.rds")))
  if (is.list(m)) m <- m$mat          # support list format
  stopifnot(is.matrix(m))
  
  # Stool samples with age for this study; first sample per subject
  pd <- md3 %>%
    filter(study_name == st,
           body_site == "stool",
           age > 0) %>%
    group_by(subject_id) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  sel <- intersect(pd$sample_id, colnames(m))
  pd  <- pd[match(sel, pd$sample_id), , drop = FALSE]
  m   <- m[, sel, drop = FALSE]
  
  # Study prefix avoids cross-study subject-ID collisions
  pd$subject_id <- paste0(st, "_", pd$subject_id)
  
  mats_list[[st]]  <- m
  infos_list[[st]] <- pd[, c("sample_id", "subject_id", "body_site",
                             "age", "study_name")]
  
  cat(st, ": ", nrow(pd), " samples | ", nrow(m), " species
", sep = "")
}

# Species intersection
common_feat <- sort(Reduce(intersect, lapply(mats_list, rownames)))
cat("Shared species:", length(common_feat), "
")

mats_list <- lapply(mats_list, function(m) m[common_feat, , drop = FALSE])

X_tr_raw <- do.call(cbind, mats_list)          # features x samples
info_tr  <- do.call(rbind, infos_list)
y_tr     <- info_tr$age

cat("Training set:", ncol(X_tr_raw), "samples x", nrow(X_tr_raw), "species
")
cat("Age:", paste(range(y_tr), collapse = "-"), "
")

## ============================================================
## 2. CLR (training set)
## ============================================================
cat("
========== 2. CLR ==========
")

X_tr <- t(X_tr_raw)                             # samples x features
X_tr <- clr_samples(X_tr)

cat("After CLR:", paste(dim(X_tr), collapse = " x "), "
")
cat("CLR range:", paste(round(range(X_tr, na.rm = TRUE), 2),
                        collapse = " - "), "
")

## ============================================================
## 3. Feature prefilter (training set only, no leakage)
## ============================================================
cat("
========== 3. Feature prefilter ==========
")

rho <- suppressWarnings(
  apply(X_tr, 2, function(x) cor(x, y_tr, method = "spearman",
                                 use = "complete.obs"))
)
rho <- rho[!is.na(rho)]
keep <- names(sort(abs(rho), decreasing = TRUE))[
  seq_len(min(150, length(rho)))
]

cat("Features retained:", length(keep), "
")

X_tr_f <- X_tr[, keep, drop = FALSE]
X_tr_f[!is.finite(X_tr_f)] <- 0

# Training-set CLR medians (impute missing external features)
clr_median_train <- apply(X_tr_f, 2, median, na.rm = TRUE)

## ============================================================
## 4. Train ENet
## ============================================================
cat("
========== 4. Train ENet ==========
")

fit <- cv.glmnet(X_tr_f, y_tr, alpha = 0.5, nfolds = 5)

cat("lambda.min:", round(fit$lambda.min, 5), "
")
cat("lambda.1se:", round(fit$lambda.1se, 5), "
")

## ============================================================
## 5. External cohort validation
## ============================================================
cat("
========== 5. External cohort validation ==========
")

metrics <- function(pred, obs) {
  ok <- !is.na(pred) & !is.na(obs)
  pred <- pred[ok]; obs <- obs[ok]
  c(
    n = length(obs),
    MAE = mean(abs(pred - obs)),
    RMSE = sqrt(mean((pred - obs)^2)),
    R2 = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
    Spearman = suppressWarnings(cor(pred, obs, method = "spearman")),
    Pearson  = suppressWarnings(cor(pred, obs, method = "pearson"))
  )
}

ext_studies <- c("LifeLinesDeep_2016", "AsnicarF_2021", "KarlssonFH_2013")

res <- list()
pred_rows <- list()  # collect per-sample predictions

for (st in ext_studies) {
  fp <- file.path(data_dir, paste0("ext_", st, ".rds"))
  if (!file.exists(fp)) {
    warning("Missing external cohort: ", st)
    next
  }
  
  d <- get_mat_meta(readRDS(fp))
  m <- d$mat; pd <- d$meta
  
  # Explicitly align metadata rows with matrix columns
  # (rownames(meta) == colnames(mat))
  if (!is.null(pd$sample_id)) {
    pd <- pd[match(colnames(m), pd$sample_id), , drop = FALSE]
  } else {
    pd <- pd[match(colnames(m), rownames(pd)), , drop = FALSE]
  }
  stopifnot(!any(is.na(pd$age)))
  
  # Samples with age (adult gut cohorts)
  keep_s <- !is.na(pd$age) & pd$age > 0
  m  <- m[, keep_s, drop = FALSE]
  pd <- pd[keep_s, , drop = FALSE]
  
  n_ext <- ncol(m)
  cat(st, ": ", n_ext, " samples | ", nrow(m), " species
", sep = "")
  
  # Feature presence
  feat_present <- intersect(keep, rownames(m))
  feat_missing <- setdiff(keep, rownames(m))
  cat("  Features present: ", length(feat_present), "/", length(keep),
      " (", round(100 * length(feat_present) / length(keep), 1), "%)
",
      sep = "")
  
  # CLR on present features only (note: CLR reference means differ from
  # the training set)
  m_present <- m[feat_present, , drop = FALSE]
  m_present[is.na(m_present)] <- 0
  
  X_ext_present <- t(m_present)                   # samples x features
  X_ext_present <- clr_samples(X_ext_present)
  
  # Assemble the full feature matrix
  X_ext <- matrix(NA_real_, nrow = nrow(X_ext_present),
                  ncol = length(keep),
                  dimnames = list(rownames(X_ext_present), keep))
  X_ext[, feat_present] <- X_ext_present
  
  # Impute missing features with training-set CLR medians
  if (length(feat_missing) > 0) {
    X_ext[, feat_missing] <- rep(clr_median_train[feat_missing],
                                 each = nrow(X_ext))
  }
  X_ext[!is.finite(X_ext)] <- 0
  
  # Predict
  pred <- as.numeric(predict(fit, X_ext, s = "lambda.min"))
  obs  <- as.numeric(pd$age)
  
  m_ext <- metrics(pred, obs)
  res[[st]] <- data.frame(
    cohort       = st,
    n            = m_ext["n"],
    n_features   = length(feat_present),
    feat_present = round(100 * length(feat_present) / length(keep), 1),
    MAE          = m_ext["MAE"],
    RMSE         = m_ext["RMSE"],
    R2           = m_ext["R2"],
    Spearman     = m_ext["Spearman"],
    Pearson      = m_ext["Pearson"],
    median_age   = median(obs, na.rm = TRUE),
    age_range    = paste(round(range(obs, na.rm = TRUE), 1), collapse = "-"),
    stringsAsFactors = FALSE
  )
  
  saveRDS(list(pred = pred, obs = obs, sample = colnames(m)),
          file.path(data_dir, paste0("ext_pred_", st, ".rds")))
  
  pred_rows[[st]] <- data.frame(cohort = st, sample_id = colnames(m),
                                age = obs, pred = pred)
  
  cat(sprintf("  MAE=%.2f rho=%.3f
", m_ext["MAE"], m_ext["Spearman"]))
}

# Fail fast when no external cohort is available (instead of crashing at
# the plotting stage)
missing_cohorts <- setdiff(ext_studies, names(res))
if (length(missing_cohorts) > 0) {
  warning("Missing external cohorts: ", paste(missing_cohorts, collapse = ", "))
}
if (length(res) == 0) {
  stop("No external cohort available. The following files are missing:
",
       paste0("  ", file.path(data_dir, paste0("ext_", ext_studies, ".rds")),
              collapse = "
"),
       "
These files are generated from curatedMetagenomicData; ",
       "run module4_prep_gut.R first.")
}

out <- do.call(rbind, res)

## ============================================================
## 6. Fisher-z pooled Spearman (inverse-variance weighting)
## ============================================================
cat("
========== 6. Pooling ==========
")

# Explicit column types (avoid character coercion after rbind)
for (cc in c("n", "n_features", "feat_present",
             "MAE", "RMSE", "R2", "Spearman", "Pearson",
             "median_age")) {
  out[[cc]] <- as.numeric(out[[cc]])
}

# Fisher-z inverse-variance weighted pooling (approximate for Spearman):
# pooled_z = sum(w_i * z_i) / sum(w_i), w_i = n_i - 3
z_weights <- out$n - 3
pooled_z  <- sum(z_weights * atanh(out$Spearman)) / sum(z_weights)
pooled_rho <- tanh(pooled_z)

# Weighted mean MAE
pooled_mae <- weighted.mean(out$MAE, out$n, na.rm = TRUE)

out_row <- data.frame(
  cohort       = "POOLED",
  n            = sum(out$n),
  n_features   = NA,
  feat_present = NA,
  MAE          = pooled_mae,
  RMSE         = NA,
  R2           = NA,
  Spearman     = pooled_rho,
  Pearson      = NA,
  median_age   = NA,
  age_range    = NA
)

out_full <- rbind(out, out_row)

cat("External validation results:
")
print(out_full, row.names = FALSE, digits = 3)

## ============================================================
## 7. Save
## ============================================================
saveRDS(out_full, file.path(data_dir, "ext_gut_results.rds"))

write.csv(out_full, file.path(res_dir, "external_validation_gut.csv"),
          row.names = FALSE)

# Per-sample prediction table (source data for Figure 5c)
pred_df <- do.call(rbind, pred_rows)
write.csv(pred_df, file.path(res_dir, "ext_gut_predictions.csv"),
          row.names = FALSE)

## ============================================================
## 8. Figure 5a — external-validation MAE / rho bars
## ============================================================
cat("
========== Figure 5a ==========
")

plot_df <- out_full %>%
  filter(cohort != "POOLED") %>%
  select(cohort, MAE, Spearman) %>%
  pivot_longer(cols = c(MAE, Spearman),
               names_to = "metric",
               values_to = "value")

plot_df$cohort <- factor(plot_df$cohort,
                         levels = c("LifeLinesDeep_2016",
                                    "AsnicarF_2021",
                                    "KarlssonFH_2013"))

pal <- c("LifeLinesDeep_2016" = "#1F77B4",
         "AsnicarF_2021"      = "#FF7F0E",
         "KarlssonFH_2013"    = "#2CA02C")

p <- ggplot(plot_df, aes(cohort, value, fill = cohort)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", value)),
            vjust = -0.4, size = 2.2) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL,
       y = NULL,
       title = "External validation of the gut species clock") +
  theme_cell() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(res_dir, "Figure5a.pdf"), p,
       width = 6.5, height = 3.2, units = "in", device = cairo_pdf)

## ============================================================
## 8b. Figure 5b + A/B composite
##     Data come from module4_oral_validation.R (<data_dir>/ext_oral.rds);
##     run order: module4_oral_validation.R first, then this script.
##     If ext_oral.rds is missing, 5b/composite are skipped.
## ============================================================
cat("
========== Figure 5b ==========
")

oral_fp <- file.path(data_dir, "ext_oral.rds")
if (file.exists(oral_fp)) {
  o <- readRDS(oral_fp)
  out_oral <- o$perf
  
  plot_df_b <- out_oral %>%
    select(clock, MAE, Spearman) %>%
    pivot_longer(cols = c(MAE, Spearman),
                 names_to = "metric",
                 values_to = "value")
  
  plot_df_b$clock <- factor(plot_df_b$clock,
                            levels = c("OMAA_RF", "ENet_ours"))
  
  pal_b <- c("ENet_ours" = "#1F77B4",
             "OMAA_RF"   = "#FF7F0E")
  
  p5b <- ggplot(plot_df_b, aes(clock, value, fill = clock)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.2f", value)),
              vjust = -0.4, size = 2.5) +
    facet_wrap(~metric, scales = "free_y") +
    scale_fill_manual(values = pal_b, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = NULL,
         title = "External validation on Huang 2020 oral cohort") +
    theme_cell()
  
  ggsave(file.path(res_dir, "Figure5b.pdf"), p5b,
         width = 4.5, height = 3.2, units = "in", device = cairo_pdf)
  
  # A/B composite (A: gut three-cohort external validation |
  # B: oral external validation vs OMAA)
  fig5 <- p / p5b +
    plot_annotation(
      tag_levels = "A",
      theme = theme(plot.tag = element_text(size = 10, face = "bold"))
    )
  
  ggsave(file.path(res_dir, "Figure5.pdf"), fig5,
         width = 6.85, height = 6.8, units = "in", device = cairo_pdf)
  
  cat("  Figure5b.pdf and Figure5.pdf generated
")
} else {
  cat("  ext_oral.rds not found — skipping Figure5b / Figure5 composite.
")
  cat("  Run module4_oral_validation.R first (creates <data_dir>/ext_oral.rds),
")
  cat("  then re-run this script to complete Figure5b.pdf and Figure5.pdf.
")
}

## ============================================================
## 9. Summary output
## ============================================================
cat("
========== module4_gut_validation DONE ==========
")
cat("Outputs:
")
cat("  ", file.path(data_dir, "ext_gut_results.rds"), "
")
cat("  ", file.path(res_dir, "external_validation_gut.csv"), "
")
cat("  ", file.path(res_dir, "ext_gut_predictions.csv"), "
")
cat("  ", file.path(res_dir, "Figure5a.pdf"), "
")
if (file.exists(file.path(res_dir, "Figure5b.pdf"))) {
  cat("  ", file.path(res_dir, "Figure5b.pdf"), "
")
  cat("  ", file.path(res_dir, "Figure5.pdf"), "
")
}
