# ============================================================================
# module5_utils.R — Shared helper functions for module 5
#
# Project: Paired oral-gut microbiome aging clock (OG-Clock)
#
# Provides:
#   get_mat_meta()        list(mat, meta) / bare-matrix compatibility
#   clr_samples()         CLR with within-sample multiplicative replacement
#   load_study_site()     Load one niche from <study>__<type>.rds
#   load_paired_view()    Load a niche from both studies, feature intersection
#   limma_age_assoc()     Per-niche age association (limma, base-R fallback)
#   theme_cell()          Cell-journal style theme
#   reverse_validation()  Shared reverse-validation logic
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  if (requireNamespace("limma", quietly = TRUE)) library(limma)
})

# ---- list(mat, meta) or bare matrix ----
get_mat_meta <- function(x) {
  if (is.list(x) && !is.null(x$mat)) {
    list(mat = x$mat, meta = x$meta)
  } else {
    list(mat = x, meta = data.frame(sample_id = colnames(x)))
  }
}

# ---- CLR (samples x features, within-sample multiplicative replacement) ----
clr_samples <- function(m, pseudocount_factor = 0.5) {
  # m: samples x features
  # Within each sample, zeros are replaced with pseudocount_factor x the
  # smallest positive value, then renormalized to sum = 1 before CLR
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

# ---- Load one study, one niche ----
# Returns list(mat = features x samples, meta = data.frame)
load_study_site <- function(study, typ, site) {
  x <- readRDS(file.path(data_dir, paste0(study, "__", typ, ".rds")))
  x <- get_mat_meta(x)
  
  m  <- x$mat
  pd <- x$meta
  
  # meta row order may not match matrix column order in list(mat, meta)
  # files; body_site filtering aligns positionally, so reorder meta by
  # sample_id first
  if (!"sample_id" %in% colnames(pd)) pd$sample_id <- colnames(m)
  if (!identical(colnames(m), pd$sample_id)) {
    j <- match(colnames(m), pd$sample_id)
    if (any(is.na(j))) {
      stop("load_study_site: meta for ", study, " is missing sample_id for ",
           sum(is.na(j)), " matrix columns; check rds integrity")
    }
    pd <- pd[j, , drop = FALSE]
    cat("  [", study, "] meta row order differed from matrix columns; ",
        "reordered by sample_id
", sep = "")
  }
  
  # Bare matrix (meta has sample_id only): pull metadata from the cMD3
  # sample_metadata.rds by sample_id
  if (!"body_site" %in% colnames(pd)) {
    md3 <- readRDS(file.path(data_dir, "sample_metadata.rds"))
    i <- match(colnames(m), md3$sample_id)
    if (any(is.na(i))) {
      cat("  [", study, "] ", sum(is.na(i)),
          " samples not in sample_metadata, dropped
", sep = "")
      keep_c <- !is.na(i)
      m <- m[, keep_c, drop = FALSE]
      i <- i[keep_c]
    }
    pd <- md3[i, , drop = FALSE]
  }
  
  if (!"sample_id" %in% colnames(pd)) pd$sample_id <- rownames(pd)
  
  # Filter to the requested niche with age available (defensive: missing
  # age/sex would put NA into the limma design matrix)
  keep <- !is.na(pd$body_site) & pd$body_site == site &
    !is.na(pd$age) & pd$age > 0
  m  <- m[, keep, drop = FALSE]
  pd <- pd[keep, , drop = FALSE]
  
  if (nrow(pd) == 0) {
    stop("load_study_site: no samples left after filtering ",
         study, " / ", typ, " / ", site, "; check body_site values")
  }
  
  # Keep one sample per subject
  if ("subject_id" %in% colnames(pd)) {
    pd$.sid <- pd$sample_id
    pd <- pd %>% group_by(subject_id) %>% slice_head(n = 1) %>% ungroup()
    m  <- m[, pd$.sid, drop = FALSE]
  }
  
  # Prefix subject ids with the study name
  prefix <- paste0(study, "_")
  if ("subject_id" %in% colnames(pd)) {
    pd$subject_id <- as.character(pd$subject_id)
    need <- !startsWith(pd$subject_id, prefix)
    pd$subject_id[need] <- paste0(prefix, pd$subject_id[need])
  }
  
  # Keep only analysis-relevant columns
  keep_cols <- c("sample_id", "subject_id", "body_site", "age",
                 "study_name", "sex", "gender", "body_subsite")
  pd <- pd[, intersect(keep_cols, colnames(pd)), drop = FALSE]
  
  if (!"sex" %in% colnames(pd) && "gender" %in% colnames(pd)) {
    pd$sex <- pd$gender
  }
  if (!"study_name" %in% colnames(pd)) pd$study_name <- study
  
  list(mat = m, meta = pd)
}

# ---- Load a niche from both studies, keep the feature intersection ----
# Returns list(X = samples x features, y, meta, feat)
load_paired_view <- function(studies, typ, site,
                             prev_thresh = 0.2,
                             pseudocount_factor = 0.5) {
  mats <- list()
  infos <- list()
  
  for (st in studies) {
    d <- load_study_site(st, typ, site)
    mats[[st]] <- d$mat
    infos[[st]] <- d$meta
  }
  
  # Feature intersection
  common_feat <- sort(Reduce(intersect, lapply(mats, rownames)))
  mats <- lapply(mats, function(m) m[common_feat, , drop = FALSE])
  
  # Studies may disagree on meta columns (e.g. sex vs gender, optional
  # body_subsite), which breaks rbind. Take the column union, fill missing
  # columns with NA, and harmonize column order before rbind.
  req_cols <- c("sample_id", "age")
  for (st in names(infos)) {
    pd <- infos[[st]]
    miss_req <- setdiff(c(req_cols, "subject_id"), colnames(pd))
    if ("subject_id" %in% miss_req) {
      stop("load_paired_view: meta for study ", st,
           " lacks subject_id, required for subject-level alignment; ",
           "columns present: ", paste(colnames(pd), collapse = ", "))
    }
    miss_req2 <- setdiff(req_cols, colnames(pd))
    if (length(miss_req2) > 0) {
      stop("load_paired_view: meta for study ", st,
           " lacks required columns: ", paste(miss_req2, collapse = ", "),
           "; columns present: ", paste(colnames(pd), collapse = ", "))
    }
    if (!"sex" %in% colnames(pd) && !"gender" %in% colnames(pd)) {
      stop("load_paired_view: meta for study ", st, " lacks sex/gender")
    }
  }
  
  all_cols <- Reduce(union, lapply(infos, colnames))
  infos <- lapply(infos, function(pd) {
    miss <- setdiff(all_cols, colnames(pd))
    if (length(miss) > 0) pd[miss] <- NA
    pd$sample_id <- as.character(pd$sample_id)
    if ("subject_id" %in% colnames(pd)) {
      pd$subject_id <- as.character(pd$subject_id)
    }
    pd[, all_cols, drop = FALSE]
  })
  
  # Merge
  X_raw <- do.call(cbind, mats)     # features x samples
  info  <- do.call(rbind, infos)
  
  if (!identical(colnames(X_raw), info$sample_id)) {
    stop("load_paired_view: meta row order does not match matrix column ",
         "order; check that meta row order equals mat column order in ",
         "each study rds")
  }
  
  # Transpose to samples x features
  X <- t(X_raw)
  
  # Prevalence filter
  prev <- colMeans(X > 0)
  keep <- prev >= prev_thresh
  X <- X[, keep, drop = FALSE]
  
  cat(sprintf("  [%s %s] %d samples x %d features (prev >= %.2f)
",
              paste(studies, collapse = "+"), site,
              nrow(X), ncol(X), prev_thresh))
  
  # CLR
  X_clr <- clr_samples(X, pseudocount_factor = pseudocount_factor)
  
  list(
    X    = X_clr,
    y    = info$age,
    meta = info,
    feat = colnames(X_clr)
  )
}

# ---- Age association (limma when available, else per-feature OLS) ----
# X: samples x features
# info: data.frame with age / sex / study_name
# The two implementations were verified equivalent: identical logFC,
# P-value Spearman > 0.999 (difference stems from eBayes variance
# shrinkage only, negligible at n ~= 227)
limma_age_assoc <- function(X, info, coef_name = "age",
                            adjust_study = TRUE) {
  # Defensive: drop samples with missing age/sex (design matrix forbids NA)
  ok <- !is.na(info$age) & !is.na(info$sex)
  if (any(!ok)) {
    cat("  dropping samples with missing age/sex:", sum(!ok), "
")
    X <- X[ok, , drop = FALSE]
    info <- info[ok, , drop = FALSE]
  }
  info$sex <- factor(info$sex)
  
  if (adjust_study && "study_name" %in% colnames(info) &&
      length(unique(info$study_name)) > 1) {
    design <- model.matrix(~ age + sex + study_name, data = info)
  } else {
    design <- model.matrix(~ age + sex, data = info)
  }
  
  has_limma <- requireNamespace("limma", quietly = TRUE)
  
  if (has_limma) {
    # limma expects features x samples
    fit <- limma::lmFit(t(X), design)
    fit <- limma::eBayes(fit)
    tt <- limma::topTable(fit, coef = coef_name, number = Inf,
                          sort.by = "none")
    tt$feature <- rownames(tt)
    tt
  } else {
    # Base-R fallback: per-feature OLS (all features fitted in one lm.fit)
    j <- match(coef_name, colnames(design))
    if (is.na(j)) stop("coefficient ", coef_name, " not in design matrix")
    
    Y <- as.matrix(X)                 # n x features
    fitq  <- lm.fit(design, Y)
    coef_mat <- fitq$coefficients[j, ]
    resid    <- fitq$residuals
    n <- nrow(design); df_res <- n - fitq$rank
    sigma2 <- colSums(resid^2) / df_res
    Rinv   <- chol2inv(qr.R(fitq$qr))         # (X'X)^{-1}
    se     <- sqrt(pmax(sigma2 * Rinv[j, j], 0))
    tstat  <- coef_mat / se
    pval   <- 2 * stats::pt(-abs(tstat), df = df_res)
    
    data.frame(
      logFC = coef_mat, t = tstat,
      P.Value = pval, adj.P.Val = p.adjust(pval, "BH"),
      feature = colnames(X),
      stringsAsFactors = FALSE
    )
  }
}

# ---- Cell-journal style theme ----
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

# ---- Shared reverse-validation logic ----
# src_tt: source-niche topTable
# dst_tt: target-niche topTable
# src_fdr / dst_fdr: significance thresholds
reverse_validation <- function(src_tt, dst_tt,
                               src_fdr = 0.05, dst_fdr = 0.05,
                               feature_col = "feature") {
  src_sig <- src_tt[src_tt$adj.P.Val < src_fdr, feature_col]
  
  m <- src_tt[src_tt[[feature_col]] %in% src_sig, ] %>%
    inner_join(
      dst_tt[, c(feature_col, "logFC", "P.Value", "adj.P.Val")],
      by = feature_col,
      suffix = c("_src", "_dst")
    )
  
  if (nrow(m) == 0) {
    return(list(table = m, n_sig = length(src_sig),
                n_testable = 0, n_same_dir = 0,
                binom_p = NA, spearman_lfc = NA))
  }
  
  m$same_dir <- sign(m$logFC_src) == sign(m$logFC_dst)
  
  bt <- tryCatch(
    binom.test(sum(m$same_dir), nrow(m), p = 0.5),
    error = function(e) list(p.value = NA)
  )
  
  sp <- suppressWarnings(cor(m$logFC_src, m$logFC_dst, method = "spearman"))
  
  list(
    table = m,
    n_sig = length(src_sig),
    n_testable = nrow(m),
    n_same_dir = sum(m$same_dir),
    binom_p = bt$p.value,
    spearman_lfc = sp,
    n_dst_sig = sum(m$adj.P.Val_dst < dst_fdr, na.rm = TRUE)
  )
}
