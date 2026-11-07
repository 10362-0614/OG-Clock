# ============================================================================
# module3_icc.R — Longitudinal ICC of clock residuals (ANOVA method-of-
# moments estimator, subject-level bootstrap CIs)
#
# Input:
#   <data_dir>/hmp2_clock_resid.rds
#   <data_dir>/hmp2_species_matrix.rds
#
# Output:
#   <data_dir>/hmp2_icc_results.rds
#   <res_dir>/hmp2_icc_results.csv
#
# Notes:
#   - Subjects with a single sample contribute NA to per-subject SD and
#     are filtered explicitly
#   - Bootstrap seeds controlled explicitly (B = 500)
#   - Stratified ICCs skipped when < 10 subjects
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(vegan)
})

set.seed(2026)

## ============================================================
## 1. Load
## ============================================================
cat("
========== 1. Load data ==========
")

r <- readRDS(file.path(data_dir, "hmp2_clock_resid.rds"))
md <- r$md
names(md)[names(md) == "Participant ID"] <- "subject_id"

cat("Samples:", nrow(md), " | subjects:", length(unique(md$subject_id)), "
")

# Shannon diversity
hmp2 <- readRDS(file.path(data_dir, "hmp2_species_matrix.rds"))
sh <- diversity(t(hmp2), index = "shannon")
md$shannon <- sh[md$sid]

# Predicted age (OOF)
md$pred_oof <- md$oof_resid + md$age_est

## ============================================================
## 2. ICC computation (ANOVA method-of-moments)
## ============================================================
icc_point <- function(x, g) {
  d <- data.frame(x = x, g = g)
  d <- d[complete.cases(d), ]
  N <- nrow(d)
  k <- table(d$g)
  m <- length(k)
  
  if (m < 2 || N <= m) return(list(icc = NA, N = N, m = m))
  
  gm  <- mean(d$x)
  msb <- sum(tapply(d$x, d$g, function(v) length(v) * (mean(v) - gm)^2)) / (m - 1)
  msw <- sum(tapply(d$x, d$g, function(v) sum((v - mean(v))^2))) / (N - m)
  
  if (msw < 1e-10) return(list(icc = NA, N = N, m = m))
  
  k0  <- (N - sum(k^2) / N) / (m - 1)
  icc <- (msb - msw) / (msb + (k0 - 1) * msw)
  
  list(icc = icc, N = N, m = m)
}

icc_one <- function(d, var, nboot = 500, seed = 2026) {
  set.seed(seed)
  
  p <- icc_point(d[[var]], d$subject_id)
  subs <- unique(d$subject_id)
  
  if (length(subs) < 10) {
    return(c(ICC = p$icc, low = NA, high = NA,
             n = p$N, nsubj = p$m))
  }
  
  bs <- replicate(nboot, {
    ss <- sample(subs, replace = TRUE)
    db <- do.call(rbind, lapply(seq_along(ss), function(i) {
      x <- d[d$subject_id == ss[i], ]
      x$bid <- i
      x
    }))
    icc_point(db[[var]], db$bid)$icc
  })
  
  # as.numeric strips quantile names ("2.5%"/"97.5%"); otherwise
  # c(low = ci[1]) concatenates names into "low.2.5%" and downstream
  # x["low"] returns NA (root cause of all-NA CIs)
  ci <- tryCatch(
    as.numeric(quantile(bs, c(0.025, 0.975), na.rm = TRUE)),
    error = function(e) c(NA_real_, NA_real_)
  )
  
  c(ICC = p$icc, low = ci[1], high = ci[2], n = p$N, nsubj = p$m)
}

## ============================================================
## 3. Primary ICC + references
## ============================================================
cat("
========== 3. Primary ICC ==========
")

rows <- list()

rows[["Clock residual (HMP2-native, OOF)"]] <- icc_one(md, "oof_resid")
rows[["Clock residual (cMD3-transferred)"]] <- icc_one(md, "resid")
rows[["Predicted age (OOF)"]]              <- icc_one(md, "pred_oof")
rows[["Shannon diversity"]]                <- icc_one(md, "shannon")

# Negative control
set.seed(2026)
md$rand <- sample(md$oof_resid)
rows[["Permuted residual (control)"]]      <- icc_one(md, "rand")

## ============================================================
## 4. Stratified ICC (diagnostic)
## ============================================================
cat("
========== 4. Stratified ICC ==========
")

for (dx in c("CD", "UC", "nonIBD")) {
  sub <- md[md$diagnosis == dx, ]
  n_subj <- length(unique(sub$subject_id))
  cat(sprintf("%s: n=%d, subjects=%d
", dx, nrow(sub), n_subj))
  
  if (n_subj >= 10) {
    rows[[sprintf("Clock residual - %s", dx)]] <- icc_one(sub, "oof_resid")
  } else {
    cat("  Too small (< 10 subjects), skipped
")
    rows[[sprintf("Clock residual - %s", dx)]] <-
      c(ICC = NA, low = NA, high = NA, n = nrow(sub), nsubj = n_subj)
  }
}

## ============================================================
## 5. Summary table
## ============================================================
icc_tab <- bind_rows(lapply(names(rows), function(nm) {
  x <- rows[[nm]]
  data.frame(
    comparison = nm,
    ICC  = unname(x["ICC"]),
    low  = unname(x["low"]),
    high = unname(x["high"]),
    n    = unname(x["n"]),
    n_subjects = unname(x["nsubj"]),
    stringsAsFactors = FALSE
  )
}))

cat("
ICC summary:
")
print(icc_tab, row.names = FALSE, digits = 3)

## ============================================================
## 6. Per-subject residual SD (single-sample subjects filtered)
## ============================================================
aux <- md %>%
  group_by(subject_id) %>%
  summarise(
    n = n(),
    sd_resid = if (n() > 1) sd(oof_resid, na.rm = TRUE) else NA_real_,
    span_days = max(interval_days, na.rm = TRUE),
    .groups = "drop"
  )

valid_aux <- aux %>% filter(n > 1, !is.na(sd_resid))

cat(sprintf("
Per-subject residual SD: median=%.2f IQR=%.2f-%.2f
",
            median(valid_aux$sd_resid, na.rm = TRUE),
            quantile(valid_aux$sd_resid, 0.25, na.rm = TRUE),
            quantile(valid_aux$sd_resid, 0.75, na.rm = TRUE)))
cat(sprintf("Follow-up span (days): median=%d
",
            as.integer(median(aux$span_days, na.rm = TRUE))))

## ============================================================
## 7. Save
## ============================================================
out <- list(
  icc_tab = icc_tab,
  aux     = aux,
  md      = md,
  n_samples   = nrow(md),
  n_subjects  = length(unique(md$subject_id)),
  native_perf = r$native_perf,
  transfer_rho = r$transfer_rho
)

saveRDS(out, file.path(data_dir, "hmp2_icc_results.rds"))

write.csv(icc_tab, file.path(res_dir, "hmp2_icc_results.csv"),
          row.names = FALSE)

cat("
========== module3_icc DONE ==========
")
cat("Outputs:
")
cat("  ", file.path(data_dir, "hmp2_icc_results.rds"), "
")
cat("  ", file.path(res_dir, "hmp2_icc_results.csv"), "
")
