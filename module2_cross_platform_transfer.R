# ============================================================================
# module2_cross_platform_transfer.R — Shotgun-to-16S cross-platform transfer
# (cMD3 oral shotgun clock -> NHANES 16S genus-level data)
#
# Input:
#   <data_dir>/BritoIL_2016__relative_abundance.rds
#   <data_dir>/HMP_2012__relative_abundance.rds
#   <data_dir>/sample_metadata.rds
#   <data_dir>/nhanes_omp/dada2rb-taxonomy-annotate.txt
#   <data_dir>/nhanes_omp/dada2rb-genus-relative.txt
#   <data_dir>/nhanes_analysis.rds
#
# Output:
#   <data_dir>/nhanes_transfer.rds = list(pred_raw, pred_z, pred_cb,
#     pred_cal, y, n_shared, shared, perf, ref_native)
#
# Harmonization schemes compared:
#   1. raw CLR projection (per-platform CLR)
#   2. z-score harmonization (post-CLR)
#   3. ComBat (age-protected via mod)
#   4. raw projection + recalibration (NHANES-native upper bound, not a
#      genuine transfer)
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(Matrix)
})

set.seed(2026)

## ============================================================
## 1. cMD3 oral: aggregate species -> genus
## ============================================================
cat("\n========== 1. cMD3 oral species -> genus ==========\n")

sp1 <- readRDS(file.path(data_dir, "BritoIL_2016__relative_abundance.rds"))
sp2 <- readRDS(file.path(data_dir, "HMP_2012__relative_abundance.rds"))

# Support both list(mat, meta) and bare matrices
get_mat <- function(x) if (is.list(x)) x$mat else x
sp1_mat <- get_mat(sp1)
sp2_mat <- get_mat(sp2)

# Species intersection
common_sp <- intersect(rownames(sp1_mat), rownames(sp2_mat))
cat("Shared species:", length(common_sp), "\n")

sp <- cbind(sp1_mat[common_sp, , drop = FALSE],
            sp2_mat[common_sp, , drop = FALSE])

# Extract genus from species names (MetaPhlAn2 convention:
# k__|p__|...|g__Genus|s__Genus_species)
genus_of_species <- function(n) {
  g <- sub(".*\\|g__([^\\|]+).*", "\\1", n)
  g[!grepl("g__", n)] <- NA
  g
}

gnames <- genus_of_species(common_sp)
cat("Species with genus extracted:", sum(!is.na(gnames)), "/", length(common_sp), "\n")
cat("Unique genera:", length(unique(gnames[!is.na(gnames)])), "\n")

# Drop species without a genus assignment
keep <- !is.na(gnames)
sp      <- sp[keep, , drop = FALSE]
gnames  <- gnames[keep]

# NA handling (rowsum does not accept NA)
if (anyNA(sp)) {
  cat("NA count:", sum(is.na(sp)), "-> filled with 0\n")
  sp[is.na(sp)] <- 0
}

# Aggregate to genus (sum of relative abundances)
agg <- rowsum(sp, group = gnames, reorder = FALSE)

cat("cMD3: species =", nrow(sp), "-> genera =", nrow(agg), "\n")
cat("First 10 genera:", paste(head(rownames(agg), 10), collapse = ", "), "\n")

# Drop unclassified genera
bad_pattern <- "uncultured|unidentified|unclassified|unknown|^NA$|^$"
keep_agg <- !grepl(bad_pattern, rownames(agg), ignore.case = TRUE)
cat("Unclassified genera dropped:", sum(!keep_agg), "\n")
agg <- agg[keep_agg, , drop = FALSE]

cat("Genera retained:", nrow(agg), "\n")

## ============================================================
## 2. cMD3 oral metadata (align sample_id prefixes)
## ============================================================
cat("\n========== 2. cMD3 oral metadata ==========\n")

md <- readRDS(file.path(data_dir, "sample_metadata.rds"))

oral <- md %>%
  filter(body_site == "oralcavity",
         study_name %in% c("BritoIL_2016", "HMP_2012"),
         age > 0) %>%
  group_by(study_name, subject_id) %>%
  slice_head(n = 1) %>%
  ungroup()

# Prefix subject IDs with the study name (idempotent)
oral$subject_id <- as.character(oral$subject_id)
need_prefix <- !startsWith(oral$subject_id, paste0(oral$study_name, "_"))
oral$subject_id[need_prefix] <- paste0(
  oral$study_name[need_prefix], "_", oral$subject_id[need_prefix]
)

cat("cMD3 oral samples:", nrow(oral), "\n")
cat("  Cohort distribution:\n"); print(table(oral$study_name))
cat("  Age range:", paste(range(oral$age), collapse = "-"), "\n")

# Align with agg column names
sel <- intersect(oral$sample_id, colnames(agg))
cat("Samples matching agg:", length(sel), "/", nrow(oral), "\n")

if (length(sel) < 50) {
  stop("Fewer than 50 matched samples; check sample_id alignment")
}

X_cmd3 <- t(agg[, sel, drop = FALSE])   # samples x genera
y_cmd3 <- oral$age[match(sel, oral$sample_id)]

cat("cMD3 training set:", nrow(X_cmd3), "samples x", ncol(X_cmd3), "genera\n")

## ============================================================
## 3. NHANES genus data + taxonomy mapping
## ============================================================
cat("\n========== 3. NHANES genus data ==========\n")

# Taxonomy annotations
tax <- read.delim(file.path(data_dir, "nhanes_omp",
                            "dada2rb-taxonomy-annotate.txt"),
                  check.names = FALSE)

cat("tax columns:", paste(colnames(tax), collapse = ", "), "\n")

# Extract genus name (last field of the SILVA lineage)
get_last_taxon <- function(x) {
  parts <- strsplit(as.character(x), ";")[[1]]
  if (length(parts) == 0) return(NA_character_)
  last <- trimws(parts[length(parts)])
  if (last %in% c("NA", "", "uncultured", "unidentified",
                  "unclassified", "unknown")) return(NA_character_)
  last
}

tax$genus <- sapply(tax[["Taxonomy in SILVA v123"]], get_last_taxon)

# Map RB_genusN -> genus name
tax_map <- setNames(tax$genus, tax$Name)

cat("First 5 tax_map entries:\n")
print(head(tax_map))

# NHANES genus-level table
gen <- read.delim(file.path(data_dir, "nhanes_omp",
                            "dada2rb-genus-relative.txt"),
                  check.names = FALSE)
gen <- gen[!is.na(gen$SEQN) & gen$SEQN != "", ]
gen$SEQN <- as.numeric(as.character(gen$SEQN))

r2 <- readRDS(file.path(data_dir, "nhanes_analysis.rds"))
keep_seqn <- r2$meta$SEQN

gen2 <- gen[gen$SEQN %in% keep_seqn, ]

X_nh_raw <- as.matrix(gen2[, -1])
rownames(X_nh_raw) <- gen2$SEQN

# Map column names
new_names <- tax_map[colnames(X_nh_raw)]
cat("Columns mapped successfully:", sum(!is.na(new_names)), "/",
    length(new_names), "\n")

# Drop columns that failed to map
keep_col <- !is.na(new_names) & nzchar(new_names)
X_nh_raw <- X_nh_raw[, keep_col, drop = FALSE]
colnames(X_nh_raw) <- new_names[keep_col]

# Aggregate duplicate genus names (distinct ASVs may map to the same
# genus); rowsum needs features x samples
X_nh_raw <- t(X_nh_raw)   # now features x samples
if (anyNA(X_nh_raw)) {
  cat("NHANES NA count:", sum(is.na(X_nh_raw)), "-> filled with 0\n")
  X_nh_raw[is.na(X_nh_raw)] <- 0
}
agg_nh <- rowsum(X_nh_raw, group = rownames(X_nh_raw), reorder = FALSE)
X_nh <- t(agg_nh)   # back to samples x genera

cat("NHANES genera:", ncol(X_nh), "\n")
cat("NHANES samples:", nrow(X_nh), "\n")

## ============================================================
## 4. Shared genera
## ============================================================
cat("\n========== 4. Shared genera ==========\n")

# Compare in lower case (genus spelling may differ in case across platforms)
cMD3_genera <- rownames(agg)
NHANES_genera <- colnames(X_nh)

shared <- intersect(tolower(cMD3_genera), tolower(NHANES_genera))
cat("Shared genera (lower case):", length(shared), "\n")
cat("cMD3 genera:", length(cMD3_genera), "\n")
cat("NHANES genera:", length(NHANES_genera), "\n")
cat("Shared fraction:",
    round(100 * length(shared) / length(cMD3_genera), 1),
    "% of cMD3 genera\n")
cat("Shared fraction:",
    round(100 * length(shared) / length(NHANES_genera), 1),
    "% of NHANES genera\n")

cat("\nFirst 20 shared genera:\n")
print(head(shared, 20))

if (length(shared) < 15) {
  warning("Fewer than 15 shared genera; transfer performance may be poor")
}

# Align casing (use the cMD3 spelling)
map_to_cmd3 <- setNames(cMD3_genera, tolower(cMD3_genera))
shared_cMD3   <- map_to_cmd3[shared]
shared_NHANES <- colnames(X_nh)[match(shared, tolower(colnames(X_nh)))]

# Check for duplicates
stopifnot(length(shared_cMD3) == length(unique(shared_cMD3)))
stopifnot(length(shared_NHANES) == length(unique(shared_NHANES)))

## ============================================================
## 5. Cross-platform harmonization
## ============================================================
cat("\n========== 5. Cross-platform harmonization ==========\n")

# ---- CLR transformation (multiplicative replacement, as in modules 1-2) ----
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
  row_sums <- rowSums(m2)
  m2 <- sweep(m2, 1, row_sums, "/")
  log_m <- log(m2)
  log_m - rowMeans(log_m)
}

# Extract shared-genus submatrices
C_raw <- X_cmd3[, shared_cMD3, drop = FALSE]
N_raw <- X_nh[,  shared_NHANES, drop = FALSE]

# ---- Feature prefilter: top-150 |Spearman| within the training set
# (cMD3 only; no NHANES information leakage) ----
rho_screen <- suppressWarnings(apply(C_raw, 2, function(x)
  cor(x, y_cmd3, method = "spearman", use = "complete.obs")))
n_keep_tr <- min(150, sum(!is.na(rho_screen)))
keep_genera <- names(sort(abs(rho_screen), decreasing = TRUE))[seq_len(n_keep_tr)]
cat("Prefilter:", n_keep_tr, "/", ncol(C_raw), "genera enter the transfer experiment\n")
C_raw <- C_raw[, keep_genera, drop = FALSE]
N_raw <- N_raw[, keep_genera, drop = FALSE]
shared_cMD3   <- keep_genera
shared_NHANES <- keep_genera  # column names unified after subsetting

# Unify column names (cMD3 spelling)
colnames(N_raw) <- shared_cMD3

# ---- Scheme 1: raw CLR (per-platform CLR) ----
C_clr <- clr_mat(C_raw)
N_clr <- clr_mat(N_raw)

cat("C_clr:", paste(dim(C_clr), collapse = " x "), "\n")
cat("N_clr:", paste(dim(N_clr), collapse = " x "), "\n")
cat("C_clr range:", paste(round(range(C_clr), 2), collapse = " - "), "\n")
cat("N_clr range:", paste(round(range(N_clr), 2), collapse = " - "), "\n")

# ---- Scheme 3: ComBat harmonization (age protected) ----
# ComBat requires a features x samples matrix
comb <- cbind(t(C_clr), t(N_clr))   # features x samples
batch <- c(rep("cmd3", nrow(C_clr)), rep("nhanes", nrow(N_clr)))

# Protect age
age_vec <- c(y_cmd3, rep(NA, nrow(N_clr)))
# Append NHANES ages
nh_seqn <- rownames(N_raw)
nh_age <- r2$meta$RIDAGEYR[match(nh_seqn, r2$meta$SEQN)]
age_vec[(nrow(C_clr) + 1):length(age_vec)] <- nh_age

cat("ComBat: batch distribution\n"); print(table(batch, useNA = "ifany"))
cat("ComBat: age non-NA:", sum(!is.na(age_vec)), "/", length(age_vec), "\n")

C_cb <- C_clr
N_cb <- N_clr

if (requireNamespace("sva", quietly = TRUE)) {
  # Build mod (age protected): model.matrix(~ age); requires age without NA
  if (anyNA(age_vec)) {
    warning("age contains NA; ComBat mod set to NULL")
    mod <- NULL
  } else {
    mod <- model.matrix(~ age_vec)
  }
  
  cb_res <- tryCatch(
    sva::ComBat(comb, batch, mod = mod, par.prior = TRUE),
    error = function(e) {
      warning("ComBat failed: ", conditionMessage(e))
      NULL
    }
  )
  
  if (!is.null(cb_res)) {
    # Split back into cMD3 / NHANES (transpose to samples x features)
    C_cb <- t(cb_res[, 1:nrow(C_clr)])
    N_cb <- t(cb_res[, (nrow(C_clr) + 1):ncol(cb_res)])
    colnames(C_cb) <- shared_cMD3
    colnames(N_cb) <- shared_cMD3
    cat("ComBat done\n")
  } else {
    warning("ComBat failed; scheme skipped")
  }
} else {
  warning("sva package not installed; ComBat skipped")
}

## ============================================================
## 6. Train cMD3 oral clock, apply to NHANES
## ============================================================
cat("\n========== 6. Transfer ==========\n")

y_nh <- nh_age

# ---- 6a. Raw projection ----
fit_cmd3 <- cv.glmnet(C_clr, y_cmd3, alpha = 0.5, nfolds = 5)
pred_nh_raw <- as.numeric(predict(fit_cmd3, N_clr, s = "lambda.min"))

m_raw <- c(
  n = length(y_nh),
  MAE = mean(abs(pred_nh_raw - y_nh), na.rm = TRUE),
  Spearman = suppressWarnings(cor(pred_nh_raw, y_nh,
                                  method = "spearman",
                                  use = "complete.obs")),
  Pearson = suppressWarnings(cor(pred_nh_raw, y_nh,
                                 method = "pearson",
                                 use = "complete.obs"))
)

cat(sprintf("Raw projection: MAE=%.2f rho=%.3f r=%.3f\n",
            m_raw["MAE"], m_raw["Spearman"], m_raw["Pearson"]))

# ---- 6b. z-score harmonization (applied post-CLR) ----
zscale_clr <- function(m) {
  mu <- colMeans(m, na.rm = TRUE)
  s  <- apply(m, 2, sd, na.rm = TRUE)
  s[s == 0 | is.na(s)] <- 1
  m2 <- sweep(m, 2, mu, "-")
  sweep(m2, 2, s, "/")
}

Cz <- zscale_clr(C_clr)
Nz <- zscale_clr(N_clr)

fit_cmd3z <- cv.glmnet(Cz, y_cmd3, alpha = 0.5, nfolds = 5)
pred_nh_z <- as.numeric(predict(fit_cmd3z, Nz, s = "lambda.min"))

m_z <- c(
  n = length(y_nh),
  MAE = mean(abs(pred_nh_z - y_nh), na.rm = TRUE),
  Spearman = suppressWarnings(cor(pred_nh_z, y_nh,
                                  method = "spearman",
                                  use = "complete.obs")),
  Pearson = suppressWarnings(cor(pred_nh_z, y_nh,
                                 method = "pearson",
                                 use = "complete.obs"))
)

cat(sprintf("z-score (post-CLR): MAE=%.2f rho=%.3f r=%.3f\n",
            m_z["MAE"], m_z["Spearman"], m_z["Pearson"]))

# ---- 6c. ComBat harmonization ----
m_cb <- c(n = length(y_nh), MAE = NA, Spearman = NA, Pearson = NA)
pred_nh_cb <- rep(NA_real_, nrow(N_cb))

if (!is.null(C_cb) && !is.null(N_cb)) {
  fit_cmd3cb <- cv.glmnet(C_cb, y_cmd3, alpha = 0.5, nfolds = 5)
  pred_nh_cb <- as.numeric(predict(fit_cmd3cb, N_cb, s = "lambda.min"))
  
  m_cb <- c(
    n = length(y_nh),
    MAE = mean(abs(pred_nh_cb - y_nh), na.rm = TRUE),
    Spearman = suppressWarnings(cor(pred_nh_cb, y_nh,
                                    method = "spearman",
                                    use = "complete.obs")),
    Pearson = suppressWarnings(cor(pred_nh_cb, y_nh,
                                   method = "pearson",
                                   use = "complete.obs"))
  )
  cat(sprintf("ComBat: MAE=%.2f rho=%.3f r=%.3f\n",
              m_cb["MAE"], m_cb["Spearman"], m_cb["Pearson"]))
}

# ---- 6d. Raw + recalibration ----
# Note: this uses NHANES's own data to recalibrate — an upper-bound
# reference, not a genuine transfer performance
ok_cal <- !is.na(pred_nh_raw) & !is.na(y_nh)
cat("Recalibration complete cases:", sum(ok_cal), "/", length(y_nh),
    "(NA predictions:", sum(is.na(pred_nh_raw)), ")\n")
cal <- lm(y_nh[ok_cal] ~ pred_nh_raw[ok_cal])
pred_cal <- rep(NA_real_, length(y_nh))
pred_cal[ok_cal] <- as.numeric(fitted(cal))

m_cal <- c(
  n = length(y_nh),
  MAE = mean(abs(pred_cal - y_nh), na.rm = TRUE),
  Spearman = suppressWarnings(cor(pred_cal, y_nh,
                                  method = "spearman",
                                  use = "complete.obs")),
  Pearson = suppressWarnings(cor(pred_cal, y_nh,
                                 method = "pearson",
                                 use = "complete.obs"))
)

cat(sprintf("Raw + recalibration: MAE=%.2f rho=%.3f r=%.3f\n",
            m_cal["MAE"], m_cal["Spearman"], m_cal["Pearson"]))

# ---- Reference: NHANES-native clock (pooled CV, from
# nhanes_clocks.rds; same evaluation set as the transfer experiments) ----
clock_rds <- file.path(data_dir, "nhanes_clocks.rds")
if (file.exists(clock_rds)) {
  ck <- readRDS(clock_rds)
  ref_row <- ck$perf[ck$perf$cohort == "pooled_CV", ]
  ref_native <- sprintf("MAE=%.2f, rho=%.3f", ref_row$MAE, ref_row$Spearman)
} else {
  warning("nhanes_clocks.rds not found; native reference missing")
  ref_native <- NA_character_
}
cat("NHANES-native reference:", ref_native, "\n")

## ============================================================
## 7. Save
## ============================================================
perf <- data.frame(
  method = c("raw projection",
             "z-score (post-CLR)",
             "ComBat harmonization",
             "raw + recalibration"),
  MAE       = c(m_raw["MAE"], m_z["MAE"], m_cb["MAE"], m_cal["MAE"]),
  Spearman  = c(m_raw["Spearman"], m_z["Spearman"],
                m_cb["Spearman"], m_cal["Spearman"]),
  Pearson   = c(m_raw["Pearson"], m_z["Pearson"],
                m_cb["Pearson"], m_cal["Pearson"]),
  n_shared  = length(shared),
  stringsAsFactors = FALSE
)

out <- list(
  pred_raw = pred_nh_raw,
  pred_z   = pred_nh_z,
  pred_cb  = pred_nh_cb,
  pred_cal = pred_cal,
  y        = y_nh,
  n_shared = length(shared),
  shared   = shared_cMD3,
  n_cmd3_genera  = nrow(agg),
  n_nhanes_genera = ncol(X_nh),
  perf     = perf,
  ref_native = ref_native
)

saveRDS(out, file.path(data_dir, "nhanes_transfer.rds"))

cat("\n========== module2_cross_platform_transfer DONE ==========\n")
cat("Output:", file.path(data_dir, "nhanes_transfer.rds"), "\n")
cat("Shared genera:", length(shared), "\n")
cat("\nPerformance summary:\n")
print(perf, row.names = FALSE, digits = 3)
