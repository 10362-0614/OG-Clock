# ============================================================================
# module4_oral_validation.R — External validation of the oral clock vs OMAA-RF
#
# Project: Paired oral-gut microbiome aging clock (OG-Clock)
#
# Input:
#   <data_dir>/nhanes_analysis.rds
#   <data_dir>/nhanes_omp/dada2rb-taxonomy-annotate.txt
#   <data_dir>/oral_2550.tsv        (Huang 2020 ASV table)
#   <data_dir>/oral_2550_map.txt
#   <data_dir>/oral_taxonomy.txt
#
# Output:
#   <data_dir>/ext_oral.rds
#   <res_dir>/external_validation_oral.csv
#   <res_dir>/Figure5b.pdf
#
# Design:
#   - Training: NHANES oral (2009-2010 + 2011-2012, full set)
#   - Two clocks: ENet (ours) vs OMAA-RF (reimplementation)
#   - External: Huang 2020 oral_2550 (ASVs aggregated to genus level)
#   - Missing genera imputed with training-set CLR medians
#
# Notes:
#   - randomForest defines its own margin generic that masks ggplot2::margin;
#     all margin calls are therefore explicitly qualified as ggplot2::margin
#   - Samples with zero total abundance at the genus level are dropped
#     before CLR (0/0 would give NaN; on the ENet side these were silently
#     zeroed by the !is.finite guard, on the RF side they errored out)
#   - GAM p-values are extracted by grepping the row names of s.table
#     instead of relying on exact row names
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(mgcv)
  library(glmnet)
  library(randomForest)
  library(ggplot2)
  library(tidyr)
})

set.seed(2026)

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
## 1. NHANES training data (genus level)
## ============================================================
cat("
========== 1. NHANES training data ==========
")

r <- readRDS(file.path(data_dir, "nhanes_analysis.rds"))
gen <- r$genus

# Align meta rows to genus columns explicitly (verified identical upstream;
# kept as a defensive check)
stopifnot(identical(rownames(gen), as.character(r$meta$SEQN)))

# Taxonomy annotations
tax <- read.delim(file.path(data_dir, "nhanes_omp",
                            "dada2rb-taxonomy-annotate.txt"),
                  check.names = FALSE)

cat("tax columns:", paste(colnames(tax), collapse = ", "), "
")

# Extract genus label (handle NA)
get_last_taxon <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  parts <- trimws(strsplit(as.character(x), ";")[[1]])
  if (length(parts) == 0) return(NA_character_)
  last <- parts[length(parts)]
  if (is.na(last) || last %in% c("NA", "", "uncultured",
                                 "unidentified", "unclassified",
                                 "unknown")) return(NA_character_)
  last
}

tax$genus <- sapply(tax[["Taxonomy in SILVA v123"]], get_last_taxon)
tax_map   <- setNames(tax$genus, tax$Name)

cat("tax_map coverage:", sum(!is.na(tax_map)), "/", length(tax_map), "
")

# Map feature names to genus labels
new_names <- tax_map[colnames(gen)]
cat("mapped features:", sum(!is.na(new_names)), "/",
    length(new_names), "
")

# Keep features with a valid genus label
keep_col <- !is.na(new_names) & nzchar(new_names)
gen_mapped <- gen[, keep_col, drop = FALSE]
colnames(gen_mapped) <- new_names[keep_col]

# Aggregate duplicate genus labels by summing (not dropping)
cat("features before aggregation:", ncol(gen_mapped),
    "| unique genera:", length(unique(colnames(gen_mapped))), "
")

gen_t <- t(gen_mapped)                        # features x samples
if (anyNA(gen_t)) gen_t[is.na(gen_t)] <- 0

agg <- rowsum(gen_t, group = rownames(gen_t), reorder = FALSE)
X_nh_all <- t(agg)                            # samples x genera

cat("NHANES genus-level:", paste(dim(X_nh_all), collapse = " x "), "
")

# Drop samples with zero total abundance at the genus level:
# CLR would give 0/0 = NaN (see Notes)
rs_nh <- rowSums(X_nh_all)
if (any(rs_nh <= 0)) {
  cat("dropping all-zero samples:", sum(rs_nh <= 0), "
")
  X_nh_all <- X_nh_all[rs_nh > 0, , drop = FALSE]
}

# CLR
X_nh_all <- clr_samples(X_nh_all)

# Align meta to samples explicitly (defensive, order-independent)
meta_nh <- r$meta[match(rownames(X_nh_all),
                        as.character(r$meta$SEQN)), , drop = FALSE]
y_nh <- meta_nh$RIDAGEYR

keep_y <- !is.na(y_nh) & y_nh > 0
X_nh_all <- X_nh_all[keep_y, , drop = FALSE]
y_nh <- y_nh[keep_y]
gender_nh <- factor(meta_nh$RIAGENDR[keep_y])

cat("NHANES training:", nrow(X_nh_all), "samples |",
    "age:", paste(range(y_nh), collapse = "-"), "
")

## ============================================================
## 2. Train the two clocks (ENet + OMAA-RF)
## ============================================================
cat("
========== 2. Training the two clocks ==========
")

# ---- 2a. ENet (top-150 by |Spearman|) ----
rho <- suppressWarnings(
  apply(X_nh_all, 2, function(x) cor(x, y_nh, method = "spearman",
                                     use = "complete.obs"))
)
rho <- rho[!is.na(rho)]
sel <- names(sort(abs(rho), decreasing = TRUE))[
  seq_len(min(150, length(rho)))
]

X_nh_enet <- X_nh_all[, sel, drop = FALSE]
X_nh_enet[!is.finite(X_nh_enet)] <- 0

fit_enet <- cv.glmnet(X_nh_enet, y_nh, alpha = 0.5, nfolds = 5)
clr_median_enet <- apply(X_nh_enet, 2, median, na.rm = TRUE)

cat("ENet training done | features:", length(sel), "
")

# ---- 2b. OMAA-RF (GAM screening + random forest) ----
cat("GAM screening...
")
df_base <- data.frame(
  age    = y_nh,
  gender = gender_nh
)

gam_res <- lapply(colnames(X_nh_all), function(g) {
  df_base$g <- X_nh_all[, g]
  fit <- tryCatch(
    gam(g ~ s(age) + gender, data = df_base, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(data.frame(Genus = g, p = NA_real_))
  }
  st <- summary(fit)$s.table
  ri <- grep("s\\(age\\)", rownames(st))[1]
  data.frame(
    Genus = g,
    p     = if (is.na(ri)) NA_real_ else st[ri, "p-value"]
  )
})
gam_res <- do.call(rbind, gam_res)

omaa_genera <- gam_res %>%
  filter(!is.na(p), p.adjust(p, method = "BH") < 0.05) %>%
  pull(Genus)

cat("OMAA genera retained:", length(omaa_genera), "
")

rf_full <- randomForest(
  x = X_nh_all[, omaa_genera, drop = FALSE],
  y = y_nh,
  ntree = 500,
  mtry  = max(1, floor(sqrt(length(omaa_genera))))
)

clr_median_omaa <- apply(X_nh_all[, omaa_genera, drop = FALSE],
                         2, median, na.rm = TRUE)

## ============================================================
## 3. Huang 2020 external cohort (ASVs -> genus level)
## ============================================================
cat("
========== 3. Huang 2020 ==========
")

asv <- read.delim(file.path(data_dir, "oral_2550.tsv"),
                  check.names = FALSE, row.names = 1)
cat("ASV table:", paste(dim(asv), collapse = " x "), "
")

taxh <- read.delim(file.path(data_dir, "oral_taxonomy.txt"),
                   check.names = FALSE)
taxh_map <- setNames(taxh$Taxon, taxh$`Feature ID`)

# Extract genus label from the SILVA-style lineage
hgenus_of <- function(seqid) {
  tx <- unname(taxh_map[seqid])
  if (is.na(tx) || !nzchar(tx)) return(NA_character_)
  parts <- trimws(strsplit(tx, ";")[[1]])
  g <- grep("^g__", parts, value = TRUE)
  if (length(g) == 0) return(NA_character_)
  gsub("^g__", "", g[1])
}

asv_genera <- sapply(colnames(asv), hgenus_of)
okh <- !is.na(asv_genera) & nzchar(asv_genera) &
  !grepl("unidentified|uncultured|unclassified|unknown",
         asv_genera, ignore.case = TRUE)

cat("ASVs with genus:", sum(okh), "/", ncol(asv), "
")

# Aggregate ASVs to genus (samples x ASVs -> samples x genera)
asv_ok <- asv[, okh, drop = FALSE]
gnames <- asv_genera[okh]

# rowsum expects features x samples
agg_h <- rowsum(t(asv_ok), group = gnames, reorder = FALSE)
X_h_all <- t(agg_h)                            # samples x genera

cat("Huang genus-level:", paste(dim(X_h_all), collapse = " x "), "
")

# Ages
map_h <- read.delim(file.path(data_dir, "oral_2550_map.txt"),
                    check.names = FALSE, row.names = 1)
age_h <- suppressWarnings(
  as.numeric(as.character(map_h[rownames(X_h_all), "age"]))
)
keep_h <- !is.na(age_h) & age_h > 0
X_h_all <- X_h_all[keep_h, , drop = FALSE]
age_h <- age_h[keep_h]

cat("Huang with age:", length(age_h),
    "| age range:", paste(range(age_h), collapse = "-"), "
")

# Drop all-zero samples here as well (CLR would give NaN)
rs_h <- rowSums(X_h_all)
if (any(rs_h <= 0)) {
  cat("dropping all-zero samples:", sum(rs_h <= 0), "
")
  X_h_all <- X_h_all[rs_h > 0, , drop = FALSE]
  age_h   <- age_h[rs_h > 0]
}

# CLR (across all genera)
X_h_all <- clr_samples(X_h_all)

## ============================================================
## 4. Apply the two clocks
## ============================================================
cat("
========== 4. Applying the two clocks ==========
")

# ---- 4a. ENet ----
common_enet <- intersect(sel, colnames(X_h_all))
miss_enet <- setdiff(sel, common_enet)

cat("ENet features present:", length(common_enet), "/", length(sel),
    " (", round(100 * length(common_enet) / length(sel), 1), "%)
", sep = "")

X_h_enet <- matrix(NA_real_, nrow = nrow(X_h_all), ncol = length(sel),
                   dimnames = list(rownames(X_h_all), sel))
X_h_enet[, common_enet] <- X_h_all[, common_enet, drop = FALSE]
if (length(miss_enet) > 0) {
  # Impute missing genera with training-set CLR medians
  X_h_enet[, miss_enet] <- rep(clr_median_enet[miss_enet],
                               each = nrow(X_h_enet))
}
X_h_enet[!is.finite(X_h_enet)] <- 0

pred_enet <- as.numeric(predict(fit_enet, X_h_enet, s = "lambda.min"))

# ---- 4b. OMAA-RF ----
common_omaa <- intersect(omaa_genera, colnames(X_h_all))
miss_omaa <- setdiff(omaa_genera, common_omaa)

cat("OMAA features present:", length(common_omaa), "/", length(omaa_genera),
    " (", round(100 * length(common_omaa) / length(omaa_genera), 1), "%)
",
    sep = "")

X_h_omaa <- matrix(NA_real_, nrow = nrow(X_h_all), ncol = length(omaa_genera),
                   dimnames = list(rownames(X_h_all), omaa_genera))
X_h_omaa[, common_omaa] <- X_h_all[, common_omaa, drop = FALSE]
if (length(miss_omaa) > 0) {
  X_h_omaa[, miss_omaa] <- rep(clr_median_omaa[miss_omaa],
                               each = nrow(X_h_omaa))
}
X_h_omaa[!is.finite(X_h_omaa)] <- 0

pred_omaa <- as.numeric(predict(rf_full, X_h_omaa))

## ============================================================
## 5. Performance evaluation
## ============================================================
cat("
========== 5. Performance evaluation ==========
")

met <- function(p, o) {
  ok <- !is.na(p) & !is.na(o)
  p <- p[ok]; o <- o[ok]
  c(
    n        = length(o),
    MAE      = mean(abs(p - o)),
    RMSE     = sqrt(mean((p - o)^2)),
    R2       = 1 - sum((o - p)^2) / sum((o - mean(o))^2),
    Spearman = suppressWarnings(cor(p, o, method = "spearman")),
    Pearson  = suppressWarnings(cor(p, o, method = "pearson"))
  )
}

# Build result rows with explicitly ordered columns
make_row <- function(name, m) {
  data.frame(
    clock    = name,
    n        = as.numeric(m["n"]),
    MAE      = as.numeric(m["MAE"]),
    RMSE     = as.numeric(m["RMSE"]),
    R2       = as.numeric(m["R2"]),
    Spearman = as.numeric(m["Spearman"]),
    Pearson  = as.numeric(m["Pearson"]),
    stringsAsFactors = FALSE
  )
}

m_enet <- met(pred_enet, age_h)
m_omaa <- met(pred_omaa, age_h)

out <- rbind(
  make_row("ENet_ours", m_enet),
  make_row("OMAA_RF",   m_omaa)
)

cat("Huang external validation results:
")
print(out, row.names = FALSE, digits = 3)

## ============================================================
## 6. Save
## ============================================================
saveRDS(
  list(
    pred_enet = pred_enet,
    pred_omaa = pred_omaa,
    age       = age_h,
    perf      = out,
    n_enet_feat = length(common_enet),
    n_omaa_feat = length(common_omaa),
    n_enet_total = length(sel),
    n_omaa_total = length(omaa_genera)
  ),
  file.path(data_dir, "ext_oral.rds")
)

write.csv(out, file.path(res_dir, "external_validation_oral.csv"),
          row.names = FALSE)

## ============================================================
## 7. Figure 5b
## ============================================================
cat("
========== Figure 5b ==========
")

plot_df <- out %>%
  select(clock, MAE, Spearman) %>%
  pivot_longer(cols = c(MAE, Spearman),
               names_to = "metric",
               values_to = "value")

plot_df$clock <- factor(plot_df$clock,
                        levels = c("OMAA_RF", "ENet_ours"))

pal <- c("ENet_ours" = "#1F77B4",
         "OMAA_RF"   = "#FF7F0E")

p <- ggplot(plot_df, aes(clock, value, fill = clock)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", value)),
            vjust = -0.4, size = 2.5) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = NULL,
       title = "External validation on Huang 2020 oral cohort") +
  theme_cell()

ggsave(file.path(res_dir, "Figure5b.pdf"), p,
       width = 4.5, height = 3.2, units = "in", device = cairo_pdf)

## ============================================================
## 8. Summary
## ============================================================
cat("
========== module4_oral_validation DONE ==========
")
cat("Outputs:
")
cat("  ", file.path(data_dir, "ext_oral.rds"), "
")
cat("  ", file.path(res_dir, "external_validation_oral.csv"), "
")
cat("  ", file.path(res_dir, "Figure5b.pdf"), "
")
