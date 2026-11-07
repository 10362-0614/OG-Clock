# ============================================================================
# module5_figure6.R — Figure 6 (panels a-e)
#
# Project: Paired oral-gut microbiome aging clock (OG-Clock), module 5
#
# Panels:
#   a  Marginal contribution bars: MAE + subject-level bootstrap 95% CI
#      + delta-MAE annotation
#   b  Incremental R2 and partial correlations
#   c  Niche pathway overlap: species-stratified (intersection 0) vs
#      de-stratified (shared 83)
#   d  De-stratified shared-pathway logFC concordance scatter
#      (same/discordant direction + Spearman)
#   e  Top-15 age-associated pathways per niche (pathway-id labels only,
#      UNINTEGRATED excluded)
#
# Input:
#   <data_dir>/marginal_contribution.rds
#   <data_dir>/pillarC_pathways.rds
#   De-stratified per-pathway table (panel d; probed in order, must
#   contain logFC_oral / logFC_gut):
#     <res_dir>/destratified_shared_pathways.csv
#     <res_dir>/TableS11_reverse_validation.csv
#
# Output:
#   <res_dir>/Figure6a.pdf ... Figure6e.pdf
#   <res_dir>/Figure6_v2.pdf (layout (a|b)/(c|d)/e
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

## ============================================================
## 0. Theme and palettes (same Cell style as module5_utils.R)
## ============================================================
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

pal_clock <- c("Oral" = "#4C72B0", "Gut" = "#55A868", "Joint" = "#C44E52")
pal_niche <- c("oral" = "#4C72B0", "gut"  = "#55A868")

## ============================================================
## 1. Load artifacts
## ============================================================
cat("========== Loading ==========
")

mc_fp <- file.path(data_dir, "marginal_contribution.rds")
pw_fp <- file.path(data_dir, "pillarC_pathways.rds")
stopifnot(file.exists(mc_fp), file.exists(pw_fp))

mc <- readRDS(mc_fp)
pw <- readRDS(pw_fp)

cat("marginal_contribution: n_subj =", mc$n_subj, "
")
cat("pillarC: oral sig =", pw$n_oral_sig,
    "| gut sig =", pw$n_gut_sig,
    "| stratified shared =", length(pw$shared_sig), "
")

## ============================================================
## Panel a — MAE bars + bootstrap CI + delta-MAE annotation
## ============================================================
cat("
========== Panel a ==========
")

y <- mc$y
n <- length(y)

# Subject-level bootstrap MAE 95% CI per clock (same B = 2000 as the
# primary analysis)
set.seed(20250101)
mae_ci <- sapply(mc$preds, function(p) {
  b <- replicate(2000, {
    idx <- sample.int(n, n, replace = TRUE)
    mean(abs(p[idx] - y[idx]), na.rm = TRUE)
  })
  quantile(b, c(0.025, 0.975), na.rm = TRUE)
})

plot_a <- data.frame(
  clock = factor(c("Oral", "Gut", "Joint"), levels = c("Oral", "Gut", "Joint")),
  MAE   = mc$res$MAE,
  lo    = mae_ci[1, ],
  hi    = mae_ci[2, ]
)

p_a <- ggplot(plot_a, aes(clock, MAE, fill = clock)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.3) +
  # Value labels above the error bars so they do not overlap the whiskers
  geom_text(aes(y = hi, label = sprintf("%.2f", MAE)),
            vjust = -0.4, size = 2.5) +
  scale_fill_manual(values = pal_clock, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL,
       y = "MAE (years), 5\u00d75 repeated CV",
       title = "Marginal contribution of each niche clock",
       subtitle = sprintf(
         "\u0394MAE joint \u2212 oral: %+.2f [%+.2f, %+.2f];  joint \u2212 gut: %+.2f [%+.2f, %+.2f]",
         mc$delta_joint_oral, mc$ci_joint_oral[1], mc$ci_joint_oral[2],
         mc$delta_joint_gut,  mc$ci_joint_gut[1],  mc$ci_joint_gut[2])) +
  theme_cell()

ggsave(file.path(res_dir, "Figure6a.pdf"), p_a,
       width = 4.2, height = 3.6, units = "in", device = cairo_pdf)

## ============================================================
## Panel b — incremental R2 and partial correlations
## ============================================================
cat("========== Panel b ==========
")

plot_b <- data.frame(
  metric = factor(c("Incremental R\u00b2
oral \u2192 joint",
                    "Incremental R\u00b2
gut \u2192 joint",
                    "Partial r
oral | gut",
                    "Partial r
gut | oral"),
                  levels = rev(c("Incremental R\u00b2
oral \u2192 joint",
                                 "Incremental R\u00b2
gut \u2192 joint",
                                 "Partial r
oral | gut",
                                 "Partial r
gut | oral"))),
  value = c(mc$r2["joint"] - mc$r2["oral"],
            mc$r2["joint"] - mc$r2["gut"],
            mc$pc_oral_given_gut,
            mc$pc_gut_given_oral),
  group = rep(c("Incremental R\u00b2", "Partial correlation"), each = 2)
)

p_b <- ggplot(plot_b, aes(metric, value, fill = group)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%+.3f", value)),
            hjust = ifelse(plot_b$value >= 0, -0.15, 1.15), size = 2.4) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Incremental R\u00b2" = "#4C72B0",
                               "Partial correlation" = "#CCB974"),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
  coord_flip() +
  labs(x = NULL, y = "Unique contribution of each niche",
       title = "Incremental R\u00b2 and partial correlation") +
  theme_cell()

ggsave(file.path(res_dir, "Figure6b.pdf"), p_b,
       width = 4.2, height = 3.6, units = "in", device = cairo_pdf)

## ============================================================
## Panel c — pathway overlap: stratified (intersection 0) vs
## de-stratified (shared 83)
## ============================================================
cat("========== Panel c ==========
")

# De-stratified counts (from the de-stratified analysis / Table S10:
# oral 128, gut 134, shared 83, EXCL_UNINTEGRATED = TRUE convention).
# Update here if the final convention changes.
n_ds_oral <- 128; n_ds_gut <- 134; n_ds_shared <- 83

plot_c <- data.frame(
  level = rep(c("Species-stratified", "De-stratified"), each = 3),
  cat = factor(rep(c("Oral only", "Shared", "Gut only"), 2),
               levels = c("Oral only", "Shared", "Gut only")),
  n = c(pw$n_oral_sig, length(pw$shared_sig), pw$n_gut_sig,
        n_ds_oral - n_ds_shared, n_ds_shared, n_ds_gut - n_ds_shared)
)

plot_c$level <- factor(plot_c$level,
                       levels = c("Species-stratified", "De-stratified"))

p_c <- ggplot(plot_c, aes(level, n, fill = cat)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.2) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.8), vjust = -0.35, size = 2.4) +
  scale_fill_manual(values = c("Oral only" = "#4C72B0",
                               "Shared"    = "#C44E52",
                               "Gut only"  = "#55A868"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Age-associated pathways (FDR < 0.05)",
       title = "Niche overlap at pathway level") +
  theme_cell() +
  theme(legend.position = "top")

ggsave(file.path(res_dir, "Figure6c.pdf"), p_c,
       width = 4.2, height = 3.6, units = "in", device = cairo_pdf)

## ============================================================
## Panel d — de-stratified shared-pathway logFC concordance scatter
## ============================================================
cat("========== Panel d ==========
")

cand <- file.path(res_dir, "destratified_shared_pathways.csv")
cand2 <- file.path(res_dir, "TableS11_reverse_validation.csv")
ds_fp <- c(cand, cand2)[file.exists(c(cand, cand2))][1]

p_d <- NULL
if (!is.na(ds_fp)) {
  ds <- read.csv(ds_fp, check.names = FALSE)
  need <- c("logFC_oral", "logFC_gut")
  if (all(need %in% colnames(ds))) {
    ds$same_dir <- sign(ds$logFC_oral) == sign(ds$logFC_gut)
    rho_ds <- suppressWarnings(cor(ds$logFC_oral, ds$logFC_gut,
                                   method = "spearman"))
    binom_p <- tryCatch(binom.test(sum(ds$same_dir), nrow(ds), p = 0.5)$p.value,
                        error = function(e) NA_real_)
    cat("panel d:", nrow(ds), "pathways | same direction", sum(ds$same_dir),
        "| Spearman", round(rho_ds, 3), "
")
    
    p_d <- ggplot(ds, aes(logFC_oral, logFC_gut, color = same_dir)) +
      geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey40") +
      geom_point(size = 1.1, alpha = 0.75) +
      scale_color_manual(values = c("TRUE" = "#4C72B0",
                                    "FALSE" = "#D55E00"),
                         labels = c("TRUE" = "Same direction",
                                    "FALSE" = "Discordant"),
                         name = NULL) +
      labs(x = "log FC per year, oral",
           y = "log FC per year, gut",
           title = "Directional concordance of shared pathways",
           subtitle = sprintf("n = %d shared pathways; Spearman \u03c1 = %.2f; same direction %d/%d (binomial P = %.1e)",
                              nrow(ds), rho_ds, sum(ds$same_dir), nrow(ds),
                              binom_p)) +
      theme_cell() +
      theme(legend.position = "top")
  } else {
    cat("  ", basename(ds_fp), "lacks logFC_oral/logFC_gut columns; ",
        "skipping panel d
")
  }
} else {
  cat("  No de-stratified per-pathway table found; skipping panel d.
")
  cat("  Run module5_destratified_export.R to generate it.
")
}

if (!is.null(p_d)) {
  ggsave(file.path(res_dir, "Figure6d.pdf"), p_d,
         width = 4.2, height = 3.6, units = "in", device = cairo_pdf)
}

## ============================================================
## Panel e — top-15 age-associated pathways per niche
## (pathway-id labels only; UNINTEGRATED excluded, matching the
## revised Table S11 convention)
## ============================================================
cat("========== Panel e ==========
")

# Labels keep the pathway id only (full names in Table S11); the facet
# strip already shows the niche. When the same pathway id enters the
# top 15 for multiple organisms, an abbreviated organism tag is appended.
mk_top <- function(tt, niche, k = 15) {
  d <- tt %>%
    filter(adj.P.Val < 0.05,
           !grepl("^UNINTEGRATED", feature)) %>%   # drop non-pathway features (Table S11 convention)
    arrange(P.Value) %>%
    slice_head(n = k)
  # Organism abbreviation: ...|g__Neisseria.s__Neisseria_cinerea -> N. cinerea
  d$org_s <- vapply(d$feature, function(f) {
    org <- sub(".*[|]", "", f)
    s <- sub(".*s__", "", org)
    s <- gsub("_", " ", s)
    if (grepl(" ", s)) {
      w <- strsplit(s, " ")[[1]]
      paste0(substr(w[1], 1, 1), ".", w[2])
    } else s
  }, character(1), USE.NAMES = FALSE)
  d$pw_id  <- trimws(sub(":.*$", "", d$feature))   # "PWY-6071: ..." -> "PWY-6071"
  d$label  <- d$pw_id
  dup <- d$pw_id %in% d$pw_id[duplicated(d$pw_id)]
  d$label[dup] <- paste0(d$pw_id[dup], " [", d$org_s[dup], "]")
  d$niche  <- niche
  d
}

top_oe <- mk_top(pw$oral, "oral")
top_ge <- mk_top(pw$gut,  "gut")

# Sort by logFC within niche (largest decrease on top); the global factor
# ordering is fixed here, and facet free_y drops the other niche's levels
# within each facet automatically
plot_e <- bind_rows(top_oe, top_ge) %>%
  group_by(niche) %>%
  arrange(logFC, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(label = factor(label, levels = unique(label)))

p_e <- ggplot(plot_e, aes(logFC, label, color = niche)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(size = 1.8, shape = 16) +
  facet_wrap(~niche, scales = "free_y") +
  scale_color_manual(values = pal_niche, guide = "none") +
  labs(x = expression("Age-associated log FC per year (natural-log CLR units)"),
       y = NULL,
       title = "Top 15 age-associated organism-stratified pathways per niche") +
  theme_cell() +
  theme(axis.text.y = element_text(size = 6))

ggsave(file.path(res_dir, "Figure6e.pdf"), p_e,
       width = 9, height = 5.2, units = "in", device = cairo_pdf)

## ============================================================
## Assemble — (a|b) / (c|d) / e (degrades to (a|b) / c / e without d)
## ============================================================
cat("
========== Assembling ==========
")

top_row <- if (!is.null(p_d)) (p_a | p_b) / (p_c | p_d) else (p_a | p_b) / p_c

fig6_v2 <- top_row / p_e +
  plot_layout(heights = if (!is.null(p_d)) c(1, 1, 1.5) else c(1, 1.4, 1.5)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 10, face = "bold"))
  )

ggsave(file.path(res_dir, "Figure6_v2.pdf"), fig6_v2,
       width = 9, height = 12, units = "in", device = cairo_pdf)

## ============================================================
## Summary
## ============================================================
cat("
========== module5_figure6 DONE ==========
")
cat("Outputs:
")
for (f in c("Figure6a.pdf", "Figure6b.pdf", "Figure6c.pdf",
            "Figure6d.pdf", "Figure6e.pdf", "Figure6_v2.pdf")) {
  fp <- file.path(res_dir, f)
  if (file.exists(fp)) cat("  ", fp, "
")
}
