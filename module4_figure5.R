# ============================================================================
# module4_figure5.R — Figure 5 (panels a-d)
#
# Project: Paired oral-gut microbiome aging clock (OG-Clock)
#
# Panels (final layout matches the Figure 5 legend):
#   a  Gut clock external validation, MAE / rho bars (3 cohorts)
#   b  Gut transfer predicted-vs-observed scatter, 1x3 (reads saved
#      per-sample predictions; no model is retrained)
#   c  Oral clock external validation vs OMAA-RF bars
#   d  Transfer gradient overview (11 experiments, 3 groups,
#      reference line at rho = 0.3)
#
# Input:
#   <res_dir>/external_validation_gut.csv       (module4_gut_validation.R)
#   <data_dir>/ext_oral.rds                     (module4_oral_validation.R)
#   <data_dir>/ext_pred_LifeLinesDeep_2016.rds  (module4_gut_validation.R)
#   <data_dir>/ext_pred_AsnicarF_2021.rds
#   <data_dir>/ext_pred_KarlssonFH_2013.rds
#
# Output:
#   <res_dir>/Figure5c.pdf
#   <res_dir>/Figure5d.pdf
#   <res_dir>/Figure5_v2.pdf
#
# Run order: module4_oral_validation.R and module4_gut_validation.R must
# have completed; this script retrains nothing.
#
# Note (panel d): "cMD3 to HMP2 (shotgun)" is shotgun-to-shotgun across
# cohorts, so it is grouped with "Same platform, cross-cohort" rather than
# with the shotgun-to-16S cross-platform rows.
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
## 0. Theme (same Cell style as module4_gut_validation.R)
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

pal3 <- c("LifeLinesDeep_2016" = "#1F77B4",
          "AsnicarF_2021"      = "#FF7F0E",
          "KarlssonFH_2013"    = "#2CA02C")

## ============================================================
## Panel a — gut clock external validation bars (3 cohorts)
## ============================================================
cat("========== Panel a ==========
")

gut_csv <- file.path(res_dir, "external_validation_gut.csv")
stopifnot(file.exists(gut_csv))
out_gut <- read.csv(gut_csv, stringsAsFactors = FALSE)

plot_a <- out_gut %>%
  filter(cohort != "POOLED") %>%
  select(cohort, MAE, Spearman) %>%
  pivot_longer(cols = c(MAE, Spearman),
               names_to = "metric", values_to = "value")

plot_a$cohort <- factor(plot_a$cohort,
                        levels = c("LifeLinesDeep_2016",
                                   "AsnicarF_2021",
                                   "KarlssonFH_2013"))

p_a <- ggplot(plot_a, aes(cohort, value, fill = cohort)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf(ifelse(metric == "Spearman",
                                       "%.3f", "%.2f"), value)),
            vjust = -0.4, size = 2.2) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = pal3, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = NULL,
       title = "Gut clock: external validation (three cohorts)") +
  theme_cell() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

## ============================================================
## Panel b — oral clock external validation vs OMAA-RF bars
## ============================================================
cat("========== Panel b ==========
")

oral_fp <- file.path(data_dir, "ext_oral.rds")
stopifnot(file.exists(oral_fp))   # if missing, run module4_oral_validation.R first
out_oral <- readRDS(oral_fp)$perf

plot_b <- out_oral %>%
  select(clock, MAE, Spearman) %>%
  pivot_longer(cols = c(MAE, Spearman),
               names_to = "metric", values_to = "value")

plot_b$clock <- factor(plot_b$clock, levels = c("OMAA_RF", "ENet_ours"))

p_b <- ggplot(plot_b, aes(clock, value, fill = clock)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", value)),
            vjust = -0.4, size = 2.5) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = c("ENet_ours" = "#1F77B4",
                               "OMAA_RF"   = "#FF7F0E"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_x_discrete(labels = c("OMAA_RF" = "OMAA-RF",
                              "ENet_ours" = "OG-Clock
(oral)")) +
  labs(x = NULL, y = NULL,
       title = "Oral clock: external validation (Huang 2020)") +
  theme_cell()

## ============================================================
## Panel c — gut transfer predicted-vs-observed scatter (1x3)
## Source: per-sample predictions saved by module4_gut_validation.R
## ============================================================
cat("========== Panel c ==========
")

studies <- c("LifeLinesDeep_2016", "AsnicarF_2021", "KarlssonFH_2013")
pred_fps <- file.path(data_dir, paste0("ext_pred_", studies, ".rds"))
if (!all(file.exists(pred_fps))) {
  stop("Missing per-sample prediction files:
",
       paste0("  ", pred_fps[!file.exists(pred_fps)], collapse = "
"),
       "
These files are written by module4_gut_validation.R inside the cohort ",
       "loop; make sure that script ran to completion.")
}

sc <- do.call(rbind, lapply(studies, function(st) {
  d <- readRDS(file.path(data_dir, paste0("ext_pred_", st, ".rds")))
  data.frame(cohort = st, sample_id = d$sample,
             age = as.numeric(d$obs), pred = as.numeric(d$pred))
}))
cat("per-sample predictions:", nrow(sc), "rows
")
sc$cohort <- factor(sc$cohort, levels = studies)   # fix panel order (matches panel a)

# Strip labels carry each cohort's MAE / rho (from the same summary table as panel a)
lab_df <- out_gut %>%
  filter(cohort != "POOLED") %>%
  mutate(lab = sprintf("%s
n = %d, MAE = %.2f, rho = %.3f",
                       cohort, n, MAE, Spearman))
names(lab_df$lab) <- lab_df$cohort

p_c <- ggplot(sc, aes(age, pred)) +
  # Grey band = training age range (10-80): shows the extrapolation region
  geom_ribbon(data = data.frame(x = c(10, 80)),
              aes(x = x, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "grey92", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey40") +
  geom_point(size = 0.6, alpha = 0.35, color = "#1F77B4") +
  facet_wrap(~cohort, labeller = labeller(cohort = lab_df$lab)) +
  # Fixed scales + equal aspect (coord_equal cannot pair with free scales
  # inside facets)
  coord_equal(xlim = c(15, 85), ylim = c(15, 85)) +
  labs(x = "Chronological age (years)",
       y = "Predicted age (years)",
       title = "Gut clock transfer: predicted vs chronological age") +
  theme_cell()

ggsave(file.path(res_dir, "Figure5c.pdf"), p_c,
       width = 9, height = 3.4, units = "in", device = cairo_pdf)

## ============================================================
## Panel d — transfer gradient overview (11 experiments, 3 groups)
## Numbers from: M1 CV (Table S2), HMP2-native, external validation
## summary, and shotgun-to-16S transfer
## ============================================================
cat("========== Panel d ==========
")

grad <- data.frame(
  group = rep(c("In-distribution CV",
                "Same platform, cross-cohort",
                "Cross-platform (shotgun to 16S)"),
              times = c(3, 5, 3)),
  label = c("Oral clock, CV (Brito+HMP)",
            "Gut clock, CV (Brito+HMP)",
            "Gut clock, HMP2-native",
            "Oral clock, to Huang_2020",
            "Gut clock, to LifeLinesDeep_2016",
            "Gut clock, to AsnicarF_2021",
            "Gut clock, to KarlssonFH_2013",
            "cMD3 to HMP2 (shotgun)",
            "16S, raw mapping",
            "16S, z-score scaling",
            "16S, ComBat"),
  rho   = c(0.416, 0.496, 0.759,
            0.407, 0.167, 0.166, -0.074, 0.210,
            0.118, 0.058, 0.112),
  stringsAsFactors = FALSE
)

# y-axis order: groups top to bottom (factor reversal puts the first group on top)
grad$label <- factor(grad$label, levels = rev(grad$label))
grad$group <- factor(grad$group, levels = c("In-distribution CV",
                                            "Same platform, cross-cohort",
                                            "Cross-platform (shotgun to 16S)"))

p_d <- ggplot(grad, aes(rho, label, color = group)) +
  geom_vline(xintercept = 0,   color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0.3, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_point(size = 2.2, shape = 16) +
  geom_text(aes(label = sprintf("%.3f", rho)),
            hjust = -0.35, size = 2.1, color = "black") +
  scale_color_manual(values = c("In-distribution CV" = "#1F77B4",
                                "Same platform, cross-cohort" = "#FF7F0E",
                                "Cross-platform (shotgun to 16S)" = "#2CA02C"),
                     name = NULL) +
  scale_x_continuous(limits = c(-0.15, 0.85),
                     breaks = seq(-0.2, 0.8, 0.2)) +
  labs(x = expression("Spearman " * rho * " (predicted vs chronological age)"),
       y = NULL,
       title = "Transfer gradient: in-distribution vs cross-cohort vs cross-platform") +
  theme_cell() +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 6.2))

ggsave(file.path(res_dir, "Figure5d.pdf"), p_d,
       width = 6.5, height = 3.6, units = "in", device = cairo_pdf)

## ============================================================
## Assemble — top row a/b, middle row c, bottom row d
## (panel order follows the Figure 5 legend: a = gut bars,
##  b = gut transfer scatter, c = oral bars, d = gradient)
## ============================================================
cat("========== Assembling ==========
")

fig5_v2 <- (p_a | p_c) / p_b / p_d +
  plot_layout(heights = c(1, 1.15, 1.35)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 10, face = "bold"))
  )

ggsave(file.path(res_dir, "Figure5_v2.pdf"), fig5_v2,
       width = 9, height = 10.5, units = "in", device = cairo_pdf)

## ============================================================
## Summary
## ============================================================
cat("
========== module4_figure5 DONE ==========
")
cat("Outputs:
")
cat("  ", file.path(res_dir, "Figure5c.pdf"), "
")
cat("  ", file.path(res_dir, "Figure5d.pdf"), "
")
cat("  ", file.path(res_dir, "Figure5_v2.pdf"), "
")
