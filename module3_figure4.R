# ============================================================================
# module3_figure4.R — Figure 4 (HMP2 longitudinal stability)
#
# Input:  <data_dir>/hmp2_icc_results.rds
#
# Output:
#   <res_dir>/Figure4a.pdf   within-person residual trajectories (spaghetti)
#   <res_dir>/Figure4b.pdf   ICC forest plot
#   <res_dir>/Figure4.pdf    A+B composite
#   <res_dir>/Figure4_source_data.xlsx, Figure4_legend.txt
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(openxlsx)
})

set.seed(2026)

## ============================================================
## 0. Cell-style theme
## ============================================================
theme_cell <- function(base_size = 7, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      text          = element_text(color = "black"),
      plot.title    = element_text(size = base_size + 1, face = "bold",
                                   hjust = 0, margin = ggplot2::margin(b = 4)),
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

pal_subject <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728",
                 "#9467BD", "#8C564B")

## ============================================================
## 1. Load
## ============================================================
r <- readRDS(file.path(data_dir, "hmp2_icc_results.rds"))
md <- r$md

cat("Samples:", nrow(md), " | subjects:", length(unique(md$subject_id)), "
")

## ============================================================
## 2. Figure 4a — spaghetti plot
## ============================================================
cat("
========== Figure 4a ==========
")

set.seed(2026)
subs <- unique(md$subject_id)
samp <- sample(subs, min(30, length(subs)))
highlight_subs <- sample(samp, min(6, length(samp)))

sp_df <- md %>%
  mutate(subject_id = as.character(subject_id)) %>%
  filter(subject_id %in% samp) %>%
  mutate(
    highlight = subject_id %in% highlight_subs,
    highlight_lab = ifelse(highlight, subject_id, NA)
  )

# Prefer interval_days; fall back to week_num * 7 when missing
if (!"interval_days" %in% colnames(sp_df) || all(is.na(sp_df$interval_days))) {
  sp_df$x_days <- sp_df$week_num * 7
  cat("  interval_days missing; using week_num * 7
")
} else {
  sp_df$x_days <- sp_df$interval_days
}

p4a <- ggplot(sp_df, aes(x_days, oof_resid, group = subject_id)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_line(data = filter(sp_df, !highlight),
            alpha = 0.25, color = "grey35", linewidth = 0.3) +
  geom_line(data = filter(sp_df, highlight),
            aes(color = highlight_lab), linewidth = 0.7) +
  geom_point(data = filter(sp_df, highlight),
             aes(color = highlight_lab), size = 1.0, stroke = 0) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE,
              color = "#8172B3", linewidth = 0.7, linetype = "dashed") +
  scale_color_manual(values = pal_subject, guide = "none") +
  labs(
    x = "Days since baseline sample",
    y = "Clock acceleration (years)",
    title = "Within-person trajectory of gut clock acceleration (HMP2)"
  ) +
  theme_cell()

ggsave(file.path(res_dir, "Figure4a.pdf"), p4a,
       width = 5.5, height = 3.5, units = "in", device = cairo_pdf)

## ============================================================
## 3. Figure 4b — ICC forest plot
## ============================================================
cat("
========== Figure 4b ==========
")

# Sort by ICC before setting the factor order
ic <- r$icc_tab %>%
  filter(!is.na(ICC)) %>%
  arrange(ICC) %>%
  mutate(comparison = factor(comparison, levels = comparison))

# Categories
ic$category <- case_when(
  grepl("Clock residual.*native", ic$comparison) ~ "Primary",
  grepl("Clock residual.*transferred", ic$comparison) ~ "Sensitivity",
  grepl("Predicted age", ic$comparison) ~ "Reference",
  grepl("Shannon", ic$comparison) ~ "Benchmark",
  grepl("Permuted", ic$comparison) ~ "Control",
  grepl("CD|UC|nonIBD", ic$comparison) ~ "Stratified",
  TRUE ~ "Other"
)

p4b <- ggplot(ic, aes(ICC, comparison, color = category)) +
  geom_vline(xintercept = 0, linetype = "solid",
             color = "grey55", linewidth = 0.3) +
  geom_pointrange(aes(xmin = low, xmax = high),
                  size = 0.4, linewidth = 0.4) +
  geom_text(aes(x = high, label = sprintf("%.2f [%.2f-%.2f]",
                                          ICC, low, high)),
            hjust = -0.05, size = 2.0, color = "grey25") +
  scale_color_manual(values = c(
    "Primary"    = "#1F77B4",
    "Sensitivity"= "#C44E52",
    "Reference"  = "#2CA02C",
    "Benchmark"  = "#8172B3",
    "Control"    = "#7F7F7F",
    "Stratified" = "#FF7F0E"
  ), name = NULL) +
  xlim(-0.1, 1.35) +
  labs(
    x = "Intraclass correlation (95% CI, subject-level bootstrap)",
    y = NULL,
    title = "Longitudinal stability of clock residual"
  ) +
  theme_cell() +
  theme(legend.position = "right",
        legend.text = element_text(size = 6))

ggsave(file.path(res_dir, "Figure4b.pdf"), p4b,
       width = 6.5, height = 3.8, units = "in", device = cairo_pdf)

## ============================================================
## 4. Composite
## ============================================================
fig4 <- p4a / p4b +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 10, face = "bold"))
  )

ggsave(file.path(res_dir, "Figure4.pdf"), fig4,
       width = 6.85, height = 7.5, units = "in", device = cairo_pdf)

## ============================================================
## 5. Source Data
## ============================================================
wb <- createWorkbook()

addWorksheet(wb, "4a_spaghetti")
writeData(wb, "4a_spaghetti",
          sp_df %>% select(subject_id, x_days, oof_resid,
                           age_est, diagnosis))

addWorksheet(wb, "4b_icc")
writeData(wb, "4b_icc", r$icc_tab)

addWorksheet(wb, "4_perf")
writeData(wb, "4_perf", data.frame(
  metric = c("HMP2-native 5x5 OOF MAE",
             "HMP2-native 5x5 OOF rho",
             "Transferred pred-age rho",
             "N samples", "N subjects",
             "Median per-subject resid SD"),
  value = c(r$native_perf["MAE"],
            r$native_perf["Spearman"],
            r$transfer_rho,
            r$n_samples,
            r$n_subjects,
            median(r$aux$sd_resid, na.rm = TRUE))
))

saveWorkbook(wb, file.path(res_dir, "Figure4_source_data.xlsx"),
             overwrite = TRUE)

## ============================================================
## 6. Figure legend
## ============================================================
legend <- sprintf(
  "Figure 4. Longitudinal stability of gut microbiome clock acceleration in HMP2.

(A) Within-person trajectory of clock acceleration (years) vs days since
    baseline for 30 randomly sampled subjects; 6 subjects highlighted.
    Dashed line: LOESS smoothing across all subjects.
(B) Intraclass correlation (ICC) of clock residuals and reference
    metrics. ICC computed via one-way random-effects ANOVA (method of
    moments, unbalanced-aware). 95%% CI from subject-level bootstrap
    (%d resamples). Primary: HMP2-native clock (out-of-fold residuals,
    5x5 CV, MAE = %.2f years, rho = %.3f). Sensitivity: cMD3-transferred
    clock (pred-age rho = %.2f; degenerate ICC). Permuted residual as
    negative control.

Data: %d samples from %d subjects.
",
  500,
  r$native_perf["MAE"], r$native_perf["Spearman"],
  r$transfer_rho,
  r$n_samples, r$n_subjects
)
writeLines(legend, file.path(res_dir, "Figure4_legend.txt"))

## ============================================================
## 7. Summary output
## ============================================================
cat("
========== module3_figure4 DONE ==========
")
for (f in c("Figure4a.pdf", "Figure4b.pdf", "Figure4.pdf",
            "Figure4_source_data.xlsx", "Figure4_legend.txt")) {
  fp <- file.path(res_dir, f)
  if (file.exists(fp)) {
    cat(sprintf("  %-30s %6.1f KB
", f, file.info(fp)$size / 1024))
  } else {
    cat(sprintf("  %-30s **MISSING**
", f))
  }
}

cat("
Key numbers:
")
cat("  HMP2-native 5x5 OOF MAE:", round(r$native_perf["MAE"], 2), "
")
cat("  HMP2-native 5x5 OOF rho:", round(r$native_perf["Spearman"], 3), "
")
cat("  Transferred pred-age rho:", round(r$transfer_rho, 3), "
")
cat("  ICC (primary):",
    round(r$icc_tab$ICC[r$icc_tab$comparison == "Clock residual (HMP2-native, OOF)"], 3),
    "
")
