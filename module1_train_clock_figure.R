# ============================================================================
# module1_train_clock_figure.R — Figure 2 (Module 1 clock performance)
#
# Panels:
#   (A) Out-of-fold predicted vs chronological age, oral/gut/joint clocks,
#       5x5 repeated CV, colored by cohort
#   (B) MAE bar chart, 5x5 CV vs LOSO
#
# Output:
#   <res_dir>/Figure2a.pdf, Figure2b.pdf, Figure2.pdf (A+B composite)
#   <res_dir>/Figure2_source_data.xlsx, Figure2_legend.txt
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(openxlsx)
})

## ============================================================
## 0. Cell-style theme
## ============================================================
theme_cell <- function(base_size = 7, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      text          = element_text(color = "black"),
      plot.title    = element_text(size = base_size + 1, face = "bold",
                                   hjust = 0, margin = margin(b = 4)),
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
      plot.margin   = margin(4, 4, 4, 4)
    )
}

# Colorblind-friendly cohort palette; names must exactly match the
# relabeled factor levels of df$study (with spaces)
cell_pal_study <- c("BritoIL 2016" = "#1F77B4",
                    "HMP 2012"     = "#FF7F0E")

cell_pal_valid <- c("5x5 CV" = "#4C72B0",
                    "LOSO"   = "#DD8452")

## ============================================================
## 1. Load model results
## ============================================================
results <- readRDS(file.path(data_dir, "03_models.rds"))

y      <- results$y
study  <- results$study
n_subj <- length(y)

## ============================================================
## 2. Assemble plotting data frame
## ============================================================
mk_df <- function(res, scheme) {
  data.frame(
    chrono     = y,
    pred       = res$pred,
    study      = study,
    clock      = res$label,
    validation = scheme
  )
}

df_cv <- bind_rows(
  mk_df(results$cv$oral,  "5x5 CV"),
  mk_df(results$cv$gut,   "5x5 CV"),
  mk_df(results$cv$joint, "5x5 CV")
)
df_loso <- bind_rows(
  mk_df(results$loso$oral,  "LOSO"),
  mk_df(results$loso$gut,   "LOSO"),
  mk_df(results$loso$joint, "LOSO")
)
df <- bind_rows(df_cv, df_loso)
df$clock <- gsub("-LOSO$", "", df$clock)
df$clock <- factor(df$clock,
                   levels = c("oral-clock", "gut-clock", "joint-clock"),
                   labels = c("Oral", "Gut", "Joint"))
df$study <- factor(df$study,
                   levels = c("BritoIL_2016", "HMP_2012"),
                   labels = c("BritoIL 2016", "HMP 2012"))

# Defensive check: palette names must match factor levels exactly
stopifnot(setequal(names(cell_pal_study), levels(df$study)))

## ============================================================
## 3. Figure 2a — predicted vs chronological age (three clocks, CV)
## ============================================================
lim <- range(c(df$chrono, df$pred), na.rm = TRUE) + c(-2, 2)

# Per-clock CV MAE / R2 for panel annotation
panel_lab <- df %>%
  filter(validation == "5x5 CV") %>%
  group_by(clock) %>%
  summarise(
    MAE = mean(abs(pred - chrono)),
    R2  = 1 - sum((chrono - pred)^2) / sum((chrono - mean(chrono))^2),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("MAE = %.1f y
R\u00B2 = %.2f", MAE, R2))

df_cv_only <- df %>% filter(validation == "5x5 CV")

p2a <- ggplot(df_cv_only, aes(chrono, pred, color = study)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(alpha = 0.6, size = 1.0, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.5, alpha = 0.9) +
  facet_wrap(~clock, nrow = 1) +
  scale_color_manual(values = cell_pal_study, name = "Cohort") +
  scale_x_continuous(breaks = seq(10, 80, 20)) +
  scale_y_continuous(breaks = seq(10, 80, 20)) +
  coord_fixed(xlim = lim, ylim = lim) +
  labs(
    x = "Chronological age (years)",
    y = "Predicted age (years)",
    title = "Microbiome aging clocks: 5\u00D75 cross-validation"
  ) +
  theme_cell() +
  # Legend on top reads better under a three-panel facet row
  theme(legend.position = "top") +
  geom_text(
    data = panel_lab,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.1,
    size = 2.2, color = "grey20", lineheight = 1.0
  )

## ============================================================
## 4. Figure 2b — MAE bar chart (CV vs LOSO)
## ============================================================
mae_bar <- df %>%
  group_by(clock, validation) %>%
  summarise(MAE = mean(abs(pred - chrono)), .groups = "drop") %>%
  mutate(validation = factor(validation, c("5x5 CV", "LOSO")))

p2b <- ggplot(mae_bar, aes(clock, MAE, fill = validation)) +
  geom_col(position = position_dodge(0.75), width = 0.65,
           color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f", MAE)),
            position = position_dodge(0.75),
            vjust = -0.4, size = 2.2) +
  scale_fill_manual(values = cell_pal_valid, name = "Validation") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(
    x = NULL,
    y = "MAE (years)",
    title = "Clock accuracy: within-distribution vs cross-cohort"
  ) +
  theme_cell() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = NA)
  )

## ============================================================
## 5. Save individual panels (Cell style: 1-col 85 mm = 3.35 in,
##    2-col 174 mm = 6.85 in; extra height reserves room for top legend)
## ============================================================
ggsave(file.path(res_dir, "Figure2a.pdf"), p2a,
       width = 6.85, height = 2.8, units = "in", device = cairo_pdf)

ggsave(file.path(res_dir, "Figure2b.pdf"), p2b,
       width = 3.35, height = 3.2, units = "in", device = cairo_pdf)

## ============================================================
## 6. Composite Figure 2 (panels A/B)
## ============================================================
fig2 <- (p2a / (plot_spacer() + p2b + plot_spacer() +
                  plot_layout(widths = c(1, 2, 1)))) +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 10, face = "bold"))
  )

ggsave(file.path(res_dir, "Figure2.pdf"), fig2,
       width = 6.85, height = 6.7, units = "in", device = cairo_pdf)

## ============================================================
## 7. Source Data (Excel)
## ============================================================
wb <- createWorkbook()

# Sheet 1: CV + LOSO predictions
addWorksheet(wb, "Predictions")
pred_out <- df %>% mutate(across(where(is.numeric), ~round(., 3)))
writeData(wb, "Predictions", pred_out)

# Sheet 2: MAE summary
addWorksheet(wb, "MAE_summary")
writeData(wb, "MAE_summary", mae_bar %>% mutate(MAE = round(MAE, 2)))

# Sheet 3: Complementarity
addWorksheet(wb, "Complementarity")
writeData(wb, "Complementarity",
          results$complementarity$table %>% mutate(value = round(value, 3)))

# Sheet 4: Residual correlations
addWorksheet(wb, "Residual_cor")
res_o <- results$cv$oral$pred - y
res_g <- results$cv$gut$pred  - y
rc <- data.frame(
  metric = c("pooled_spearman", "BritoIL_2016", "HMP_2012"),
  value  = c(
    cor(res_o, res_g, method = "spearman"),
    cor(res_o[study == "BritoIL_2016"],
        res_g[study == "BritoIL_2016"], method = "spearman"),
    cor(res_o[study == "HMP_2012"],
        res_g[study == "HMP_2012"], method = "spearman")
  )
)
writeData(wb, "Residual_cor", rc %>% mutate(value = round(value, 3)))

saveWorkbook(wb, file.path(res_dir, "Figure2_source_data.xlsx"),
             overwrite = TRUE)

## ============================================================
## 8. Figure legend
## ============================================================
legend <- sprintf(
  "Figure 2. Paired oral-gut microbiome aging clocks.
(A) Out-of-fold predicted vs chronological age for oral, gut, and joint
    clocks under 5\u00D75 repeated cross-validation (n = %d subjects, %d
    from BritoIL_2016, %d from HMP_2012). Dashed line: y = x. Solid lines:
    linear fit per cohort. Colors denote cohort.
(B) Mean absolute error (MAE) of each clock under 5\u00D75 CV
    (within-distribution) and leave-one-study-out (LOSO; cross-cohort).

Source Data: Figure2_source_data.xlsx contains all underlying numerical
values, including out-of-fold predictions, MAE, complementarity metrics,
and residual correlations.
",
  n_subj,
  sum(study == "BritoIL_2016"),
  sum(study == "HMP_2012")
)
writeLines(legend, file.path(res_dir, "Figure2_legend.txt"))

## ============================================================
## 9. Summary output
## ============================================================
cat("
========== module1_train_clock_figure DONE ==========
")
cat("Output directory:", res_dir, "

")
cat("Files:
")
for (f in c("Figure2a.pdf", "Figure2b.pdf", "Figure2.pdf",
            "Figure2_source_data.xlsx", "Figure2_legend.txt")) {
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
MAE summary:
")
print(mae_bar, row.names = FALSE, digits = 2)

cat("
Residual correlations:
")
print(rc, row.names = FALSE, digits = 3)
