# ============================================================================
# module2_figure3.R — Figure 3 (NHANES clinical validation)
#
# Input:
#   <data_dir>/nhanes_clocks.rds
#   <data_dir>/nhanes_cox.rds
#   <data_dir>/nhanes_transfer.rds
#
# Output:
#   <res_dir>/Figure3a.pdf   clock performance scatter (ENet vs OMAA,
#                            discovery/validation)
#   <res_dir>/Figure3b.pdf   survey-weighted KM curves (acceleration quartiles)
#   <res_dir>/Figure3c.pdf   Cox HR forest plot
#   <res_dir>/Figure3d.pdf   cross-platform transfer performance
#   <res_dir>/Figure3.pdf    A+B+C+D composite
#   <res_dir>/Figure3_source_data.xlsx, Figure3_legend.txt
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(survival)
  library(survey)
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

# Colorblind-friendly palettes
pal_method  <- c("Elastic net (ours)" = "#1F77B4",
                 "OMAA-RF (replicated)" = "#FF7F0E")
pal_cycle   <- c("Discovery (2009-2010)" = "#4C72B0",
                 "Validation (2011-2012)" = "#DD8452")
pal_km      <- c("Q1" = "#4C72B0", "Q2" = "#55A868",
                 "Q3" = "#C44E52", "Q4" = "#8172B3")
pal_panel   <- c("All-cause: ENet clock" = "#1F77B4",
                 "All-cause: OMAA-RF" = "#FF7F0E",
                 "CVD mortality" = "#C44E52")

## ============================================================
## 1. Load data
## ============================================================
ck <- readRDS(file.path(data_dir, "nhanes_clocks.rds"))
cx <- readRDS(file.path(data_dir, "nhanes_cox.rds"))
tr <- readRDS(file.path(data_dir, "nhanes_transfer.rds"))

## ============================================================
## 2. Figure 3a — clock performance scatter
## ============================================================
cat("\n========== Figure 3a ==========\n")

plot_3a <- bind_rows(
  data.frame(
    chrono = ck$preds$y_d09,
    pred   = ck$preds$enet_d09_cv,
    method = "Elastic net (ours)",
    cycle  = "Discovery (2009-2010)"
  ),
  data.frame(
    chrono = ck$preds$y_d11,
    pred   = ck$preds$enet_d11,
    method = "Elastic net (ours)",
    cycle  = "Validation (2011-2012)"
  ),
  data.frame(
    chrono = ck$preds$y_d09,
    pred   = ck$preds$omaa_d09_cv,
    method = "OMAA-RF (replicated)",
    cycle  = "Discovery (2009-2010)"
  ),
  data.frame(
    chrono = ck$preds$y_d11,
    pred   = ck$preds$omaa_d11,
    method = "OMAA-RF (replicated)",
    cycle  = "Validation (2011-2012)"
  )
)

plot_3a$method <- factor(plot_3a$method,
                         levels = c("OMAA-RF (replicated)",
                                    "Elastic net (ours)"))
plot_3a$cycle  <- factor(plot_3a$cycle,
                         levels = c("Discovery (2009-2010)",
                                    "Validation (2011-2012)"))

# Per-panel MAE / R2 for annotation
lab_3a <- plot_3a %>%
  group_by(method, cycle) %>%
  summarise(
    MAE = mean(abs(pred - chrono), na.rm = TRUE),
    R2  = 1 - sum((chrono - pred)^2, na.rm = TRUE) /
      sum((chrono - mean(chrono, na.rm = TRUE))^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("MAE = %.1f y\nR\u00B2 = %.2f", MAE, R2))

lim_3a <- range(c(plot_3a$chrono, plot_3a$pred), na.rm = TRUE) + c(-3, 3)

p3a <- ggplot(plot_3a, aes(chrono, pred, color = cycle)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(alpha = 0.35, size = 0.7, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.5) +
  facet_grid(method ~ cycle) +
  scale_color_manual(values = pal_cycle, guide = "none") +
  scale_x_continuous(breaks = seq(20, 80, 20)) +
  scale_y_continuous(breaks = seq(20, 80, 20)) +
  coord_fixed(xlim = lim_3a, ylim = lim_3a) +
  labs(
    x = "Chronological age (years)",
    y = "Predicted age (years)",
    title = "Oral microbiome clocks: discovery vs validation"
  ) +
  theme_cell() +
  geom_text(
    data = lab_3a,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.1,
    size = 2.0, color = "grey20", lineheight = 1.0
  )

ggsave(file.path(res_dir, "Figure3a.pdf"), p3a,
       width = 5.0, height = 4.2, units = "in", device = cairo_pdf)

## ============================================================
## 3. Figure 3b — survey-weighted KM (quartiles)
## ============================================================
cat("\n========== Figure 3b ==========\n")

cc <- cx$cc

cc$q <- factor(ntile(cc$acc_enet_sd, 4),
               levels = 1:4, labels = paste0("Q", 1:4))

des_q <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
  data = cc, nest = TRUE
)

km <- tryCatch(
  svykm(Surv(time, event) ~ q, design = des_q, se = FALSE),
  error = function(e) {
    warning("svykm failed: ", conditionMessage(e)); NULL
  }
)

km_df <- bind_rows(lapply(names(km), function(nm) {
  # nm is either "q=Q1" or "q1"; handle both
  q_lab <- sub(".*=", "", nm)
  if (!q_lab %in% levels(cc$q)) {
    q_lab <- sub(".*q", "Q", nm)
  }
  data.frame(
    time = km[[nm]]$time,
    surv = km[[nm]]$surv,
    q    = q_lab,
    stringsAsFactors = FALSE
  )
}))
km_df$q <- factor(km_df$q, levels = paste0("Q", 1:4))

# Quartile labels
q_lab <- cc %>%
  group_by(q) %>%
  summarise(
    med = median(acc_enet, na.rm = TRUE),
    n   = n(),
    .groups = "drop"
  ) %>%
  mutate(q_lab = sprintf("%s (median acc %+.1f y, n=%d)", q, med, n))

km_df <- km_df %>%
  left_join(q_lab[, c("q", "q_lab")], by = "q")

p3b <- ggplot(km_df, aes(time, surv, color = q_lab)) +
  geom_step(linewidth = 0.7) +
  scale_color_manual(values = setNames(pal_km,
                                       q_lab$q_lab),
                     name = "Clock acceleration quartile") +
  labs(
    x = "Follow-up time (months)",
    y = "Survival probability (survey-weighted)",
    title = "Oral clock acceleration and mortality"
  ) +
  theme_cell() +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white",
                                         color = "grey80",
                                         linewidth = 0.2),
        legend.text = element_text(size = 5.5))

ggsave(file.path(res_dir, "Figure3b.pdf"), p3b,
       width = 4.5, height = 3.8, units = "in", device = cairo_pdf)

## ============================================================
## 4. Figure 3c — Cox HR forest plot
## ============================================================
cat("\n========== Figure 3c ==========\n")

forest_df <- bind_rows(
  cx$cox_tab %>% mutate(panel = case_when(
    model == "ENet-clock"   ~ "All-cause: ENet clock",
    model == "OMAA-RF"      ~ "All-cause: OMAA-RF",
    TRUE ~ model
  )),
  cx$cvd_row %>% mutate(panel = "CVD mortality")
)

forest_df$panel <- factor(forest_df$panel,
                          levels = c("All-cause: ENet clock",
                                     "All-cause: OMAA-RF",
                                     "CVD mortality"))
forest_df$adj_lab <- ifelse(forest_df$adj == "basic",
                            "Basic adj.", "Full adj.")

forest_df <- forest_df %>%
  filter(!is.na(HR), !is.na(HR_low), !is.na(HR_high))

p3c <- ggplot(forest_df, aes(HR, adj_lab, color = panel)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_pointrange(aes(xmin = HR_low, xmax = HR_high),
                  size = 0.4, linewidth = 0.4) +
  facet_wrap(~panel, ncol = 3) +
  scale_color_manual(values = pal_panel, guide = "none") +
  scale_x_log10() +
  labs(
    x = "Hazard ratio per SD of clock acceleration (95% CI)",
    y = NULL,
    title = "Mortality hazard by oral clock acceleration"
  ) +
  theme_cell()

ggsave(file.path(res_dir, "Figure3c.pdf"), p3c,
       width = 6.85, height = 2.6, units = "in", device = cairo_pdf)

## ============================================================
## 5. Figure 3d — cross-platform transfer performance
## ============================================================
cat("\n========== Figure 3d ==========\n")

# Native reference from nhanes_clocks.rds (pooled CV; same evaluation
# set as the transfer experiments)
ck_perf_ref <- ck$perf[ck$perf$cohort == "pooled_CV", ]
ref_rho <- ck_perf_ref$Spearman
tp <- tr$perf %>%
  select(method, MAE, Spearman) %>%
  bind_rows(data.frame(
    method = "NHANES-native (reference)",
    MAE = ck_perf_ref$MAE,
    Spearman = ref_rho
  ))

# Sort by Spearman ascending, then set factor order
tp <- tp %>%
  arrange(Spearman) %>%
  mutate(method = factor(method, levels = method))

# Reference flag
tp$is_ref <- grepl("reference", tp$method)

p3d <- ggplot(tp, aes(method, Spearman, fill = is_ref)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.2) +
  geom_hline(yintercept = ref_rho, linetype = "dashed",
             color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", Spearman)),
            vjust = -0.4, size = 2.2) +
  scale_fill_manual(values = c("FALSE" = "#1F77B4",
                               "TRUE" = "#FF7F0E"),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Spearman rho (predicted vs chronological age)",
    title = "cMD3 shotgun clock transferred to NHANES 16S"
  ) +
  theme_cell()

ggsave(file.path(res_dir, "Figure3d.pdf"), p3d,
       width = 5.5, height = 3.0, units = "in", device = cairo_pdf)

## ============================================================
## 6. Composite Figure 3
## ============================================================
cat("\n========== Composite Figure 3 ==========\n")

# Layout: A top-left, B top-right, C full-width, D bottom
fig3 <- (p3a | p3b) / p3c / p3d +
  plot_layout(heights = c(1.0, 0.7, 0.7)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 10, face = "bold"))
  )

ggsave(file.path(res_dir, "Figure3.pdf"), fig3,
       width = 6.85, height = 10.5, units = "in", device = cairo_pdf)

## ============================================================
## 7. Source Data (Excel)
## ============================================================
cat("\n========== Source Data ==========\n")

wb <- createWorkbook()

# Sheet 1: clock performance
addWorksheet(wb, "3a_clock_perf")
sd_3a <- plot_3a %>%
  mutate(across(where(is.numeric), ~round(., 2)))
writeData(wb, "3a_clock_perf", sd_3a)

# Sheet 2: KM data
addWorksheet(wb, "3b_km")
writeData(wb, "3b_km", km_df)

# Sheet 3: Cox
addWorksheet(wb, "3c_cox")
writeData(wb, "3c_cox", cx$cox_tab)
writeData(wb, "3c_cox", cx$cvd_row, startRow = nrow(cx$cox_tab) + 3)

# Sheet 4: Transfer
addWorksheet(wb, "3d_transfer")
writeData(wb, "3d_transfer", tr$perf)

# Sheet 5: C-index
addWorksheet(wb, "Cindex_PH")
writeData(wb, "Cindex_PH", cx$cindex)

saveWorkbook(wb, file.path(res_dir, "Figure3_source_data.xlsx"),
             overwrite = TRUE)

## ============================================================
## 8. Figure legend
## ============================================================
legend <- "
Figure 3. NHANES oral microbiome clock: performance, mortality association, and cross-platform transfer.

(A) Out-of-fold predicted vs chronological age for two methods (our
    elastic net and the replicated OMAA random forest), across discovery
    (NHANES 2009-2010) and validation (NHANES 2011-2012) cycles. Dashed
    line: y = x. Solid line: linear fit. Colors denote cycle.
(B) Survey-weighted Kaplan-Meier survival curves by quartile of ENet
    clock acceleration. Quartiles defined by ntile() on standardized
    acceleration.
(C) Hazard ratios per SD of clock acceleration for all-cause mortality
    (ENet and OMAA clocks) and CVD mortality (ENet clock), from
    survey-weighted Cox models. Basic adjustment: age, sex, race. Full
    adjustment: + education, income-to-poverty ratio, BMI, smoking.
(D) Cross-platform transfer performance (cMD3 shotgun -> NHANES 16S),
    measured by Spearman rho between predicted and chronological age.
    Reference: NHANES-native clock. Dashed line: native performance.

All P values are two-sided. Source Data: Figure3_source_data.xlsx.
"
writeLines(legend, file.path(res_dir, "Figure3_legend.txt"))

## ============================================================
## 9. Summary output
## ============================================================
cat("\n========== module2_figure3 DONE ==========\n")
cat("Output directory:", res_dir, "\n\n")
for (f in c("Figure3a.pdf", "Figure3b.pdf", "Figure3c.pdf",
            "Figure3d.pdf", "Figure3.pdf",
            "Figure3_source_data.xlsx", "Figure3_legend.txt")) {
  fp <- file.path(res_dir, f)
  if (file.exists(fp)) {
    cat(sprintf("  %-30s %6.1f KB\n", f, file.info(fp)$size / 1024))
  } else {
    cat(sprintf("  %-30s **MISSING**\n", f))
  }
}

cat("\n--- Key numbers ---\n")
cat("Cox N:", cx$n_cc, " | deaths:", cx$n_deaths, "\n")
cat("C-index base:", round(cx$cindex$base, 4),
    " | +ENet:", round(cx$cindex$enet, 4),
    " | +OMAA:", round(cx$cindex$omaa, 4), "\n")
cat("Shared genera:", tr$n_shared, "\n")
