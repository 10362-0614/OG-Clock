# ============================================================================
# module2_summary_table.R — NHANES results summary table
#
# Input:
#   <data_dir>/nhanes_clocks.rds
#   <data_dir>/nhanes_transfer.rds
#   <data_dir>/nhanes_cox.rds
#
# Output:
#   <res_dir>/nhanes_summary.csv
#   <res_dir>/nhanes_summary.xlsx (per-block sheets)
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

## ============================================================
## 1. Check input files
## ============================================================
cat("\n========== 1. Check inputs ==========\n")

need <- c(
  clock    = file.path(data_dir, "nhanes_clocks.rds"),
  transfer = file.path(data_dir, "nhanes_transfer.rds"),
  cox      = file.path(data_dir, "nhanes_cox.rds")
)

missing <- need[!file.exists(need)]
if (length(missing) > 0) {
  stop("Missing input files:\n",
       paste(sprintf("  - %s", missing), collapse = "\n"),
       "\nRun module2_clock_benchmark.R / module2_cross_platform_transfer.R / ",
       "module2_cox_mortality.R first")
}

clock <- readRDS(need["clock"])
trf   <- readRDS(need["transfer"])
cox   <- readRDS(need["cox"])

cat("clock:    ", nrow(clock$perf), "rows\n")
cat("transfer: ", nrow(trf$perf),   "rows\n")
cat("cox:      ", nrow(cox$cox_tab), "rows\n")

## ============================================================
## 2. Normalize column structure per block
## ============================================================
# Unified columns: block, model, adj, exposure, HR, HR_low, HR_high, p,
# MAE, Spearman, R2, Pearson, n, deaths; missing columns filled with NA

make_block <- function(df, block_name, cols_order) {
  for (cc in cols_order) {
    if (!cc %in% colnames(df)) df[[cc]] <- NA
  }
  df <- df[, cols_order, drop = FALSE]
  df$block <- block_name
  df[, c("block", cols_order), drop = FALSE]
}

cols_order <- c("model", "adj", "exposure",
                "HR", "HR_low", "HR_high", "p",
                "MAE", "Spearman", "R2", "Pearson",
                "n", "deaths")

## ============================================================
## 3. Block 1: clock performance (ENet vs OMAA-RF)
## ============================================================
cat("\n========== 3. Clock performance ==========\n")

perf_clock <- clock$perf %>%
  mutate(
    model = method,
    adj   = cohort,
    HR = NA_real_, HR_low = NA_real_, HR_high = NA_real_, p = NA_real_,
    Pearson = as.numeric(Pearson),
    n = NA_integer_, deaths = NA_integer_
  )

block1 <- make_block(perf_clock, "1_clock_performance", cols_order)

print(block1, row.names = FALSE, digits = 3)

## ============================================================
## 4. Block 2: cross-platform transfer
## ============================================================
cat("\n========== 4. Cross-platform transfer ==========\n")

perf_trf <- trf$perf %>%
  mutate(
    model = method,
    adj   = paste0("shared_genera=", n_shared),
    exposure = "cMD3 -> NHANES",
    HR = NA_real_, HR_low = NA_real_, HR_high = NA_real_, p = NA_real_,
    R2 = NA_real_,
    n = NA_integer_, deaths = NA_integer_
  )

block2 <- make_block(perf_trf, "2_transfer", cols_order)

print(block2, row.names = FALSE, digits = 3)

## ============================================================
## 5. Block 3: Cox (all-cause / CVD mortality)
## ============================================================
cat("\n========== 5. Cox results ==========\n")

cx <- bind_rows(cox$cox_tab, cox$cvd_row) %>%
  mutate(
    MAE = NA_real_, Spearman = NA_real_, R2 = NA_real_, Pearson = NA_real_
  )

block3 <- make_block(cx, "3_cox", cols_order)

print(block3, row.names = FALSE, digits = 4)

## ============================================================
## 6. Block 3 addendum: C-index increments + PH tests
## ============================================================
cat("\n========== 6. C-index + PH ==========\n")

cidx <- data.frame(
  model = c("C-index: base",
            "C-index: + ENet",
            "C-index: + OMAA",
            "dC-index: + ENet",
            "dC-index: + OMAA",
            "PH p: ENet accel",
            "PH p: OMAA accel"),
  adj   = c("age + covariates", "age + covariates", "age + covariates",
            "vs base", "vs base",
            "full model", "full model"),
  value = c(cox$cindex$base,
            cox$cindex$enet,
            cox$cindex$omaa,
            cox$cindex$delta_enet,
            cox$cindex$delta_omaa,
            unname(cox$ph_p["ENet"]),
            unname(cox$ph_p["OMAA"])),
  stringsAsFactors = FALSE
)

print(cidx, row.names = FALSE, digits = 4)

# Append C-index rows to block3
block3_extra <- data.frame(
  model = paste0(cidx$model, " (", cidx$adj, ")"),
  adj   = cidx$adj,
  exposure = "c-index / PH",
  HR = NA_real_, HR_low = NA_real_, HR_high = NA_real_,
  p = NA_real_,
  MAE = NA_real_, Spearman = NA_real_, R2 = NA_real_, Pearson = NA_real_,
  n = NA_integer_, deaths = NA_integer_
)

block3 <- bind_rows(block3, make_block(block3_extra, "3_cox", cols_order))

## ============================================================
## 7. Block 4: descriptives (population flow + acceleration spread)
## ============================================================
cat("\n========== 7. Descriptives ==========\n")

cc <- cox$cc
n_cc <- nrow(cc)
n_deaths <- sum(cc$event, na.rm = TRUE)
n_cvd <- sum(cc$cvd, na.rm = TRUE)
median_fu <- median(cc$time[cc$event == 0], na.rm = TRUE)

sd_enet <- sd(cc$acc_enet, na.rm = TRUE)
sd_omaa <- sd(cc$acc_omaa, na.rm = TRUE)

ref_native <- trf$ref_native

block4 <- data.frame(
  block   = "4_descriptives",
  model   = c("N (Cox population)",
              "Deaths (all-cause)",
              "Deaths (CVD)",
              "Median follow-up (months, alive)",
              "ENet accel SD (years)",
              "OMAA accel SD (years)",
              "Shared genera (cMD3 vs NHANES)",
              "NHANES-native clock (reference)"),
  adj     = NA_character_,
  exposure = NA_character_,
  HR = NA_real_, HR_low = NA_real_, HR_high = NA_real_,
  p = NA_real_,
  MAE = NA_real_, Spearman = NA_real_, R2 = NA_real_, Pearson = NA_real_,
  n = NA_integer_, deaths = NA_integer_,
  value = c(n_cc, n_deaths, n_cvd, median_fu,
            sd_enet, sd_omaa,
            trf$n_shared,
            NA),
  label = c(as.character(n_cc), as.character(n_deaths),
            as.character(n_cvd), sprintf("%.1f", median_fu),
            sprintf("%.2f", sd_enet), sprintf("%.2f", sd_omaa),
            as.character(trf$n_shared),
            ref_native),
  stringsAsFactors = FALSE
)

print(block4[, c("model", "label")], row.names = FALSE)

## ============================================================
## 8. Combine blocks
## ============================================================
summary_tab <- bind_rows(
  block1 %>% mutate(value = NA_real_, label = NA_character_),
  block2 %>% mutate(value = NA_real_, label = NA_character_),
  block3 %>% mutate(value = NA_real_, label = NA_character_),
  block4
)

final_cols <- c("block", "model", "adj", "exposure",
                "HR", "HR_low", "HR_high", "p",
                "MAE", "Spearman", "R2", "Pearson",
                "n", "deaths",
                "value", "label")
summary_tab <- summary_tab[, final_cols, drop = FALSE]

## ============================================================
## 9. Format values for display
## ============================================================
fmt_tab <- summary_tab %>%
  mutate(
    HR        = ifelse(is.na(HR), NA, sprintf("%.3f", HR)),
    HR_low    = ifelse(is.na(HR_low), NA, sprintf("%.3f", HR_low)),
    HR_high   = ifelse(is.na(HR_high), NA, sprintf("%.3f", HR_high)),
    p         = ifelse(is.na(p), NA, formatC(p, format = "g", digits = 3)),
    MAE       = ifelse(is.na(MAE), NA, sprintf("%.2f", MAE)),
    Spearman  = ifelse(is.na(Spearman), NA, sprintf("%.3f", Spearman)),
    R2        = ifelse(is.na(R2), NA, sprintf("%.3f", R2)),
    Pearson   = ifelse(is.na(Pearson), NA, sprintf("%.3f", Pearson))
  )

## ============================================================
## 10. Output
## ============================================================
csv_path <- file.path(res_dir, "nhanes_summary.csv")
write.csv(summary_tab, csv_path, row.names = FALSE, na = "")
cat("\nCSV saved:", csv_path, "\n")

xlsx_path <- file.path(res_dir, "nhanes_summary.xlsx")
wb <- createWorkbook()

addWorksheet(wb, "All")
writeData(wb, "All", fmt_tab)
freezePane(wb, "All", firstRow = TRUE)

addWorksheet(wb, "Clock_performance")
writeData(wb, "Clock_performance", block1)

addWorksheet(wb, "Transfer")
writeData(wb, "Transfer", block2)

addWorksheet(wb, "Cox")
writeData(wb, "Cox", block3)

addWorksheet(wb, "Descriptives")
writeData(wb, "Descriptives", block4)

addWorksheet(wb, "Cindex_PH")
writeData(wb, "Cindex_PH", cidx)

saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat("Excel saved:", xlsx_path, "\n")

## ============================================================
## 11. Print full summary
## ============================================================
cat("\n========== NHANES summary ==========\n")
print(fmt_tab[, c("block", "model", "adj",
                  "HR", "p", "MAE", "Spearman", "R2")],
      row.names = FALSE)

cat("\n========== module2_summary_table DONE ==========\n")
