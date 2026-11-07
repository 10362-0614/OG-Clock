# ============================================================================
# module2_cox_mortality.R — All-cause and CVD mortality Cox models
#
# Input:  <data_dir>/nhanes_clocks.rds (from module2_clock_benchmark.R)
# Output: <data_dir>/nhanes_cox.rds = list(cox_tab, cvd_row, cindex, cc,
#                                          km_tab, km, ph_p, ...)
#         <res_dir>/nhanes_cox_results.csv
#
# Design:
#   - Survey-weighted (SDMVPSU / SDMVSTRA / WTMEC2YR) svycoxph primary
#     analysis
#   - Unweighted coxph (weighted) used for C-index and PH tests, explicitly
#     labeled as unweighted approximations (svycoxph does not support
#     concordance() / cox.zph)
#   - Acceleration variables SD-standardized (HR reported per SD)
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(survey)
  library(survival)
})

## ============================================================
## 1. Load clock results
## ============================================================
cat("\n========== 1. Load data ==========\n")

r <- readRDS(file.path(data_dir, "nhanes_clocks.rds"))

dat <- r$dat
dat$acc_enet <- r$acc_enet
dat$acc_omaa <- r$acc_omaa

cat("Total samples:", nrow(dat), "\n")
cat("ENet acceleration non-NA:", sum(!is.na(dat$acc_enet)), "\n")
cat("OMAA acceleration non-NA:", sum(!is.na(dat$acc_omaa)), "\n")

## ============================================================
## 2. Cox analysis population
## ============================================================
cat("\n========== 2. Cox analysis population ==========\n")

cc <- dat %>%
  filter(
    cc_complete,                       # complete covariates
    !is.na(acc_enet),
    !is.na(acc_omaa),
    !is.na(time_m), time_m > 0,
    !is.na(MORTSTAT),
    ELIGSTAT == 1
  )

cat("Cox analysis population:", nrow(cc), "\n")
cat("Deaths:", sum(cc$MORTSTAT == 1, na.rm = TRUE), "\n")
cat("Median follow-up (alive):",
    median(cc$time_m[cc$MORTSTAT == 0], na.rm = TRUE), "months\n")
cat("Median follow-up (deceased):",
    median(cc$time_m[cc$MORTSTAT == 1], na.rm = TRUE), "months\n")

# SD-standardize acceleration variables
cc <- cc %>%
  mutate(
    acc_enet_sd = as.numeric(scale(acc_enet)),
    acc_omaa_sd = as.numeric(scale(acc_omaa)),
    event       = as.numeric(MORTSTAT == 1),
    time        = time_m,
    sex         = factor(RIAGENDR),
    race        = factor(RIDRETH1),
    educ        = factor(DMDEDUC2),
    smoke       = factor(smoke, c("never", "former", "current"))
  )

cat("\nAcceleration SD (years):\n")
cat("  ENet:", round(sd(cc$acc_enet), 2), "\n")
cat("  OMAA:", round(sd(cc$acc_omaa), 2), "\n")

cat("\nCovariate missingness:\n")
for (v in c("sex", "race", "educ", "smoke", "INDFMPIR", "BMXBMI")) {
  cat(sprintf("  %-10s missing: %d (%.1f%%)\n",
              v, sum(is.na(cc[[v]])), 100 * mean(is.na(cc[[v]]))))
}

## ============================================================
## 3. Unweighted Cox (crude diagnostics)
## ============================================================
cat("\n========== 3. Unweighted Cox (crude diagnostics) ==========\n")

cat("--- ENet acceleration (per SD) ---\n")
fit_crude_e <- coxph(Surv(time, event) ~ acc_enet_sd, data = cc)
print(summary(fit_crude_e)$coefficients[1, c("coef", "exp(coef)", "Pr(>|z|)")])

cat("\n--- OMAA acceleration (per SD) ---\n")
fit_crude_o <- coxph(Surv(time, event) ~ acc_omaa_sd, data = cc)
print(summary(fit_crude_o)$coefficients[1, c("coef", "exp(coef)", "Pr(>|z|)")])

## ============================================================
## 4. Survey-weighted design
## ============================================================
cat("\n========== 4. Survey design ==========\n")

des <- svydesign(
  ids     = ~SDMVPSU,
  strata  = ~SDMVSTRA,
  weights = ~WTMEC2YR,
  data    = cc,
  nest    = TRUE
)

cat("Survey design created\n")

## ============================================================
## 5. svycoxph primary analysis
## ============================================================
cat("\n========== 5. svycoxph primary analysis ==========\n")

fit_cox2 <- function(design, exposure, tag) {
  f_basic <- as.formula(sprintf(
    "Surv(time, event) ~ %s + RIDAGEYR + sex + race", exposure
  ))
  f_full <- as.formula(sprintf(
    "Surv(time, event) ~ %s + RIDAGEYR + sex + race + educ + INDFMPIR + BMXBMI + smoke",
    exposure
  ))
  
  res <- list()
  
  for (nm in c("basic", "full")) {
    f <- if (nm == "basic") f_basic else f_full
    
    fit <- tryCatch(
      svycoxph(f, design = design),
      error = function(e) {
        warning(tag, " ", nm, ": ", conditionMessage(e))
        NULL
      }
    )
    
    if (is.null(fit)) {
      res[[nm]] <- data.frame(
        model = tag, adj = nm, exposure = exposure,
        HR = NA, HR_low = NA, HR_high = NA, p = NA,
        n = NA, deaths = NA
      )
      next
    }
    
    s  <- summary(fit)$coefficients[1, ]
    ci <- tryCatch(confint(fit)[1, ], error = function(e) c(NA, NA))
    pcol <- names(s)[grepl("Pr", names(s))][1]
    
    res[[nm]] <- data.frame(
      model    = tag,
      adj      = nm,
      exposure = exposure,
      HR       = exp(s[1]),
      HR_low   = exp(ci[1]),
      HR_high  = exp(ci[2]),
      p        = s[pcol],
      n        = nrow(design$variables),
      deaths   = sum(design$variables$event, na.rm = TRUE)
    )
  }
  
  # PH test via unweighted coxph (svycoxph does not support cox.zph);
  # explicitly an unweighted approximation
  fit_u <- tryCatch(
    coxph(f_full,
          data = as.data.frame(design$variables),
          weights = WTMEC2YR),
    error = function(e) NULL
  )
  
  if (!is.null(fit_u)) {
    res$ph_p <- tryCatch(
      cox.zph(fit_u)$table[1, "p"],
      error = function(e) NA
    )
  } else {
    res$ph_p <- NA
  }
  
  res
}

res_e <- fit_cox2(des, "acc_enet_sd", "ENet-clock")
res_o <- fit_cox2(des, "acc_omaa_sd", "OMAA-RF")

cox_tab <- bind_rows(res_e$basic, res_e$full, res_o$basic, res_o$full)

cat("\nCox results:\n")
print(cox_tab, row.names = FALSE, digits = 4)

cat(sprintf("\nPH assumption p-values (full model, unweighted approximation):\n"))
cat(sprintf("  ENet: %.4f\n", res_e$ph_p))
cat(sprintf("  OMAA: %.4f\n", res_o$ph_p))

if (!is.na(res_e$ph_p) && res_e$ph_p < 0.05) {
  warning("ENet acceleration violates the PH assumption (p < 0.05); ",
          "consider a time-varying term or stratification")
}
if (!is.na(res_o$ph_p) && res_o$ph_p < 0.05) {
  warning("OMAA acceleration violates the PH assumption (p < 0.05); ",
          "consider a time-varying term or stratification")
}

## ============================================================
## 6. CVD mortality (secondary analysis)
## ============================================================
cat("\n========== 6. CVD mortality ==========\n")

# UCOD codes: 001 = Diseases of heart, 005 = Cerebrovascular diseases
# (113 cause-of-death codes, 2019 public-release data dictionary)
cc$cvd <- as.numeric(cc$MORTSTAT == 1 & cc$UCOD %in% c("001", "005"))

cat("CVD deaths:", sum(cc$cvd, na.rm = TRUE), "\n")
cat("CVD share of deaths:", round(100 * sum(cc$cvd, na.rm = TRUE) /
                                    sum(cc$event, na.rm = TRUE), 1), "%\n")

cvd_row <- data.frame(
  model = "ENet-clock", adj = "full", exposure = "CVD mortality",
  HR = NA, HR_low = NA, HR_high = NA, p = NA, n = NA, deaths = NA
)

if (sum(cc$cvd, na.rm = TRUE) >= 30) {
  des_cvd <- svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
    data = cc, nest = TRUE
  )
  
  fit_cvd <- tryCatch(
    svycoxph(
      Surv(time, cvd) ~ acc_enet_sd + RIDAGEYR + sex + race +
        educ + INDFMPIR + BMXBMI + smoke,
      design = des_cvd
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit_cvd)) {
    s  <- summary(fit_cvd)$coefficients[1, ]
    ci <- tryCatch(confint(fit_cvd)[1, ], error = function(e) c(NA, NA))
    pcol <- names(s)[grepl("Pr", names(s))][1]
    
    cvd_row <- data.frame(
      model = "ENet-clock", adj = "full", exposure = "CVD mortality",
      HR = exp(s[1]), HR_low = exp(ci[1]), HR_high = exp(ci[2]),
      p = s[pcol], n = nrow(cc), deaths = sum(cc$cvd, na.rm = TRUE)
    )
    
    cat("\nCVD Cox result:\n")
    print(cvd_row, row.names = FALSE, digits = 4)
  }
} else {
  cat("CVD deaths < 30; underpowered, Cox skipped\n")
}

## ============================================================
## 7. C-index increment (unweighted approximation)
## ============================================================
cat("\n========== 7. C-index increment ==========\n")
cat("Note: C-index based on unweighted coxph (svycoxph does not support ",
    "concordance); a survey-weighted C-index would require svycontrast ",
    "or bootstrap and is not implemented here.\n\n")

# Unweighted coxph with x=TRUE for concordance
fit1 <- coxph(
  Surv(time, event) ~ RIDAGEYR + sex + race + educ + INDFMPIR + BMXBMI + smoke,
  data = cc, weights = WTMEC2YR, x = TRUE
)
fit2 <- update(fit1, . ~ . + acc_enet_sd)
fit3 <- update(fit1, . ~ . + acc_omaa_sd)

c1 <- concordance(fit1)$concordance
c2 <- concordance(fit2)$concordance
c3 <- concordance(fit3)$concordance

cat(sprintf("C-index:\n"))
cat(sprintf("  base (age + covariates):  %.4f\n", c1))
cat(sprintf("  base + ENet acceleration: %.4f (delta = %+.4f)\n", c2, c2 - c1))
cat(sprintf("  base + OMAA acceleration: %.4f (delta = %+.4f)\n", c3, c3 - c1))

cindex <- data.frame(base = c1, enet = c2, omaa = c3,
                     delta_enet = c2 - c1, delta_omaa = c3 - c1)

## ============================================================
## 8. KM quartiles (svykm; quartile factor)
## ============================================================
cat("\n========== 8. KM quartiles ==========\n")

cc$q <- factor(ntile(cc$acc_enet_sd, 4),
               levels = 1:4, labels = paste0("Q", 1:4))

des_q <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
  data = cc, nest = TRUE
)

km <- tryCatch(
  svykm(Surv(time, event) ~ q, design = des_q, se = FALSE),
  error = function(e) {
    warning("svykm failed: ", conditionMessage(e))
    NULL
  }
)

km_tab <- cc %>%
  group_by(q) %>%
  summarise(
    n          = n(),
    deaths     = sum(event, na.rm = TRUE),
    median_acc = median(acc_enet, na.rm = TRUE),
    .groups = "drop"
  )

cat("KM quartile table:\n")
print(km_tab, row.names = FALSE, digits = 3)

## ============================================================
## 9. Save
## ============================================================
out <- list(
  cox_tab   = cox_tab,
  cvd_row   = cvd_row,
  cindex    = cindex,
  cc        = cc,
  km_tab    = km_tab,
  km        = km,
  ph_p      = c(ENet = res_e$ph_p, OMAA = res_o$ph_p),
  n_cc      = nrow(cc),
  n_deaths  = sum(cc$event, na.rm = TRUE)
)

saveRDS(out, file.path(data_dir, "nhanes_cox.rds"))

write.csv(bind_rows(cox_tab, cvd_row),
          file.path(res_dir, "nhanes_cox_results.csv"),
          row.names = FALSE)

cat("\n========== module2_cox_mortality DONE ==========\n")
cat("Output:", file.path(data_dir, "nhanes_cox.rds"), "\n")
cat("Cox population:", nrow(cc), " | deaths:", sum(cc$event), "\n")
