# ============================================================================
# module2_preprocess.R — NHANES preprocessing
#
# Input:  <data_dir>/nhanes_omp/  (OMP genus table + XPT files + mortality .dat)
# Output: <data_dir>/nhanes_analysis.rds = list(meta, genus, flow)
#
# Steps: unclassified genera removed; complete-case audit; time_m > 0
# filter; mortality column-position diagnostics; robust UCOD coding;
# participant-flow accounting.
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

d <- file.path(data_dir, "nhanes_omp")
stopifnot(dir.exists(d))

## ============================================================
## 1. OMP genus table (drop unclassified genera)
## ============================================================
cat("\n========== 1. OMP genus table ==========\n")

gen <- read.delim(file.path(d, "dada2rb-genus-relative.txt"),
                  check.names = FALSE)
gen <- gen[!is.na(gen$SEQN) & gen$SEQN != "", ]
gen$SEQN <- as.numeric(as.character(gen$SEQN))

cat("Raw genus table:", nrow(gen), "samples x", ncol(gen) - 1, "genera\n")

# Drop unclassified genera (consistent with module 5)
genus_cols <- setdiff(colnames(gen), "SEQN")
bad_pattern <- "uncultured|unidentified|unclassified|unknown|^NA$|^$"
keep_genus <- !grepl(bad_pattern, genus_cols, ignore.case = TRUE)
cat("Unclassified genera dropped:", sum(!keep_genus), "\n")
cat("Dropped genera (first 10):\n")
print(head(genus_cols[!keep_genus], 10))

gen <- gen[, c("SEQN", genus_cols[keep_genus]), drop = FALSE]
cat("Genera retained:", ncol(gen) - 1, "\n")

## ============================================================
## 2. Demographics (two cycles)
## ============================================================
cat("\n========== 2. Demographics ==========\n")

d09 <- read_xpt(file.path(d, "DEMO_F.xpt"))
d11 <- read_xpt(file.path(d, "DEMO_G.xpt"))

vars <- c("SEQN", "RIDSTATR", "RIDAGEYR", "RIAGENDR",
          "RIDRETH1", "DMDEDUC2", "INDFMPIR",
          "WTMEC2YR", "SDMVPSU", "SDMVSTRA", "RIDEXPRG")

# Keep only columns present in both cycles (column names may differ
# slightly between NHANES releases)
vars_09 <- intersect(vars, colnames(d09))
vars_11 <- intersect(vars, colnames(d11))
vars_common <- intersect(vars_09, vars_11)

if (length(setdiff(vars, vars_common)) > 0) {
  warning("Missing columns: ", paste(setdiff(vars, vars_common), collapse = ", "))
}

demo <- bind_rows(d09[vars_common], d11[vars_common]) %>%
  filter(!is.na(SEQN))

# Cycle label determined directly from the XPT files (no SEQN hard-coding)
demo$cycle <- ifelse(demo$SEQN %in% d09$SEQN, "d09",
                     ifelse(demo$SEQN %in% d11$SEQN, "d11", NA))

cat("demo rows:", nrow(demo), "\n")
cat("Cycle distribution:\n"); print(table(demo$cycle, useNA = "ifany"))

## ============================================================
## 3. Covariates: BMI, smoking
## ============================================================
cat("\n========== 3. Covariates ==========\n")

bmx <- bind_rows(read_xpt(file.path(d, "BMX_F.xpt")),
                 read_xpt(file.path(d, "BMX_G.xpt"))) %>%
  select(SEQN, BMXBMI)

smq <- bind_rows(read_xpt(file.path(d, "SMQ_F.xpt")),
                 read_xpt(file.path(d, "SMQ_G.xpt"))) %>%
  select(SEQN, SMQ020, SMQ040) %>%
  mutate(smoke = case_when(
    SMQ020 == 2 ~ "never",
    SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "current",
    SMQ020 == 1 & SMQ040 == 3 ~ "former",
    TRUE ~ NA_character_
  ))

cat("BMX rows:", nrow(bmx), " | SMQ rows:", nrow(smq), "\n")
cat("Smoking distribution:\n"); print(table(smq$smoke, useNA = "ifany"))

## ============================================================
## 4. Mortality linkage (fixed-width files)
## ============================================================
cat("\n========== 4. Mortality linkage ==========\n")

# Column positions from the CDC 2019 public-release data dictionary:
#   SEQN        1-14
#   ELIGSTAT    15       (1 = eligible)
#   MORTSTAT    16       (0 = alive, 1 = dead)
#   UCOD        17-19    (113 cause-of-death codes)
#   PERMTH_INT  43-45    (months from interview to death/end of follow-up)
#   PERMTH_EXM  46-48    (months from MEC exam)
# Positions should be identical across years, but diagnose first.

mort_files <- c("NHANES_2009_2010_MORT_2019_PUBLIC.dat",
                "NHANES_2011_2012_MORT_2019_PUBLIC.dat")

parse_mort <- function(f) {
  l <- readLines(file.path(d, f))
  
  # Diagnostic: line count and character widths of first 3 lines
  cat("  ", f, " lines:", length(l),
      " first-3 line widths:", paste(nchar(l[1:3]), collapse = ","), "\n")
  
  data.frame(
    SEQN       = as.numeric(trimws(substr(l, 1, 14))),
    ELIGSTAT   = as.numeric(substr(l, 15, 15)),
    MORTSTAT   = suppressWarnings(as.numeric(substr(l, 16, 16))),
    UCOD       = trimws(substr(l, 17, 19)),
    PERMTH_INT = suppressWarnings(as.numeric(substr(l, 43, 45))),
    PERMTH_EXM = suppressWarnings(as.numeric(substr(l, 46, 48))),
    stringsAsFactors = FALSE
  )
}

mort <- bind_rows(lapply(mort_files, parse_mort))

# ---- Sanity checks ----
cat("\nMortality data sanity checks:\n")
cat("  Total rows:", nrow(mort), "\n")
cat("  ELIGSTAT distribution:\n"); print(table(mort$ELIGSTAT, useNA = "ifany"))
cat("  MORTSTAT distribution:\n"); print(table(mort$MORTSTAT, useNA = "ifany"))
cat("  UCOD non-empty fraction:", mean(nzchar(mort$UCOD)), "\n")
cat("  PERMTH_INT range:", paste(range(mort$PERMTH_INT, na.rm = TRUE),
                                 collapse = "-"), "\n")

stopifnot(all(mort$MORTSTAT[mort$ELIGSTAT != 1] %in% c(NA)))
stopifnot(all(mort$UCOD[mort$MORTSTAT == 0 & !is.na(mort$MORTSTAT)] == ""))

# ---- Robustness: follow-up time must be positive ----
mort <- mort %>%
  mutate(
    # PERMTH_INT is measured from the interview date; used as the main
    # time axis
    time_m = ifelse(!is.na(PERMTH_INT) & PERMTH_INT > 0,
                    PERMTH_INT, NA),
    # CVD death: 001 = Diseases of heart, 005 = Cerebrovascular diseases
    # (113 cause-of-death codes, 2019 public-release data dictionary)
    cvd_death = as.numeric(MORTSTAT == 1 & UCOD %in% c("001", "005"))
  )

cat("  time_m > 0 fraction:", mean(!is.na(mort$time_m)), "\n")
cat("  CVD deaths:", sum(mort$cvd_death, na.rm = TRUE), "\n")

## ============================================================
## 5. Merge & participant flow
## ============================================================
cat("\n========== 5. Merge & Flow ==========\n")

flow <- c(genus_table = nrow(gen))

dat <- gen %>% select(SEQN) %>% inner_join(demo, by = "SEQN")
flow["with_demographics"] <- nrow(dat)

dat <- dat %>% inner_join(mort, by = "SEQN")
flow["with_mortality"] <- nrow(dat)

dat <- dat %>% filter(RIDAGEYR >= 18)
flow["age_ge18"] <- nrow(dat)

dat <- dat %>% filter(ELIGSTAT == 1)
flow["mort_eligible"] <- nrow(dat)

# time_m > 0 (required by Cox)
dat <- dat %>% filter(!is.na(time_m) & time_m > 0)
flow["time_m_positive"] <- nrow(dat)

# Left-join covariates
dat <- dat %>%
  left_join(bmx, by = "SEQN") %>%
  left_join(smq[, c("SEQN", "smoke")], by = "SEQN") %>%
  mutate(DMDEDUC2 = ifelse(DMDEDUC2 %in% 1:5, DMDEDUC2, NA))

# ---- Complete cases (all covariates) ----
cc_vars <- c("RIDAGEYR", "RIAGENDR", "RIDRETH1", "DMDEDUC2",
             "INDFMPIR", "BMXBMI", "smoke")
dat$cc_complete <- complete.cases(dat[cc_vars])

flow["complete_case_full"] <- sum(dat$cc_complete)
flow["deaths_cc"] <- sum(dat$MORTSTAT == 1 & dat$cc_complete, na.rm = TRUE)

cat("\nAnalysis population flow:\n"); print(flow)

cat("\n--- Key numbers ---\n")
cat("Age >= 18 and mortality-linkable:", flow["mort_eligible"], "\n")
cat("time_m > 0:", flow["time_m_positive"], "\n")
cat("Complete cases:", flow["complete_case_full"], "\n")
cat("Deaths among complete cases:", flow["deaths_cc"], "\n")
cat("CVD deaths (complete cases):",
    sum(dat$cvd_death == 1 & dat$cc_complete, na.rm = TRUE), "\n")

cat("Median follow-up (alive):",
    median(dat$time_m[dat$MORTSTAT == 0], na.rm = TRUE), "months\n")

# ---- Sample-size warnings ----
if (flow["complete_case_full"] < 1500) {
  warning("Complete cases < 1500; Cox power may be insufficient")
}
if (flow["deaths_cc"] < 100) {
  warning("Deaths among complete cases < 100; Cox power may be insufficient")
}

## ============================================================
## 6. Feature matrix (5% prevalence filter)
## ============================================================
cat("\n========== 6. Feature matrix ==========\n")

gen2 <- gen %>% filter(SEQN %in% dat$SEQN)
mat <- as.matrix(gen2[, -1])
rownames(mat) <- gen2$SEQN

cat("Before filtering:", nrow(mat), "samples x", ncol(mat), "genera\n")

# Drop all-zero samples (no valid reads; CLR would produce NaN and
# downstream randomForest rejects them)
zero_rows <- rowSums(mat, na.rm = TRUE) == 0
if (any(zero_rows)) {
  cat("All-zero samples dropped:", sum(zero_rows), "(SEQN:",
      paste(head(rownames(mat)[zero_rows], 10), collapse = ", "), ")\n")
  mat <- mat[!zero_rows, , drop = FALSE]
  dat <- dat[dat$SEQN %in% as.numeric(rownames(mat)), ]
}

# 5% prevalence filter (as in the original OMAA publication)
keep <- colMeans(mat > 0) >= 0.05
mat <- mat[, keep, drop = FALSE]

cat("After 5% prevalence filter:", ncol(mat), "genera\n")
cat("Max relative abundance:", max(mat, na.rm = TRUE), "\n")

# ---- Sanity ----
stopifnot(all(mat >= 0, na.rm = TRUE))
stopifnot(max(mat, na.rm = TRUE) < 1.5)   # values near 1.0 allowed, not >1.5

## ============================================================
## 7. Save
## ============================================================
saveRDS(list(meta = dat, genus = mat, flow = flow),
        file.path(data_dir, "nhanes_analysis.rds"))

cat("\n========== module2_preprocess DONE ==========\n")
cat("Output:", file.path(data_dir, "nhanes_analysis.rds"), "\n")
cat("Analysis population:", nrow(dat), " | complete cases:", sum(dat$cc_complete), "\n")
cat("Genus features:", ncol(mat), "\n")
cat("Deaths (complete cases):", sum(dat$MORTSTAT == 1 & dat$cc_complete, na.rm = TRUE), "\n")
