# ============================================================================
# module3_parse.R — Parse HMP2 (IBDMDB) MetaPhlAn2 taxonomic profiles
#
# Input (two modes, auto-detected):
#   Mode A (recommended): <data_dir>/hmp2_taxonomic_profiles.tsv.gz
#     — official IBDMDB merged table (1,638 MGX samples):
#     https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/2018-05-04/taxonomic_profiles.tsv.gz
#   Mode B (fallback): <data_dir>/hmp2_tax/*_taxonomic_profile.tsv (per-file)
#
# Output:
#   <data_dir>/hmp2_species_matrix.rds (species x samples, relative abundance)
#   <data_dir>/hmp2_sample_meta.rds    (sample metadata)
#
# Key parsing rules: strain rows (t__) are EXCLUDED — in MetaPhlAn2
# hierarchical output strain abundance is already contained in the species
# row, so keeping them would double-count abundance and pollute species
# names.
# ============================================================================

source(file.path(Sys.getenv("OG_CLOCK_ROOT", unset = "C:/workspace"),
                 "module0_setup.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
})

## ============================================================
## 0. Detect input mode
## ============================================================
merged_file <- file.path(data_dir, "hmp2_taxonomic_profiles.tsv.gz")
if (!file.exists(merged_file)) {
  alt <- file.path(data_dir, "hmp2_taxonomic_profiles.tsv")
  if (file.exists(alt)) merged_file <- alt else merged_file <- NULL
}

if (is.null(merged_file) &&
    !dir.exists(file.path(data_dir, "hmp2_tax"))) {
  stop("No HMP2 profile input found. Place one of the following in data_dir:
",
       "  (1) hmp2_taxonomic_profiles.tsv.gz (official merged table, recommended)
",
       "  (2) or an hmp2_tax/ directory (per-file .tsv)
",
       "Merged-table download:
",
       "  download.file('https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/2018-05-04/taxonomic_profiles.tsv.gz',
",
       "               file.path(data_dir, 'hmp2_taxonomic_profiles.tsv.gz'), mode = 'wb')")
}

## ============================================================
## 1. Parse -> per-sample species abundance lists
## ============================================================
cat("
========== 1. Parse profiles ==========
")

# Lineage rows -> species name + abundance (strain rows excluded)
clade_to_species <- function(cl, v) {
  # s__ rows without t__ (a "...|s__Sp|t__XYZ" row is a strain whose
  # abundance is already included in the s__Sp row; keeping it would
  # double-count)
  is_sp <- grepl("(\\||^)s__", cl) & !grepl("\\|t__", cl)
  if (!any(is_sp)) return(NULL)
  sp <- sub(".*\\|s__", "", cl[is_sp])
  sp <- sub("^s__", "", sp)
  out <- tapply(v[is_sp], sp, sum)   # aggregate duplicate species names within a sample
  out[out > 0]
}

if (!is.null(merged_file)) {
  cat("Mode A: official merged table", basename(merged_file), "
")
  x <- read.delim(merged_file, check.names = FALSE, stringsAsFactors = FALSE)
  sid_all <- colnames(x)[-1]
  cl <- as.character(x[[1]])
  lists <- lapply(seq_along(sid_all), function(j) {
    v <- as.numeric(x[[j + 1]])
    clade_to_species(cl, v)
  })
  names(lists) <- sid_all
  sample_ids <- sid_all
  ok <- !sapply(lists, is.null)
  cat("Samples:", length(sample_ids), " | parsed:", sum(ok), "
")
} else {
  cat("Mode B: per-file directory
")
  d <- file.path(data_dir, "hmp2_tax")
  files <- list.files(d, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE)
  cat("files:", length(files), "
")
  if (length(files) == 0) {
    # Diagnostics: does the directory exist, and what does it contain?
    cat("Diagnostics: directory exists?", dir.exists(d), "-> ", d, "
")
    if (dir.exists(d)) {
      all_f <- list.files(d, recursive = TRUE)
      cat("Directory contains", length(all_f), "files; first 10:
")
      print(head(all_f, 10))
    } else {
      hits <- list.files(data_dir, pattern = "taxonomic", recursive = TRUE)
      cat("Files/dirs matching 'taxonomic' under data_dir:
")
      print(head(hits, 10))
    }
    stop("No .tsv/.tsv.gz files found. Recommended: download the official merged table (666 KB, 1,638 samples):
",
         "  download.file('https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/2018-05-04/taxonomic_profiles.tsv.gz',
",
         "               file.path(data_dir, 'hmp2_taxonomic_profiles.tsv.gz'), mode = 'wb')
",
         "Or check the hmp2_tax directory path.")
  }
  
  sid_of <- function(f) sub("_taxonomic_profile\\.tsv(\\.gz)?$", "", basename(f))
  sample_ids <- sid_of(files)
  stopifnot(!any(duplicated(sample_ids)))
  
  parse_one <- function(f) {
    x <- tryCatch(read.delim(f, check.names = FALSE), error = function(e) NULL)
    if (is.null(x) || ncol(x) < 2) return(NULL)
    cl <- as.character(x[[1]])
    v  <- as.numeric(x[[ncol(x)]])
    clade_to_species(cl, v)
  }
  
  n_cores <- max(1, parallel::detectCores() - 1)
  cat("Using", n_cores, "cores
")
  t0 <- Sys.time()
  lists <- parallel::mclapply(files, parse_one, mc.cores = n_cores)
  cat("Parsing time:",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min
")
  ok <- !sapply(lists, is.null)
  cat("Parsed:", sum(ok), " | failed:", sum(!ok), "
")
  if (any(!ok)) {
    cat("First 10 failed files:
")
    print(head(basename(files[!ok]), 10))
  }
}

## ============================================================
## 2. Build species matrix
## ============================================================
cat("
========== 2. Build matrix ==========
")

all_sp <- sort(unique(unlist(lapply(lists[ok], names))))
cat("Species:", length(all_sp), "
")

mat <- matrix(0,
              nrow = length(all_sp),
              ncol = sum(ok),
              dimnames = list(all_sp, sample_ids[ok]))

for (i in which(ok)) {
  sp_i <- names(lists[[i]])
  mat[sp_i, sample_ids[i]] <- lists[[i]]
}

## ============================================================
## 3. Convert to relative abundance (guard against zero colSums)
## ============================================================
cat("
========== 3. Relative abundance ==========
")

cs <- colSums(mat)
cat("colSums range:", paste(round(range(cs), 3), collapse = " - "), "
")

if (any(cs <= 0)) {
  cat("Warning:", sum(cs <= 0), "samples with zero total abundance, dropped
")
  keep <- cs > 0
  mat <- mat[, keep, drop = FALSE]
  cs  <- cs[keep]
}

mat <- sweep(mat, 2, cs, "/")
cat("Matrix:", paste(dim(mat), collapse = " x "), "
")
cat("colSums check (should be ~1):",
    paste(round(range(colSums(mat)), 3), collapse = " - "), "
")

stopifnot(all(abs(colSums(mat) - 1) < 0.05))

saveRDS(mat, file.path(data_dir, "hmp2_species_matrix.rds"))

## ============================================================
## 4. Metadata merge (uniqueness check)
## ============================================================
cat("
========== 4. Metadata ==========
")

m <- read.csv(file.path(data_dir, "hmp2_metadata.csv"),
              check.names = FALSE, stringsAsFactors = FALSE)

# Metagenomics rows only
mg <- m[m[["data_type"]] == "metagenomics", ]
mg$sid <- sub("_P$", "", mg[["External ID"]])

cat("mg rows:", nrow(mg), "
")
cat("mg$sid unique:", !any(duplicated(mg$sid)), "
")

if (any(duplicated(mg$sid))) {
  dup_sids <- mg$sid[duplicated(mg$sid)]
  cat("First 5 duplicated sids:", paste(head(unique(dup_sids), 5), collapse = ", "), "
")
  mg <- mg[!duplicated(mg$sid), ]   # keep first row per sid
  cat("mg rows after dedup:", nrow(mg), "
")
}

# Preserve the column order of mat
idx <- match(colnames(mat), mg$sid)
md <- mg[idx, c("sid", "Participant ID", "week_num", "interval_days",
                "consent_age", "diagnosis", "sex"), drop = FALSE]

cat("Metadata matched:", sum(!is.na(md[["Participant ID"]])), "/", ncol(mat), "
")
cat("Subjects:", length(unique(na.omit(md[["Participant ID"]]))), "
")
cat("With consent_age:", sum(!is.na(md$consent_age)), "
")

if (any(is.na(md[["Participant ID"]]))) {
  miss_idx <- which(is.na(md[["Participant ID"]]))
  cat("First 5 samples with missing metadata:",
      paste(colnames(mat)[miss_idx[1:min(5, length(miss_idx))]],
            collapse = ", "), "
")
}

saveRDS(md, file.path(data_dir, "hmp2_sample_meta.rds"))

## ============================================================
## 5. Summary
## ============================================================
cat("
========== module3_parse DONE ==========
")
cat("Species matrix:", paste(dim(mat), collapse = " x "), "
")
cat("Subjects:", length(unique(na.omit(md[["Participant ID"]]))), "
")
cat("Diagnosis distribution:
"); print(table(md$diagnosis, useNA = "ifany"))
cat("Outputs:
")
cat("  ", file.path(data_dir, "hmp2_species_matrix.rds"), "
")
cat("  ", file.path(data_dir, "hmp2_sample_meta.rds"), "
")
