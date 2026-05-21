# Observed quasibinomial GLM fits for Old vs New comparisons
# Following Wiberg et al. (2017) model structure (2): y ~ treatment
# AvA: founder A (ACO, AO) vs reversed A (NACO, ANCO)
# CvC: founder C (CO, NCO) vs reversed C (CACO, CAO)
#
# Model: cbind(minor, major) ~ treatment, family = quasibinomial
# Significance: t-test on the treatment coefficient from summary(fit)
#
# Reads only the columns needed for each comparison to minimize memory use

# libraries
library(dplyr)
library(furrr)
library(future)
library(purrr)

# ── PATHS ─────────────────────────────────────────────────────────────────────
# Set data_dir to the path of your local data directory
data_dir <- "."
snp_path <- file.path(data_dir, "SNP_Table/AC_trajectory_modified_snp_table.csv")
out_dir  <- file.path(data_dir, "GLMM")
# ──────────────────────────────────────────────────────────────────────────────

# ── SETTINGS ──────────────────────────────────────────────────────────────────
n_workers <- 13
# ──────────────────────────────────────────────────────────────────────────────

# STEP 1: Read only the header to identify available columns
message("Reading column names from SNP table...")
all_col_names <- colnames(read.csv(snp_path, header = TRUE, nrows = 0))

# STEP 2: Quasibinomial GLM fitting function
# following Wiberg et al. (2017) model structure (2)
fit_glm <- function(snp_data) {
  tryCatch({
    fit <- suppressWarnings(glm(
      cbind(A1, A2) ~ treatment,
      data   = snp_data,
      family = quasibinomial()
    ))
    # second row is the treatment coefficient; "Pr(>|t|)" for quasibinomial
    pval <- summary(fit)$coefficients[2, "Pr(>|t|)"]
    tibble(
      chr  = unique(snp_data$chr),
      pos  = unique(snp_data$pos),
      pval = pval
    )
  }, error = function(e) {
    tibble(
      chr  = unique(snp_data$chr),
      pos  = unique(snp_data$pos),
      pval = NA_real_
    )
  })
}

# STEP 3: Function to run one comparison
run_comparison <- function(comparison_name, col_pattern, out_file) {

  message("=== Starting ", comparison_name, " ===")

  data_cols   <- grep(col_pattern, all_col_names, value = TRUE)
  needed_cols <- c("chr", "pos", data_cols)

  col_classes <- rep("NULL", length(all_col_names))
  col_classes[all_col_names %in% needed_cols] <- NA

  message("  Reading ", length(needed_cols), " columns from SNP table...")
  sub <- read.csv(snp_path, header = TRUE, colClasses = col_classes)

  # polymorphism filter
  sum_maf <- rowSums(sub[, seq(3, ncol(sub), by = 2)] / sub[, seq(4, ncol(sub), by = 2)])
  sub     <- cbind(sub, sum_maf)
  sub     <- subset(sub, sum_maf > 0 & sum_maf < 20)
  sub$sum_maf <- NULL

  coverage <- sub[, grep("_total_coverage$", names(sub), value = TRUE)]
  minor    <- sub[, grep("_minor_count$",    names(sub), value = TRUE)]
  major    <- coverage - minor

  # first 10 cols are "old" (founders), last 10 are "new" (reversed)
  treatment_vec <- c(rep("old", 10), rep("new", 10))

  message("  SNPs after polymorphism filter: ", nrow(minor))
  message("  Building SNP data list...")

  snp_split <- lapply(seq_len(nrow(minor)), function(i) {
    data.frame(
      chr       = sub$chr[i],
      pos       = sub$pos[i],
      A1        = as.numeric(minor[i, ]),
      A2        = as.numeric(major[i, ]),
      treatment = treatment_vec
    )
  })

  # free intermediates — no longer needed now that snp_split is built
  rm(sub, coverage, minor, major)
  gc()

  message("  Starting parallel model fitting...")

  results <- future_map_dfr(snp_split, fit_glm, .progress = TRUE)

  rm(snp_split)
  gc()

  write.table(results, out_file, quote = FALSE, row.names = FALSE, sep = "\t")

  message("=== ", comparison_name, " complete — saved to ", out_file, " ===")
}

# STEP 4: Set up parallelization and run both comparisons
plan(multisession, workers = n_workers)

run_comparison(
  comparison_name = "AvA",
  col_pattern     = "^(ACO|AO|NACO|ANCO)[0-9]*_24(_minor_count|_total_coverage)$",
  out_file        = file.path(out_dir, "GLMM_AvA.txt")
)

run_comparison(
  comparison_name = "CvC",
  col_pattern     = "^(CO|NCO|CACO|CAO)[0-9]*_24(_minor_count|_total_coverage)$",
  out_file        = file.path(out_dir, "GLMM_CvC.txt")
)

message("All observed GLM fits complete.")
