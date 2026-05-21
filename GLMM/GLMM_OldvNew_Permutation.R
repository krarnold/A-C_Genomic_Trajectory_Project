# Permutation-based null distribution for the Old vs New comparisons
# Uses quasibinomial GLM matching the observed-data approach
# Following Wiberg et al. (2017) model structure (2): y ~ treatment

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
n_permutations <- 5
n_workers      <- 13
# ──────────────────────────────────────────────────────────────────────────────

# STEP 1: Read only the header to identify available columns
message("Reading column names from SNP table...")
all_col_names <- colnames(read.csv(snp_path, header = TRUE, nrows = 0))

# STEP 2: Quasibinomial GLM fitting function
fit_glm <- function(snp_data) {
  tryCatch({
    fit <- suppressWarnings(glm(
      cbind(A1, A2) ~ treatment,
      data   = snp_data,
      family = quasibinomial()
    ))
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

# STEP 3: Function to run all permutations for a given comparison
run_comparison_perms <- function(comparison_name, col_pattern, out_file) {

  message("=== Starting ", comparison_name, " permutations ===")

  data_cols   <- grep(col_pattern, all_col_names, value = TRUE)
  needed_cols <- c("chr", "pos", data_cols)

  col_classes <- rep("NULL", length(all_col_names))
  col_classes[all_col_names %in% needed_cols] <- NA

  message("  Reading ", length(needed_cols), " columns from SNP table...")
  sub <- read.csv(snp_path, header = TRUE, colClasses = col_classes)

  sum_maf <- rowSums(sub[, seq(3, ncol(sub), by = 2)] / sub[, seq(4, ncol(sub), by = 2)])
  sub     <- cbind(sub, sum_maf)
  sub     <- subset(sub, sum_maf > 0 & sum_maf < 20)
  sub$sum_maf <- NULL

  coverage <- sub[, grep("_total_coverage$", names(sub), value = TRUE)]
  minor    <- sub[, grep("_minor_count$",    names(sub), value = TRUE)]
  major    <- coverage - minor

  true_treatment <- c(rep("old", 10), rep("new", 10))

  message("  SNPs after polymorphism filter: ", nrow(minor))
  message("  Building SNP data list once (reused across all permutations)...")

  # build one SNP data frame per SNP with TRUE treatment labels
  # we overwrite the treatment column inside each permutation
  snp_split <- lapply(seq_len(nrow(minor)), function(k) {
    data.frame(
      chr       = sub$chr[k],
      pos       = sub$pos[k],
      A1        = as.numeric(minor[k, ]),
      A2        = as.numeric(major[k, ]),
      treatment = true_treatment
    )
  })

  # free intermediates — no longer needed now that snp_split is built
  rm(sub, coverage, minor, major)
  gc()

  message("  SNP list built. Starting permutations.")

  # write header to output file
  write.table(
    data.frame(perm = integer(), chr = character(), pos = integer(), pval = numeric()),
    out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )

  for (i in seq_len(n_permutations)) {

    message("  Running permutation ", i, " of ", n_permutations, "...")

    # shuffle treatment labels across the 20 populations
    shuffled_treatment <- sample(true_treatment)

    # update the treatment column within each SNP data frame
    snp_split_perm <- lapply(snp_split, function(snp_df) {
      snp_df$treatment <- shuffled_treatment
      snp_df
    })

    # fit GLM across SNPs in parallel
    perm_result <- future_map_dfr(snp_split_perm, fit_glm, .progress = TRUE) %>%
      mutate(perm = i) %>%
      select(perm, chr, pos, pval)

    write.table(
      perm_result,
      out_file,
      sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = FALSE, append = TRUE
    )

    rm(snp_split_perm, perm_result)
    gc()

    message("  Permutation ", i, " complete.")
  }

  # free the SNP split before returning
  rm(snp_split)
  gc()

  message("=== ", comparison_name, " complete — saved to ", out_file, " ===")
}

# STEP 4: Set up parallelization and run both comparisons
plan(multisession, workers = n_workers)

run_comparison_perms(
  comparison_name = "AvA",
  col_pattern     = "^(ACO|AO|NACO|ANCO)[0-9]*_24(_minor_count|_total_coverage)$",
  out_file        = file.path(out_dir, "GLMM_AvA_Permuted.txt")
)

run_comparison_perms(
  comparison_name = "CvC",
  col_pattern     = "^(CO|NCO|CACO|CAO)[0-9]*_24(_minor_count|_total_coverage)$",
  out_file        = file.path(out_dir, "GLMM_CvC_Permuted.txt")
)

message("All permutations complete.")
