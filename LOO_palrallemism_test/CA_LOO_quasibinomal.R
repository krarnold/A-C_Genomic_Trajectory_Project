# ===============================================================
#  C>A Leave-One-Out Parallelism Analysis
#  (FDR < 0.05 and |Δp| >= 0.02, random controls)
#
#  Inputs:
#     Long_format_SNP_Data.csv (SNP counts)
#  Outputs:
#     CA_LOO_medians.csv  (median Δp for targets vs random controls)
# ===============================================================

# --- Libraries ---
library(dplyr)
library(tidyr)
library(stringr)
library(broom)
library(furrr)
library(purrr)

# ---------------------------------------------------------------
# 1. Read and prepare C>A data
# ---------------------------------------------------------------
data_long <- read.csv("Long_format_SNP_Data.csv", header = TRUE) %>%
  filter(treatment == "C>A", chr != "4") %>%  # exclude chr 4
  mutate(
    replicate  = str_extract(population, "C>A\\d+"),
    generation = as.numeric(generation)
  )

replicates <- sort(unique(data_long$replicate))

# ---------------------------------------------------------------
# 2. Helper functions
# ---------------------------------------------------------------

# Fit quasibinomial GLM for allele frequency change over generations
fit_glm_gen <- function(snp_data) {
  if (n_distinct(snp_data$generation) < 2) return(NULL)
  tryCatch({
    model <- glm(
      cbind(Minor_Count, Major_Count) ~ generation + replicate,
      data   = snp_data,
      family = quasibinomial()
    )
    tidy(model) %>%
      filter(term == "generation") %>%
      mutate(
        chr = unique(snp_data$chr),
        pos = unique(snp_data$pos)
      )
  }, error = function(e) NULL)
}

# Compute allele frequency change (Δp)
compute_deltas <- function(df) {
  df %>%
    group_by(chr, pos, generation) %>%
    summarise(
      p = sum(Minor_Count) / sum(Minor_Count + Major_Count),
      .groups = "drop_last"
    ) %>%
    summarise(
      delta_p = p[which.max(generation)] - p[which.min(generation)],
      .groups = "drop"
    )
}

# ---------------------------------------------------------------
# 3. Parallel setup
# ---------------------------------------------------------------
plan(multisession, workers = 10)

# ---------------------------------------------------------------
# 4. Leave-One-Out analysis (FDR + effect-size filter, random controls)
# ---------------------------------------------------------------
median_list <- list()

for (rep in replicates) {
  message("=== Leaving out replicate ", rep, " ===")
  
  # Split data into training and test sets
  train_data <- data_long %>% filter(replicate != rep)
  test_data  <- data_long %>% filter(replicate == rep)
  
  # Run per-SNP GLMs on training data
  snp_split <- train_data %>% group_by(chr, pos) %>% group_split()
  
  results_train <- future_map_dfr(
    snp_split,
    fit_glm_gen,
    .progress = TRUE,
    .options  = furrr_options(seed = TRUE)
  )
  
  rm(snp_split); gc()
  
  # Compute Δp in training data
  delta_train <- compute_deltas(train_data)
  
  # Merge GLM results with Δp estimates
  results_train <- results_train %>%
    left_join(delta_train, by = c("chr", "pos")) %>%
    filter(!is.na(delta_p), !is.na(p.value)) %>%
    mutate(FDR = p.adjust(p.value, method = "fdr"))
  
  # Define significant (target) and nonsignificant (control) SNPs
  targets <- results_train %>% filter(FDR < 0.05, abs(delta_p) >= 0.02)
  nonsig  <- results_train %>% filter(!(FDR < 0.05 & abs(delta_p) >= 0.02))
  
  if (nrow(targets) == 0 | nrow(nonsig) == 0) {
    message("No targets or controls for replicate ", rep, " — skipping.")
    next
  }
  
  # Randomly sample equal number of controls
  matched_controls <- nonsig %>%
    slice_sample(n = min(nrow(targets), nrow(nonsig)))
  
  sig_table <- bind_rows(
    targets %>% mutate(sig = "target"),
    matched_controls %>% mutate(sig = "control")
  ) %>%
    select(chr, pos, sig) %>%
    distinct()
  
  # Compute Δp in left-out test replicate
  delta_test <- compute_deltas(test_data)
  
  merged <- delta_test %>% inner_join(sig_table, by = c("chr", "pos"))
  
  if (nrow(merged) == 0) {
    message("No overlapping SNPs in ", rep, " — skipping.")
    next
  }
  
  # Compute medians by chromosome
  med_chr <- merged %>%
    group_by(chr, sig) %>%
    summarise(
      median_delta = median(delta_p, na.rm = TRUE),
      n_snps       = n(),
      .groups      = "drop"
    ) %>%
    mutate(replicate = rep)
  
  # Genome-wide median
  med_genome <- merged %>%
    group_by(sig) %>%
    summarise(
      median_delta = median(delta_p, na.rm = TRUE),
      n_snps       = n(),
      .groups      = "drop"
    ) %>%
    mutate(replicate = rep, chr = "genome")
  
  median_list[[rep]] <- bind_rows(med_genome, med_chr)
  
  rm(train_data, test_data, results_train, merged,
     delta_train, targets, nonsig, matched_controls,
     sig_table, delta_test)
  gc()
}

# ---------------------------------------------------------------
# 5. Combine and write results
# ---------------------------------------------------------------
if (length(median_list) == 0) {
  stop("No LOO results generated — check thresholds or input data.")
}

median_table <- bind_rows(median_list)
write.csv(median_table, "CA_LOO_medians.csv", row.names = FALSE)
message("Finished C>A LOO analysis (FDR<0.05, |Δp|>=0.02, chr 4 excluded, random controls). Wrote CA_LOO_medians.csv")
