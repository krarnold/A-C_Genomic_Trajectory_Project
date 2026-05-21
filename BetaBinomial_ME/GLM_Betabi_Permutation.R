# libraries
library(tidyr)
library(dplyr)
library(stringr)
library(glmmTMB)
library(broom.mixed)
library(furrr)
library(future)
library(purrr)

# ── SETTINGS ──────────────────────────────────────────────────────────────────
n_permutations <- 5       # number of permutations (following Kawecki et al.)
n_workers      <- 13      # number of cores for SNP-level parallelisation
# ──────────────────────────────────────────────────────────────────────────────

# ── PATHS ─────────────────────────────────────────────────────────────────────
# Set data_dir to the path of your local data directory
data_dir <- "."
out_file <- file.path(data_dir, "BetaBinomial_ME/Permuted_GLM_Betabi_Results.csv")
# ──────────────────────────────────────────────────────────────────────────────

# STEP 1: Load preprocessed long-format data saved by script 1
data_long <- read.csv(file.path(data_dir, "BetaBinomial_ME/Long_format_SNP_Data.csv"), header = TRUE) %>%
  mutate(
    treatment  = factor(treatment),
    population = factor(population),
    generation = as.numeric(generation)
  )

# get the 20 unique population identities
pop_treatments <- data_long %>%
  distinct(population, treatment)

# STEP 2: Split by SNP once — reused across all permutations
# This is the expensive structural operation and should not be repeated
message("Splitting data by SNP...")
data_split <- data_long %>%
  group_by(chr, pos) %>%
  group_split()
message("Split complete. ", length(data_split), " SNPs.")

# free the unsplit data from memory — no longer needed
rm(data_long)
gc()

# STEP 3: Define the GLMM fitting function
fit_glmm <- function(snp_data) {
  if (n_distinct(snp_data$treatment) < 2 ||
      n_distinct(snp_data$generation) < 2) return(NULL)

  tryCatch({
    model <- suppressMessages(glmmTMB(
      cbind(Minor_Count, Major_Count) ~ treatment * generation + (1 | population),
      data   = snp_data,
      family = betabinomial()
    ))

    broom.mixed::tidy(model, effects = "fixed") %>%
      filter(term %in% c("treatmentC>A", "generation", "treatmentC>A:generation")) %>%
      mutate(
        chr = unique(snp_data$chr),
        pos = unique(snp_data$pos)
      )
  }, error = function(e) NULL)
}

# STEP 4: Set up parallelisation across SNPs
plan(multisession, workers = n_workers)

# STEP 5: Run permutations sequentially
# For each permutation:
#   - reassign treatment labels within the already-split list (fast)
#   - fit GLMMs in parallel across SNPs
#   - write results to disk immediately to free memory
#   - append to output file so results are never all in memory at once

# write header to output file
write.csv(
  data.frame(perm = integer(), chr = character(), pos = integer(),
             term = character(), p.value = numeric()),
  out_file, row.names = FALSE
)

for (i in seq_len(n_permutations)) {

  message("Running permutation ", i, " of ", n_permutations, " ...")

  # shuffle treatment labels across the 20 populations
  shuffled <- pop_treatments %>%
    mutate(treatment = sample(treatment))

  # build a named lookup vector for fast replacement
  treatment_lookup <- setNames(
    as.character(shuffled$treatment),
    as.character(shuffled$population)
  )

  # reassign treatment labels within the already-split list
  # lapply is fast here because each element is a small data frame
  data_split_perm <- lapply(data_split, function(snp_data) {
    snp_data$treatment <- factor(treatment_lookup[as.character(snp_data$population)])
    snp_data
  })

  # fit GLMM across SNPs in parallel
  perm_result <- future_map_dfr(
    data_split_perm,
    fit_glmm,
    .progress = TRUE,
    .options  = furrr_options(seed = TRUE)
  ) %>%
    mutate(perm = i) %>%
    select(perm, chr, pos, term, p.value)

  # append to output file immediately and free memory
  write.table(
    perm_result,
    out_file,
    sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
  )

  rm(data_split_perm, perm_result)
  gc()

  message("Permutation ", i, " complete.")
}

message("All permutations complete. Results saved to ", out_file)
