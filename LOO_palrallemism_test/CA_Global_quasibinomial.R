# ===============================================================
#  C>A Global Parallelism Scan (Bitter-style, quasibinomial)
#  Output:
#     1) CA_results_all.csv  (per-SNP GLM results, chr 4 excluded)
# ===============================================================

setwd("~/Library/CloudStorage/Dropbox/AC Trajectory Genomic Data/LOO_palrallemism_test")
set.seed(123)

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
  filter(treatment == "C>A", chr != "4") %>%
  mutate(
    replicate = str_extract(population, "C>A\\d+"),
    generation = as.numeric(generation),
    population = factor(population)
  )

# ---------------------------------------------------------------
# 2. Define GLM function
# ---------------------------------------------------------------
fit_glm_gen <- function(snp_data) {
  if (n_distinct(snp_data$generation) < 2) return(NULL)
  tryCatch({
    model <- glm(
      cbind(Minor_Count, Major_Count) ~ generation + replicate,
      data = snp_data,
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

# ---------------------------------------------------------------
# 3. Parallel setup
# ---------------------------------------------------------------
plan(multisession, workers = 10)

# ---------------------------------------------------------------
# 4. Run GLMs across all SNPs (chr != 4)
# ---------------------------------------------------------------
message("Running global C>A per-SNP GLMs (excluding chr 4)...")

snp_split_all <- data_long %>%
  group_by(chr, pos) %>%
  group_split()

results_all <- future_map_dfr(
  snp_split_all,
  fit_glm_gen,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

rm(snp_split_all); gc()

# ---------------------------------------------------------------
# 5. Write full per-SNP results
# ---------------------------------------------------------------
write.csv(results_all, "CA_results_all.csv", row.names = FALSE)
message(paste("Finished C>A per-SNP scan. Wrote", nrow(results_all), "rows to CA_results_all.csv"))

