library(dplyr)
library(tidyr)
library(stringr)

# ============================================================
# GENOME-WIDE HETEROZYGOSITY SUMMARY STATISTICS SCRIPT
# Computes H = 2f(1-f) for all population groups across
# all timepoints from the raw SNP table
# Covers:
#   - A-type founders (CACO/CAO at t18, ACO/AO at t24)
#   - C-type founders (NACO/ANCO at t18, CO/NCO at t24)
#   - A->C trajectory (CACO/CAO at t18 as gen 0, t19, t20, t24)
#   - C->A trajectory (NACO/ANCO at t18 as gen 0, t19, t20, t24)
# ============================================================

# ---- PATHS ----
# Set data_dir to the path of your local data directory
data_dir       <- "."
snp_table_path <- file.path(data_dir, "SNP_Table/arnold25_genomic_trajectory_SNP_table_counts.csv")
output_dir     <- file.path(data_dir, "Heterozygosity")

# ---- LOAD DATA ----
cat("Reading raw SNP table...\n")
raw <- read.csv(snp_table_path, header = TRUE)
cat("Total sites:", nrow(raw), "\n")
cat("Total columns:", ncol(raw), "\n\n")

# ---- DEFINE POPULATION GROUPS ----
sample_map <- tribble(
  ~prefix,  ~year, ~group,            ~generation,

  # A-type founders sampled at t18
  "CACO1",  "18",  "A-type Founders", 0,
  "CACO2",  "18",  "A-type Founders", 0,
  "CACO3",  "18",  "A-type Founders", 0,
  "CACO4",  "18",  "A-type Founders", 0,
  "CACO5",  "18",  "A-type Founders", 0,
  "CAO1",   "18",  "A-type Founders", 0,
  "CAO2",   "18",  "A-type Founders", 0,
  "CAO3",   "18",  "A-type Founders", 0,
  "CAO4",   "18",  "A-type Founders", 0,
  "CAO5",   "18",  "A-type Founders", 0,

  # A-type founders sampled at t24
  "ACO1",   "24",  "A-type Founders", 219,
  "ACO2",   "24",  "A-type Founders", 219,
  "ACO3",   "24",  "A-type Founders", 219,
  "ACO4",   "24",  "A-type Founders", 219,
  "ACO5",   "24",  "A-type Founders", 219,
  "AO1",    "24",  "A-type Founders", 219,
  "AO2",    "24",  "A-type Founders", 219,
  "AO3",    "24",  "A-type Founders", 219,
  "AO4",    "24",  "A-type Founders", 219,
  "AO5",    "24",  "A-type Founders", 219,

  # C-type founders sampled at t18
  "NACO1",  "18",  "C-type Founders", 0,
  "NACO2",  "18",  "C-type Founders", 0,
  "NACO3",  "18",  "C-type Founders", 0,
  "NACO4",  "18",  "C-type Founders", 0,
  "NACO5",  "18",  "C-type Founders", 0,
  "ANCO1",  "18",  "C-type Founders", 0,
  "ANCO2",  "18",  "C-type Founders", 0,
  "ANCO3",  "18",  "C-type Founders", 0,
  "ANCO4",  "18",  "C-type Founders", 0,
  "ANCO5",  "18",  "C-type Founders", 0,

  # C-type founders sampled at t24
  "CO1",    "24",  "C-type Founders", 219,
  "CO2",    "24",  "C-type Founders", 219,
  "CO3",    "24",  "C-type Founders", 219,
  "CO4",    "24",  "C-type Founders", 219,
  "CO5",    "24",  "C-type Founders", 219,
  "NCO1",   "24",  "C-type Founders", 219,
  "NCO2",   "24",  "C-type Founders", 219,
  "NCO3",   "24",  "C-type Founders", 219,
  "NCO4",   "24",  "C-type Founders", 219,
  "NCO5",   "24",  "C-type Founders", 219,

  # A->C trajectory baseline at t18 (generation 0)
  "CACO1",  "18",  "A->C Trajectory", 0,
  "CACO2",  "18",  "A->C Trajectory", 0,
  "CACO3",  "18",  "A->C Trajectory", 0,
  "CACO4",  "18",  "A->C Trajectory", 0,
  "CACO5",  "18",  "A->C Trajectory", 0,
  "CAO1",   "18",  "A->C Trajectory", 0,
  "CAO2",   "18",  "A->C Trajectory", 0,
  "CAO3",   "18",  "A->C Trajectory", 0,
  "CAO4",   "18",  "A->C Trajectory", 0,
  "CAO5",   "18",  "A->C Trajectory", 0,

  # A->C trajectory at t19
  "CACO1",  "19",  "A->C Trajectory", 6,
  "CACO2",  "19",  "A->C Trajectory", 6,
  "CACO3",  "19",  "A->C Trajectory", 6,
  "CACO4",  "19",  "A->C Trajectory", 6,
  "CACO5",  "19",  "A->C Trajectory", 6,
  "CAO1",   "19",  "A->C Trajectory", 6,
  "CAO2",   "19",  "A->C Trajectory", 6,
  "CAO3",   "19",  "A->C Trajectory", 6,
  "CAO4",   "19",  "A->C Trajectory", 6,
  "CAO5",   "19",  "A->C Trajectory", 6,

  # A->C trajectory at t20
  "CACO1",  "20",  "A->C Trajectory", 10,
  "CACO2",  "20",  "A->C Trajectory", 10,
  "CACO3",  "20",  "A->C Trajectory", 10,
  "CACO4",  "20",  "A->C Trajectory", 10,
  "CACO5",  "20",  "A->C Trajectory", 10,
  "CAO1",   "20",  "A->C Trajectory", 10,
  "CAO2",   "20",  "A->C Trajectory", 10,
  "CAO3",   "20",  "A->C Trajectory", 10,
  "CAO4",   "20",  "A->C Trajectory", 10,
  "CAO5",   "20",  "A->C Trajectory", 10,

  # A->C trajectory at t24
  "CACO1",  "24",  "A->C Trajectory", 65,
  "CACO2",  "24",  "A->C Trajectory", 65,
  "CACO3",  "24",  "A->C Trajectory", 65,
  "CACO4",  "24",  "A->C Trajectory", 65,
  "CACO5",  "24",  "A->C Trajectory", 65,
  "CAO1",   "24",  "A->C Trajectory", 65,
  "CAO2",   "24",  "A->C Trajectory", 65,
  "CAO3",   "24",  "A->C Trajectory", 65,
  "CAO4",   "24",  "A->C Trajectory", 65,
  "CAO5",   "24",  "A->C Trajectory", 65,

  # C->A trajectory baseline at t18 (generation 0)
  "NACO1",  "18",  "C->A Trajectory", 0,
  "NACO2",  "18",  "C->A Trajectory", 0,
  "NACO3",  "18",  "C->A Trajectory", 0,
  "NACO4",  "18",  "C->A Trajectory", 0,
  "NACO5",  "18",  "C->A Trajectory", 0,
  "ANCO1",  "18",  "C->A Trajectory", 0,
  "ANCO2",  "18",  "C->A Trajectory", 0,
  "ANCO3",  "18",  "C->A Trajectory", 0,
  "ANCO4",  "18",  "C->A Trajectory", 0,
  "ANCO5",  "18",  "C->A Trajectory", 0,

  # C->A trajectory at t19
  "NACO1",  "19",  "C->A Trajectory", 16,
  "NACO2",  "19",  "C->A Trajectory", 16,
  "NACO3",  "19",  "C->A Trajectory", 16,
  "NACO4",  "19",  "C->A Trajectory", 16,
  "NACO5",  "19",  "C->A Trajectory", 16,
  "ANCO1",  "19",  "C->A Trajectory", 16,
  "ANCO2",  "19",  "C->A Trajectory", 16,
  "ANCO3",  "19",  "C->A Trajectory", 16,
  "ANCO4",  "19",  "C->A Trajectory", 16,
  "ANCO5",  "19",  "C->A Trajectory", 16,

  # C->A trajectory at t20
  "NACO1",  "20",  "C->A Trajectory", 29,
  "NACO2",  "20",  "C->A Trajectory", 29,
  "NACO3",  "20",  "C->A Trajectory", 29,
  "NACO4",  "20",  "C->A Trajectory", 29,
  "NACO5",  "20",  "C->A Trajectory", 29,
  "ANCO1",  "20",  "C->A Trajectory", 29,
  "ANCO2",  "20",  "C->A Trajectory", 29,
  "ANCO3",  "20",  "C->A Trajectory", 29,
  "ANCO4",  "20",  "C->A Trajectory", 29,
  "ANCO5",  "20",  "C->A Trajectory", 29,

  # C->A trajectory at t24
  "NACO1",  "24",  "C->A Trajectory", 182,
  "NACO2",  "24",  "C->A Trajectory", 182,
  "NACO3",  "24",  "C->A Trajectory", 182,
  "NACO4",  "24",  "C->A Trajectory", 182,
  "NACO5",  "24",  "C->A Trajectory", 182,
  "ANCO1",  "24",  "C->A Trajectory", 182,
  "ANCO2",  "24",  "C->A Trajectory", 182,
  "ANCO3",  "24",  "C->A Trajectory", 182,
  "ANCO4",  "24",  "C->A Trajectory", 182,
  "ANCO5",  "24",  "C->A Trajectory", 182
)

cat("Total sample x timepoint combinations to process:", nrow(sample_map), "\n\n")

# ---- COMPUTE HETEROZYGOSITY PER SITE PER SAMPLE ----
cat("Computing heterozygosity per site per sample...\n")

het_per_site <- bind_rows(lapply(1:nrow(sample_map), function(i) {

  prefix     <- sample_map$prefix[i]
  year       <- sample_map$year[i]
  group      <- sample_map$group[i]
  generation <- sample_map$generation[i]

  mj_col <- paste0(prefix, "_", year, "_mj")
  mi_col <- paste0(prefix, "_", year, "_mi")

  if (!all(c(mj_col, mi_col) %in% names(raw))) {
    cat("  WARNING: Missing columns for", prefix, year, "\n")
    return(NULL)
  }

  mj    <- raw[[mj_col]]
  mi    <- raw[[mi_col]]
  total <- mj + mi

  # Apply minimum coverage filter of 20x
  keep <- !is.na(total) & total >= 20

  freq <- mi[keep] / total[keep]
  het  <- 2 * freq * (1 - freq)

  data.frame(
    sample     = paste0(prefix, "_", year),
    population = prefix,
    year       = year,
    group      = group,
    generation = generation,
    chr        = raw$chr[keep],
    pos        = raw$pos[keep],
    coverage   = total[keep],
    freq_minor = freq,
    H          = het,
    stringsAsFactors = FALSE
  )
}))

cat("Done. Total site x sample observations:", nrow(het_per_site), "\n\n")

# ---- SUMMARY 1: PER SAMPLE ----
cat("Generating per-sample summary...\n")

het_per_sample <- het_per_site %>%
  group_by(group, population, year, generation, sample) %>%
  summarise(
    n_sites       = n(),
    mean_H        = round(mean(H, na.rm = TRUE), 6),
    median_H      = round(median(H, na.rm = TRUE), 6),
    sd_H          = round(sd(H, na.rm = TRUE), 6),
    mean_coverage = round(mean(coverage, na.rm = TRUE), 1),
    .groups = "drop"
  )

cat("Per-sample summary:\n")
print(as.data.frame(het_per_sample), row.names = FALSE)

# ---- SUMMARY 2: PER GROUP x TIMEPOINT ----
cat("\nGenerating per-group x timepoint summary...\n")

het_per_group <- as.data.frame(ungroup(het_per_sample)) %>%
  rename(rep_mean_H = mean_H) %>%
  group_by(group, generation) %>%
  summarise(
    n_populations = n(),
    mean_H        = mean(rep_mean_H, na.rm = TRUE),
    median_H      = median(rep_mean_H, na.rm = TRUE),
    sd_H          = sd(rep_mean_H, na.rm = TRUE),
    min_H         = min(rep_mean_H, na.rm = TRUE),
    max_H         = max(rep_mean_H, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(mean_H, median_H, sd_H, min_H, max_H),
                ~ round(.x, 6))) %>%
  arrange(group, generation)

cat("Per-group x timepoint summary:\n")
print(as.data.frame(het_per_group), row.names = FALSE)

# ---- SUMMARY 3: STATISTICAL COMPARISONS ----
# Holm-Bonferroni correction applied separately within each group:
# (1) each trajectory/founder group constitutes its own family of comparisons
#     against a common baseline, so corrections are not pooled across groups;
# (2) Holm controls the familywise error rate as stringently as Bonferroni
#     while gaining power by adapting the rejection threshold to the rank
#     of each p-value (Holm 1979).
cat("\nGenerating statistical comparisons...\n")

run_comparison <- function(df, grp, gen_before, gen_after) {
  before <- df %>%
    filter(group == grp, generation == gen_before) %>%
    pull(mean_H)
  after <- df %>%
    filter(group == grp, generation == gen_after) %>%
    pull(mean_H)

  if (length(before) < 2 | length(after) < 2) {
    return(data.frame(
      group             = grp,
      generation_before = gen_before,
      generation_after  = gen_after,
      mean_H_before     = NA,
      mean_H_after      = NA,
      delta_H           = NA,
      p_value           = NA,
      stringsAsFactors  = FALSE
    ))
  }

  t_res <- t.test(after, before)

  data.frame(
    group             = grp,
    generation_before = gen_before,
    generation_after  = gen_after,
    mean_H_before     = round(mean(before), 6),
    mean_H_after      = round(mean(after), 6),
    delta_H           = round(mean(after) - mean(before), 6),
    p_value           = t_res$p.value,
    stringsAsFactors  = FALSE
  )
}

comparisons <- bind_rows(
  # A-type founders: t18 vs t24
  run_comparison(het_per_sample, "A-type Founders", 0, 219),
  # C-type founders: t18 vs t24
  run_comparison(het_per_sample, "C-type Founders", 0, 219),
  # A->C trajectory: gen 0 vs each later timepoint
  run_comparison(het_per_sample, "A->C Trajectory", 0, 6),
  run_comparison(het_per_sample, "A->C Trajectory", 0, 10),
  run_comparison(het_per_sample, "A->C Trajectory", 0, 65),
  # C->A trajectory: gen 0 vs each later timepoint
  run_comparison(het_per_sample, "C->A Trajectory", 0, 16),
  run_comparison(het_per_sample, "C->A Trajectory", 0, 29),
  run_comparison(het_per_sample, "C->A Trajectory", 0, 182)
)

# Apply Holm correction within each group
comparisons <- comparisons %>%
  group_by(group) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "holm"),
    significance = case_when(
      is.na(p_adj)   ~ "NA",
      p_adj < 0.001  ~ "***",
      p_adj < 0.01   ~ "**",
      p_adj < 0.05   ~ "*",
      TRUE           ~ "ns"
    )
  ) %>%
  ungroup()

cat("Statistical comparisons (Holm corrected within group):\n")
print(as.data.frame(comparisons), row.names = FALSE)

# ---- SUMMARY 4: OVERALL GRAND MEANS PER GROUP ----
cat("\nGenerating overall grand means per group...\n")

het_grand <- as.data.frame(ungroup(het_per_sample)) %>%
  rename(rep_mean_H = mean_H) %>%
  group_by(group) %>%
  summarise(
    n_total_observations = n(),
    overall_mean_H       = round(mean(rep_mean_H, na.rm = TRUE), 6),
    overall_sd_H         = round(sd(rep_mean_H, na.rm = TRUE), 6),
    overall_min_H        = round(min(rep_mean_H, na.rm = TRUE), 6),
    overall_max_H        = round(max(rep_mean_H, na.rm = TRUE), 6),
    .groups = "drop"
  )

cat("Grand means per group:\n")
print(as.data.frame(het_grand), row.names = FALSE)

# ---- SAVE ALL OUTPUTS ----
cat("\nSaving outputs...\n")

write.csv(het_per_sample,
          file.path(output_dir, "heterozygosity_per_sample.csv"),
          row.names = FALSE)

write.csv(het_per_group,
          file.path(output_dir, "heterozygosity_per_group_timepoint.csv"),
          row.names = FALSE)

write.csv(comparisons,
          file.path(output_dir, "heterozygosity_statistical_comparisons.csv"),
          row.names = FALSE)

write.csv(het_grand,
          file.path(output_dir, "heterozygosity_grand_means_per_group.csv"),
          row.names = FALSE)

cat("\nAll files saved to:", output_dir, "\n")
cat("\nOutput files:\n")
cat("  heterozygosity_per_sample.csv - mean H per individual population per timepoint\n")
cat("  heterozygosity_per_group_timepoint.csv - mean H per group per timepoint with SD\n")
cat("  heterozygosity_statistical_comparisons.csv - t-tests with Holm correction within group\n")
cat("  heterozygosity_grand_means_per_group.csv - overall summary per group across all timepoints\n")
