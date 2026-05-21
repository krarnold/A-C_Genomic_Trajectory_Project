# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(UpSetR)
library(patchwork)

# ── PATHS ─────────────────────────────────────────────────────────────────────
# Set data_dir to the path of your local data directory
data_dir <- "."
# ──────────────────────────────────────────────────────────────────────────────

# ── SETTINGS ──────────────────────────────────────────────────────────────────
# SNPs must pass both FDR <= fdr_threshold AND effect size >= effect_size_threshold
fdr_threshold         <- 0.01
effect_size_threshold <- 0.3
# ──────────────────────────────────────────────────────────────────────────────

# STEP 1: Load observed and permuted results
observed <- read.csv(file.path(data_dir, "Raw_GLM_Betabi_Results.csv"),      header = TRUE)
permuted <- read.csv(file.path(data_dir, "Permuted_GLM_Betabi_Results.csv"), header = TRUE)

# replace p-values of exactly 0 with .Machine$double.xmin in both datasets
observed <- observed %>%
  mutate(p.value = ifelse(p.value == 0, .Machine$double.xmin, p.value))

permuted <- permuted %>%
  mutate(p.value = ifelse(p.value == 0, .Machine$double.xmin, p.value))

terms <- c("treatmentC>A", "generation", "treatmentC>A:generation")

# STEP 2: FDR function
calculate_fdr <- function(obs_pvals, perm_df) {

  perm_sorted <- perm_df %>%
    group_by(perm) %>%
    arrange(p.value, .by_group = TRUE) %>%
    summarise(p.value = list(p.value), .groups = "drop") %>%
    pull(p.value)

  perm_matrix <- do.call(cbind, perm_sorted)
  null_dist   <- sort(apply(perm_matrix, 1, min))

  obs_order  <- order(obs_pvals)
  sorted_obs <- obs_pvals[obs_order]
  n_obs      <- length(sorted_obs)

  i_vec   <- findInterval(sorted_obs, null_dist)
  fdr_raw <- i_vec / seq_len(n_obs)

  qvals <- rev(cummin(rev(fdr_raw)))
  qvals <- pmin(qvals, 1)

  fdr_unsorted            <- numeric(n_obs)
  fdr_unsorted[obs_order] <- qvals
  fdr_unsorted
}

# STEP 3: Apply FDR calculation to each term
fdr_list <- list()
sig_list_fdr_only <- list()

for (term_name in terms) {

  message("Calculating FDR for term: ", term_name)

  obs_term  <- observed %>% filter(term == term_name, !is.na(p.value))
  perm_term <- permuted %>% filter(term == term_name, !is.na(p.value))

  obs_snps <- obs_term %>%
    mutate(SNP = paste(chr, pos, sep = ":")) %>%
    pull(SNP)

  perm_snps <- perm_term %>%
    mutate(SNP = paste(chr, pos, sep = ":")) %>%
    group_by(SNP) %>%
    summarise(n_perms = n_distinct(perm), .groups = "drop") %>%
    filter(n_perms == max(perm_term$perm)) %>%
    pull(SNP)

  common_snps <- intersect(obs_snps, perm_snps)
  message("  SNPs in common set: ", length(common_snps))

  obs_term <- obs_term %>%
    mutate(SNP = paste(chr, pos, sep = ":")) %>%
    filter(SNP %in% common_snps) %>%
    select(-SNP)

  perm_term <- perm_term %>%
    mutate(SNP = paste(chr, pos, sep = ":")) %>%
    filter(SNP %in% common_snps) %>%
    select(-SNP)

  obs_term$fdr <- calculate_fdr(obs_term$p.value, perm_term)

  fdr_list[[term_name]] <- obs_term

  sig_list_fdr_only[[term_name]] <- obs_term %>%
    filter(fdr <= fdr_threshold) %>%
    mutate(SNP = paste(chr, pos, sep = ":"))

  message("  Significant SNPs (FDR <= ", fdr_threshold, "): ",
          nrow(sig_list_fdr_only[[term_name]]))
}

# ── CALCULATE TERM-SPECIFIC EFFECT SIZES ──────────────────────────────────────

message("Loading long format data and calculating effect sizes...")

data_long <- read.csv(file.path(data_dir, "Long_format_SNP_Data.csv"),
                      header = TRUE) %>%
  mutate(
    treatment  = factor(treatment),
    population = factor(population),
    generation = as.numeric(generation)
  )

# filter to union of FDR-significant SNPs for effect size calculation
all_fdr_sig_snps <- unique(c(
  sig_list_fdr_only[["treatmentC>A"]]$SNP,
  sig_list_fdr_only[["generation"]]$SNP,
  sig_list_fdr_only[["treatmentC>A:generation"]]$SNP
))

# keep full long format for the antiparallel plot later — filter to sig for effect calcs
data_sig <- data_long %>%
  mutate(SNP = paste(chr, pos, sep = ":")) %>%
  filter(SNP %in% all_fdr_sig_snps)

data_freq <- data_sig %>%
  mutate(freq = Minor_Count / (Minor_Count + Major_Count)) %>%
  group_by(SNP, treatment, generation) %>%
  summarise(mean_freq = mean(freq, na.rm = TRUE), .groups = "drop")

# treatment effect size: |mean freq A>C - mean freq C>A| averaged across timepoints
effect_treatment <- data_freq %>%
  group_by(SNP, treatment) %>%
  summarise(mean_freq_overall = mean(mean_freq, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = treatment, values_from = mean_freq_overall) %>%
  mutate(effect_size = abs(`A>C` - `C>A`)) %>%
  select(SNP, effect_size)

# per-trajectory start/end frequencies — reused for generation, interaction, antiparallel plot
freq_AC_gen <- data_freq %>%
  filter(treatment == "A>C", generation %in% c(0, 65)) %>%
  pivot_wider(names_from = generation, values_from = mean_freq,
              names_prefix = "gen_") %>%
  mutate(delta_AC = gen_65 - gen_0) %>%
  select(SNP, delta_AC)

freq_CA_gen <- data_freq %>%
  filter(treatment == "C>A", generation %in% c(0, 182)) %>%
  pivot_wider(names_from = generation, values_from = mean_freq,
              names_prefix = "gen_") %>%
  mutate(delta_CA = gen_182 - gen_0) %>%
  select(SNP, delta_CA)

# generation effect size: mean |freq change over time| averaged across trajectories
effect_generation <- freq_AC_gen %>%
  inner_join(freq_CA_gen, by = "SNP") %>%
  mutate(effect_size = (abs(delta_AC) + abs(delta_CA)) / 2) %>%
  select(SNP, effect_size)

# interaction effect size: |delta A>C - delta C>A| in terms of A-allele
a_allele_freq <- data_freq %>%
  filter(treatment == "A>C", generation == 0) %>%
  select(SNP, freq_AC_start = mean_freq)

# also store oriented deltas for the antiparallel plot
deltas_oriented <- freq_AC_gen %>%
  inner_join(freq_CA_gen, by = "SNP") %>%
  inner_join(a_allele_freq, by = "SNP") %>%
  mutate(
    flip     = freq_AC_start < 0.5,
    delta_AC = ifelse(flip, -delta_AC, delta_AC),
    delta_CA = ifelse(flip, -delta_CA, delta_CA)
  ) %>%
  select(SNP, delta_AC, delta_CA)

effect_interaction <- deltas_oriented %>%
  mutate(effect_size = abs(delta_AC - delta_CA)) %>%
  select(SNP, effect_size)

effect_sizes_list <- list(
  "treatmentC>A"            = effect_treatment,
  "generation"              = effect_generation,
  "treatmentC>A:generation" = effect_interaction
)

# free memory
rm(data_long, data_sig)
gc()

# ── APPLY COMBINED FILTER ─────────────────────────────────────────────────────

sig_list <- list()

for (term_name in terms) {

  sig_list[[term_name]] <- sig_list_fdr_only[[term_name]] %>%
    inner_join(effect_sizes_list[[term_name]], by = "SNP") %>%
    filter(effect_size >= effect_size_threshold)

  message("Term: ", term_name,
          " | FDR-only: ", nrow(sig_list_fdr_only[[term_name]]),
          " | FDR + effect size: ", nrow(sig_list[[term_name]]))
}

# ── SAVE SIGNIFICANT SNPs ─────────────────────────────────────────────────────

for (term_name in terms) {

  clean_name <- gsub("[^a-zA-Z0-9]", "_", term_name)
  out_name   <- paste0("Significant_SNPs_", clean_name, ".csv")

  sig_list[[term_name]] %>%
    select(chr, pos, estimate, std.error, statistic, p.value, fdr, effect_size) %>%
    arrange(chr, pos) %>%
    write.csv(file.path(data_dir, out_name), row.names = FALSE)

  message("Saved: ", out_name)
}

# ── UPSET PLOT ────────────────────────────────────────────────────────────────
# UpSetR requires at least two non-empty sets — skip if any set is empty

snps <- list(
  Generation  = sig_list[["generation"]]$SNP,
  Interaction = sig_list[["treatmentC>A:generation"]]$SNP,
  Treatment   = sig_list[["treatmentC>A"]]$SNP
)

set_sizes <- sapply(snps, length)
message("Set sizes for UpSet plot: ",
        paste(names(set_sizes), "=", set_sizes, collapse = ", "))

if (all(set_sizes > 0)) {
  tiff(file.path(data_dir, "UpSet_Plot.tiff"),
       width = 12, height = 10, units = "in", res = 900, compression = "lzw")
  upset(fromList(snps),
        order.by = "freq",
        text.scale = c(1.5, 1.5, 1.5, 1.2, 1.5, 1.3))
  dev.off()
  message("UpSet plot saved.")
} else {
  message("UpSet plot skipped — at least one set is empty under current thresholds.")
}

# ── MANHATTAN PLOT ────────────────────────────────────────────────────────────
# plots all SNPs with y-axis = -log10(p-value)
# threshold line is drawn at the raw p-value corresponding to FDR <= fdr_threshold
# p-values below 1e-30 are excluded from the plot (floor artefacts)

interaction_obs <- fdr_list[["treatmentC>A:generation"]] %>%
  filter(chr %in% c("X", "2L", "2R", "3L", "3R")) %>%
  filter(p.value >= 1e-30) %>%
  mutate(
    p.value_plot = p.value,
    pos_mb       = pos / 1e6
  )

perm_threshold <- fdr_list[["treatmentC>A:generation"]] %>%
  filter(chr %in% c("X", "2L", "2R", "3L", "3R"),
         fdr <= fdr_threshold) %>%
  summarise(threshold = max(p.value)) %>%
  pull(threshold)

chr_colors <- c(
  "X"  = "#E69F00",
  "2L" = "#56B4E9",
  "2R" = "#009E73",
  "3L" = "#B22222",
  "3R" = "#0072B2"
)

manhattan_plot <- ggplot(interaction_obs,
                         aes(x = pos_mb, y = -log10(p.value_plot), color = chr)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(perm_threshold),
             linetype = "dashed", color = "black") +
  facet_wrap(~ chr, scales = "free_x", nrow = 1) +
  scale_color_manual(values = chr_colors) +
  labs(
    y = "-log10(p-value)",
    x = "Position (Mb)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.spacing    = unit(0.5, "lines"),
    legend.position  = "none"
  )

ggsave(
  filename = file.path(data_dir, "Manhattan_Plot.tiff"),
  plot     = manhattan_plot,
  width    = 12,
  height   = 4,
  dpi      = 900
)

message("Manhattan plot saved.")

# also save individual Manhattan object for multipanel figure below

# ── ANTIPARALLEL PLOT ─────────────────────────────────────────────────────────
# hexbin density plot of allele frequency change in A>C vs C>A, A-allele oriented
# uses SNPs passing FDR <= fdr_threshold only (no effect size filter)
# this plot is a diagnostic showing the interaction test detects anti-parallel
# responses — applying the effect size filter here would bias the visualization
# toward the most anti-parallel SNPs by construction

sig_interaction_snps <- sig_list_fdr_only[["treatmentC>A:generation"]]$SNP

plot_data <- deltas_oriented %>%
  filter(SNP %in% sig_interaction_snps)

message("Antiparallel plot SNPs: ", nrow(plot_data))

antiparallel_plot <- ggplot(plot_data, aes(x = delta_CA, y = delta_AC)) +
  geom_hex(bins = 80) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_viridis_c(option = "magma", trans = "log10",
                       name = "SNP count\n(log10)") +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  coord_fixed() +
  labs(
    x = "Allele frequency change C→A",
    y = "Allele frequency change A→C"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(data_dir, "Antiparallel_Plot.tiff"),
  plot     = antiparallel_plot,
  width    = 6,
  height   = 6,
  dpi      = 900
)

message("Antiparallel plot saved.")

# ── MULTIPANEL FIGURE ─────────────────────────────────────────────────────────
# stacked layout: Manhattan (A) on top, Antiparallel (B) below
# antiparallel is centered beneath the wide Manhattan

combined_plot <- manhattan_plot / antiparallel_plot +
  plot_layout(heights = c(1, 3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave(
  filename = file.path(data_dir, "Combined_Figure.tiff"),
  plot     = combined_plot,
  width    = 12,
  height   = 16,
  dpi      = 900
)

message("Combined multipanel figure saved.")
