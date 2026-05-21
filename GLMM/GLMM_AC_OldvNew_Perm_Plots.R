# Permutation-based FDR significance and plotting for Old vs New comparisons
# AvA: founder A (ACO, AO) vs reversed A (NACO, ANCO)
# CvC: founder C (CO, NCO) vs reversed C (CACO, CAO)
#
# SNPs are called significant if they pass both:
#   1. permutation-based FDR <= fdr_threshold
#   2. |mean allele frequency old - mean allele frequency new| >= effect_size_threshold

# libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# ── PATHS ─────────────────────────────────────────────────────────────────────
# Set data_dir to the path of your local data directory
data_dir <- "."
glmm_dir <- file.path(data_dir, "GLMM")
snp_path <- file.path(data_dir, "SNP_Table/AC_trajectory_modified_snp_table.csv")
# ──────────────────────────────────────────────────────────────────────────────

# ── SETTINGS ──────────────────────────────────────────────────────────────────
fdr_threshold         <- 0.01
effect_size_threshold <- 0.3
# ──────────────────────────────────────────────────────────────────────────────

# STEP 1: FDR function — identical logic to the main trajectory analysis
calculate_fdr <- function(obs_pvals, perm_df) {

  perm_sorted <- perm_df %>%
    group_by(perm) %>%
    arrange(pval, .by_group = TRUE) %>%
    summarise(pval = list(pval), .groups = "drop") %>%
    pull(pval)

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

# STEP 2: Read only header to identify columns (for effect size calculation later)
all_col_names <- colnames(read.csv(snp_path, header = TRUE, nrows = 0))

# STEP 3: Process one comparison — FDR + effect size
process_comparison <- function(obs_file, perm_file, col_pattern, comparison_name) {

  message("=== Processing ", comparison_name, " ===")

  obs <- read.table(file.path(glmm_dir, obs_file), header = TRUE) %>%
    filter(pval < 1, chr != "4", !is.na(pval))

  perm <- read.table(file.path(glmm_dir, perm_file), header = TRUE) %>%
    filter(pval < 1, chr != "4", !is.na(pval))

  obs <- obs %>%
    mutate(pval = ifelse(pval == 0, .Machine$double.xmin, pval))
  perm <- perm %>%
    mutate(pval = ifelse(pval == 0, .Machine$double.xmin, pval))

  obs_snps <- obs %>%
    mutate(SNP = paste(chr, pos, sep = ":")) %>%
    pull(SNP)

  perm_snps <- perm %>%
    mutate(SNP = paste(chr, pos, sep = ":")) %>%
    group_by(SNP) %>%
    summarise(n_perms = n_distinct(perm), .groups = "drop") %>%
    filter(n_perms == max(perm$perm)) %>%
    pull(SNP)

  common <- intersect(obs_snps, perm_snps)
  message("  SNPs in common set: ", length(common))

  obs <- obs %>%
    mutate(SNP = paste(chr, pos, sep = ":")) %>%
    filter(SNP %in% common)

  perm <- perm %>%
    mutate(SNP = paste(chr, pos, sep = ":")) %>%
    filter(SNP %in% common)

  obs$fdr <- calculate_fdr(obs$pval, perm)
  message("  FDR computed.")

  # ── effect size from raw counts, reading only needed columns ──
  data_cols   <- grep(col_pattern, all_col_names, value = TRUE)
  needed_cols <- c("chr", "pos", data_cols)

  col_classes <- rep("NULL", length(all_col_names))
  col_classes[all_col_names %in% needed_cols] <- NA

  snp_data <- read.csv(snp_path, header = TRUE, colClasses = col_classes)

  coverage <- snp_data[, grep("_total_coverage$", names(snp_data), value = TRUE)]
  minor    <- snp_data[, grep("_minor_count$",    names(snp_data), value = TRUE)]
  freq     <- minor / coverage

  # first 10 cols are "old" (founder), last 10 are "new" (reversed)
  old_mean <- rowMeans(freq[, 1:10],  na.rm = TRUE)
  new_mean <- rowMeans(freq[, 11:20], na.rm = TRUE)

  effect_df <- data.frame(
    chr         = snp_data$chr,
    pos         = snp_data$pos,
    effect_size = abs(old_mean - new_mean)
  ) %>%
    mutate(SNP = paste(chr, pos, sep = ":")) %>%
    select(SNP, effect_size)

  rm(snp_data, coverage, minor, freq)
  gc()

  obs <- obs %>%
    left_join(effect_df, by = "SNP")

  sig <- obs %>%
    filter(fdr <= fdr_threshold, effect_size >= effect_size_threshold)

  message("  ", comparison_name,
          " | FDR <= ", fdr_threshold, ": ", sum(obs$fdr <= fdr_threshold, na.rm = TRUE),
          " | FDR + effect size >= ", effect_size_threshold, ": ", nrow(sig))

  list(all = obs, sig = sig)
}

# STEP 4: Process both comparisons
AvA <- process_comparison(
  obs_file        = "GLMM_AvA.txt",
  perm_file       = "GLMM_AvA_Permuted.txt",
  col_pattern     = "^(ACO|AO|NACO|ANCO)[0-9]*_24(_minor_count|_total_coverage)$",
  comparison_name = "AvA"
)

CvC <- process_comparison(
  obs_file        = "GLMM_CvC.txt",
  perm_file       = "GLMM_CvC_Permuted.txt",
  col_pattern     = "^(CO|NCO|CACO|CAO)[0-9]*_24(_minor_count|_total_coverage)$",
  comparison_name = "CvC"
)

# STEP 5: Save significant SNP CSVs
AvA$sig %>%
  select(chr, pos, pval, fdr, effect_size) %>%
  arrange(chr, pos) %>%
  write.csv(file.path(glmm_dir, "Significant_SNPs_AvA.csv"), row.names = FALSE)

CvC$sig %>%
  select(chr, pos, pval, fdr, effect_size) %>%
  arrange(chr, pos) %>%
  write.csv(file.path(glmm_dir, "Significant_SNPs_CvC.csv"), row.names = FALSE)

message("Significant SNP files saved.")

# STEP 6: p-value threshold corresponding to FDR = fdr_threshold
# If observed SNPs meet FDR <= threshold, return the largest such p-value.
# Otherwise compute the theoretical p-value where FDR would equal the threshold
# given the null distribution — this places the threshold line above the data
# so the plot visually communicates that no SNPs reach significance.
get_perm_threshold <- function(obs_df, perm_file) {

  # first try: any observed SNP meeting the FDR cutoff
  below <- obs_df %>% filter(fdr <= fdr_threshold)
  if (nrow(below) > 0) {
    return(max(below$pval))
  }

  # fall back to theoretical threshold from null distribution
  # FDR(p) = (null p-values <= p) / (observed p-values <= p)
  # find largest p where this ratio crosses fdr_threshold
  message("  No observed SNPs at FDR <= ", fdr_threshold,
          ". Computing theoretical threshold from null distribution.")

  perm <- read.table(file.path(glmm_dir, perm_file), header = TRUE) %>%
    filter(pval < 1, chr != "4", !is.na(pval)) %>%
    mutate(pval = ifelse(pval == 0, .Machine$double.xmin, pval))

  perm_sorted <- perm %>%
    group_by(perm) %>%
    arrange(pval, .by_group = TRUE) %>%
    summarise(pval = list(pval), .groups = "drop") %>%
    pull(pval)

  perm_matrix <- do.call(cbind, perm_sorted)
  null_dist   <- sort(apply(perm_matrix, 1, min))
  obs_sorted  <- sort(obs_df$pval)

  # scan candidate p-values from very small to moderate; find largest where FDR <= target
  candidates <- 10 ^ seq(-20, -1, length.out = 2000)
  fdrs <- sapply(candidates, function(p) {
    obs_below <- sum(obs_sorted <= p)
    if (obs_below == 0) return(NA)
    sum(null_dist <= p) / obs_below
  })

  valid <- !is.na(fdrs) & fdrs <= fdr_threshold
  if (!any(valid)) return(NA_real_)
  max(candidates[valid])
}

AvA_threshold <- get_perm_threshold(AvA$all, "GLMM_AvA_Permuted.txt")
CvC_threshold <- get_perm_threshold(CvC$all, "GLMM_CvC_Permuted.txt")

message("AvA p-value threshold at FDR ", fdr_threshold, ": ", AvA_threshold)
message("CvC p-value threshold at FDR ", fdr_threshold, ": ", CvC_threshold)

# STEP 7: Manhattan plots

chr_colors <- c(
  "X"  = "#E69F00",
  "2L" = "#56B4E9",
  "2R" = "#009E73",
  "3L" = "#B22222",
  "3R" = "#0072B2"
)

color_scale <- c(chr_colors, nonsig = "gray80")

flag_significance <- function(df) {
  df %>%
    mutate(
      color_group = ifelse(
        fdr <= fdr_threshold & effect_size < effect_size_threshold,
        "nonsig",
        as.character(chr)
      )
    )
}

AvA_data <- flag_significance(AvA$all)
CvC_data <- flag_significance(CvC$all)

y_max <- max(-log10(c(AvA$all$pval, CvC$all$pval)), na.rm = TRUE)

make_manhattan <- function(data, threshold, y_max) {
  p <- ggplot(data, aes(x = pos / 1e6, y = -log10(pval), color = color_group)) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_manual(values = color_scale) +
    facet_wrap(~ chr, scales = "free_x", nrow = 1) +
    labs(y = "-log10(p-value)", x = "Position (Mb)") +
    ylim(0, y_max) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(size = 12, face = "plain"),
      axis.title       = element_text(size = 12, face = "plain"),
      panel.spacing    = unit(0.5, "lines"),
      legend.position  = "none",
      plot.title       = element_text(hjust = 0.5)
    )

  # only draw threshold line if we have a valid threshold
  if (!is.na(threshold)) {
    p <- p + geom_hline(yintercept = -log10(threshold),
                        linetype = "dashed", color = "black")
  }
  p
}

AvA_plot <- make_manhattan(AvA_data, AvA_threshold, y_max)
CvC_plot <- make_manhattan(CvC_data, CvC_threshold, y_max)

combined <- ggarrange(
  AvA_plot,
  CvC_plot,
  ncol = 1, nrow = 2,
  labels = c("A", "B"),
  heights = c(1, 1)
)

ggsave(
  file.path(glmm_dir, "GLMM_AvA_CvC_Perm_faceted.tiff"),
  plot = combined,
  width = 14, height = 10, units = "in",
  dpi = 900, compression = "lzw"
)

ggsave(
  file.path(glmm_dir, "GLMM_AvA_CvC_Perm_faceted.png"),
  plot = combined,
  width = 14, height = 10, units = "in",
  dpi = 900
)

message("Combined plot saved.")
