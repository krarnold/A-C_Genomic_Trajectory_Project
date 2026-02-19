# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Read in data
data <- read.csv("SNP_Table/AC_trajectory_modified_snp_table.csv", header = TRUE)

# Remove columns corresponding to 24th timepoint if present (you later remap _24 anyway)
# NOTE: If you truly want to keep endpoints, comment this block out.
cols_to_remove <- grep("^(ACO|AO|CO|NCO)[0-9]_24", names(data), value = TRUE)
data <- data[, !names(data) %in% cols_to_remove]

# Rename map for samples
rename_map <- c(
  CACO1="A>C1", CACO2="A>C2", CACO3="A>C3", CACO4="A>C4", CACO5="A>C5",
  CAO1="A>C6",  CAO2="A>C7",  CAO3="A>C8",  CAO4="A>C9",  CAO5="A>C10",
  NACO1="C>A1", NACO2="C>A2", NACO3="C>A3", NACO4="C>A4", NACO5="C>A5",
  ANCO1="C>A6", ANCO2="C>A7", ANCO3="C>A8", ANCO4="C>A9", ANCO5="C>A10"
)

names(data) <- names(data) |>
  sapply(\(x) {
    for (old in names(rename_map)) {
      if (startsWith(x, old)) {
        x <- sub(old, rename_map[[old]], x)
        break
      }
    }
    x
  }, USE.NAMES = FALSE)

# Replace total coverage with major count
total_cols <- grep("_total_coverage$", names(data), value = TRUE)

for (col in total_cols) {
  minor_col <- sub("_total_coverage$", "_minor_count", col)
  new_colname <- sub("_total_coverage$", "_major_count", col)
  
  data[[col]] <- data[[col]] - data[[minor_col]]
  names(data)[names(data) == col] <- new_colname
}

# Change year designations to generations
# A>C samples
names(data) <- sub("A>C([0-9]+)_18", "A>C\\1_0",  names(data))
names(data) <- sub("A>C([0-9]+)_19", "A>C\\1_6",  names(data))
names(data) <- sub("A>C([0-9]+)_20", "A>C\\1_10", names(data))
names(data) <- sub("A>C([0-9]+)_24", "A>C\\1_65", names(data))

# C>A samples
names(data) <- sub("C>A([0-9]+)_18", "C>A\\1_0",   names(data))
names(data) <- sub("C>A([0-9]+)_19", "C>A\\1_16",  names(data))
names(data) <- sub("C>A([0-9]+)_20", "C>A\\1_29",  names(data))
names(data) <- sub("C>A([0-9]+)_24", "C>A\\1_182", names(data))

# Make table long format
data_long <- data %>%
  pivot_longer(
    cols = -c(chr, pos, Major_Nuc, Minor_Nuc),
    names_to = c("popdir", "replicate", "generation", ".value"),
    names_pattern = "^(A>C|C>A)(\\d+)_([0-9]+)_(minor|major)_count$"
  ) %>%
  mutate(
    population = paste0(popdir, replicate)  # e.g., "A>C1"
  ) %>%
  transmute(
    population,
    treatment = popdir,
    generation = as.numeric(generation),
    chr, pos, Major_Nuc, Minor_Nuc,
    Minor_Count = minor,
    Major_Count = major
  )

rm(data)

# Add heterozygosity columns
data_long <- data_long %>%
  mutate(
    total_count = Minor_Count + Major_Count,
    freq_minor  = Minor_Count / total_count,
    Heterozygosity = 2 * freq_minor * (1 - freq_minor)
  )

# Summarize mean heterozygosity per population x generation x treatment
het_summary <- data_long %>%
  group_by(treatment, population, generation) %>%
  summarise(
    mean_heterozygosity = mean(Heterozygosity, na.rm = TRUE),
    n_sites = n(),
    .groups = "drop"
  )

#Compare each intermediate to Generation 0

# Fixed function: treatment_name argument is not shadowed
get_pval_label <- function(df, treatment_name, before_gen, after_gen) {
  het_before <- df %>%
    filter(treatment == treatment_name, generation == before_gen) %>%
    pull(mean_heterozygosity)
  
  het_after <- df %>%
    filter(treatment == treatment_name, generation == after_gen) %>%
    pull(mean_heterozygosity)
  
  # Guardrails: if missing timepoints, return NA row
  if (length(het_before) < 2 || length(het_after) < 2) {
    return(data.frame(
      treatment = treatment_name,
      before = before_gen,
      after  = after_gen,
      p_value = NA_real_,
      label = "NA",
      stringsAsFactors = FALSE
    ))
  }
  
  t_res <- t.test(het_after, het_before)
  
  p_val <- t_res$p.value
  
  label <- if (p_val < 0.001) {
    "***"
  } else if (p_val < 0.01) {
    "**"
  } else if (p_val < 0.05) {
    "*"
  } else {
    "ns"
  }
  
  data.frame(
    treatment = treatment_name,
    before = before_gen,
    after  = after_gen,
    p_value = p_val,
    label = label,
    stringsAsFactors = FALSE
  )
}

# Build comparisons dynamically: for each treatment, compare all generations (except 0) vs 0
comparisons <- het_summary %>%
  distinct(treatment, generation) %>%
  group_by(treatment) %>%
  summarise(
    after = list(sort(unique(generation[generation != 0]))),
    .groups = "drop"
  ) %>%
  tidyr::unnest(after) %>%
  transmute(treatment, before = 0, after)

labels_df <- comparisons %>%
  rowwise() %>%
  do(get_pval_label(het_summary, .$treatment, .$before, .$after)) %>%
  ungroup()

# Multiple-testing correction within each treatment (Holm)
labels_df <- labels_df %>%
  group_by(treatment) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "holm"),
    label = case_when(
      is.na(p_adj) ~ "NA",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  ) %>%
  ungroup()


# Plot helpers

# Ensure consistent ordering within each facet
gen_levels <- sort(unique(het_summary$generation))
het_summary <- het_summary %>%
  mutate(generation_f = factor(generation, levels = gen_levels))

# Place labels above the top box in each facet
y_pos <- het_summary %>%
  group_by(treatment) %>%
  summarise(y = max(mean_heterozygosity, na.rm = TRUE) + 0.01, .groups = "drop")

labels_df <- labels_df %>%
  left_join(y_pos, by = "treatment") %>%
  mutate(after_f = factor(after, levels = gen_levels))


# Plot 1: quick view (all timepoints)
het_boxplot_quick <- ggplot(
  het_summary,
  aes(x = generation_f, y = mean_heterozygosity, fill = treatment)
) +
  geom_boxplot(outlier.shape = 21, color = "black") +
  facet_wrap(~ treatment, scales = "free_x") +
  scale_fill_manual(values = c("A>C" = "#377eb8", "C>A" = "#e41a1c")) +
  scale_y_continuous(limits = c(0.15, 0.3)) +
  labs(
    x = "Generation",
    y = "Mean Genome-wide Heterozygosity"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "none"
  )

print(het_boxplot_quick)


# Plot 2: publication SVG (all timepoints + stars)
het_boxplot_svg <- ggplot(
  het_summary,
  aes(x = generation_f, y = mean_heterozygosity, fill = treatment)
) +
  geom_boxplot(outlier.shape = 21, color = "black", size = 1) +
  facet_wrap(~ treatment, scales = "free_x") +
  scale_fill_manual(values = c("A>C" = "#377eb8", "C>A" = "#e41a1c")) +
  scale_y_continuous(limits = c(0.15, 0.3)) +
  labs(
    x = "Generation",
    y = "Mean Genome-wide Heterozygosity"
  ) +
  theme_bw(base_size = 20) +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1.2),
    strip.text = element_text(color = "black", face = "bold", size = 20),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    legend.position = "none"
  ) +
  geom_text(
    data = labels_df,
    aes(x = after_f, y = y, label = label),
    inherit.aes = FALSE,
    size = 8
  )

print(het_boxplot_svg)

ggsave(
  "Mean_Heterozygosity_Boxplot.svg",
  plot = het_boxplot_svg,
  device = "svg",
  width = 14,
  height = 8,
  units = "in"
)
