# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(UpSetR)
library(dplyr)

data <- read.csv("Raw_GLM_Betabi_Results.csv", header = TRUE)

#“To avoid numerical underflow from sites with near-fixed allele frequencies, p-values were capped at a minimum of 1 × 10⁻³⁰ prior to false discovery rate correction.”
p_floor <- 1e-30

data <- data %>%
  mutate(
    p.value_capped = pmax(p.value, p_floor)
  )

# --- interaction term ---
interaction <- data %>%
  filter(term == "treatmentC>A:generation") %>%
  mutate(
    fdr  = p.adjust(p.value_capped, method = "BH"),
    p_bonf = p.adjust(p.value_capped, method = "bonferroni")
  )

# --- treatment term ---
treatment <- data %>%
  filter(term == "treatmentC>A") %>%
  mutate(
    fdr  = p.adjust(p.value_capped, method = "BH"),
    p_bonf = p.adjust(p.value_capped, method = "bonferroni")
  )

# --- generation term ---
generation <- data %>%
  filter(term == "generation") %>%
  mutate(
    fdr  = p.adjust(p.value_capped, method = "BH"),
    p_bonf = p.adjust(p.value_capped, method = "bonferroni")
  )


# --- significant subsets for Bonferroni (α = 0.05) ---
sig_interaction <- interaction %>%
  filter(p_bonf < 0.05)

sig_treatment <- treatment %>%
  filter(p_bonf < 0.05)

sig_generation <- generation %>%
  filter(p_bonf < 0.05)


#Make Upset Plot

# Assuming your data frames are:
# sig_generation, sig_interaction, sig_treatment
# Each has columns: chr and pos

# Create unique SNP IDs
sig_generation <- sig_generation %>%
  mutate(SNP = paste(chr, pos, sep = ":"))

sig_interaction <- sig_interaction %>%
  mutate(SNP = paste(chr, pos, sep = ":"))

sig_treatment <- sig_treatment %>%
  mutate(SNP = paste(chr, pos, sep = ":"))

# Create lists
snps <- list(
  Generation = sig_generation$SNP,
  Interaction = sig_interaction$SNP,
  Treatment = sig_treatment$SNP
)

# Create UpSet plot
upset(fromList(snps), order.by = "freq")

# Save UpSet plot as SVG
svg("UpSet_Plot.svg", width = 8, height = 6)  # specify size in inches
upset(fromList(snps), order.by = "freq")
dev.off()

# Filter out chr 4
interaction_plot <- interaction %>%
  filter(chr %in% c("X", "2L", "2R", "3L", "3R"))

# Bonferroni threshold for α = 0.05
bonf_threshold <- 0.05 / nrow(interaction_plot)

# Define custom colors
chr_colors <- c(
  "X"  = "#E69F00",  # orange
  "2L" = "#56B4E9",  # sky blue
  "2R" = "#009E73",  # bluish green
  "3L" = "#B22222",  # firebrick
  "3R" = "#0072B2"   # blue
)

# Convert positions to megabases
interaction_plot <- interaction_plot %>%
  mutate(pos_mb = pos / 1e6)

# Create the Manhattan plot
manhattan_plot <- ggplot(interaction_plot,
                         aes(x = pos_mb, y = -log10(p_bonf), color = chr)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  facet_wrap(~ chr, scales = "free_x", nrow = 1) +
  scale_color_manual(values = chr_colors) +
  labs(
    y = "-log10(Bonferroni adjusted p-value)",
    x = "Position (Mb)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none"
  ) +
  ylim(0,23)

# Save as PNG
ggsave(
  filename = "Manhattan_Plot.png",
  plot = manhattan_plot,
  width = 12,
  height = 4,
  dpi = 900
)
