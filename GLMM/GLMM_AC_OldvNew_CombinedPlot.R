# Load required packages
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)

# --- Load and prepare AvA data ---
newA_A <- read.table("GLMM/GLMM_AvA.txt", header = TRUE)
newA_A <- subset(newA_A, pval < 1)
newA_A <- newA_A %>% filter(chr != "4")
newA_A$p_bonf_a <- p.adjust(newA_A$pval, method = "bonferroni")
sigA_A <- subset(newA_A, p_bonf_a < 0.05)

# --- Load and prepare CvC data ---
newC_C <- read.table("GLMM/GLMM_CVC.txt", header = TRUE)
newC_C <- subset(newC_C, pval < 1)
newC_C <- newC_C %>% filter(chr != "4")
newC_C$p_bonf_c <- p.adjust(newC_C$pval, method = "bonferroni")
sigC_C <- subset(newC_C, p_bonf_c < 0.05)

# --- Define colorblind-friendly palette for chromosomes ---
chr_colors <- c(
  "X"  = "#E69F00",  # orange
  "2L" = "#56B4E9",  # sky blue
  "2R" = "#009E73",  # bluish green
  "3L" = "#B22222",  # firebrick
  "3R" = "#0072B2"   # blue
)

# --- AvA plot ---
AvA_faceted <- ggplot(newA_A, aes(x = pos/1e6, y = -log10(pval), color = chr)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(values = chr_colors) +
  geom_hline(yintercept = -log10(0.05/nrow(newA_A)), linetype = "dashed", color = "black") +
  facet_wrap(~ chr, scales = "free_x", nrow = 1) +
  labs(y = "-log10(p-value)", x = "Position (Mb)") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "plain"), # chromosome labels
    axis.title = element_text(size = 12, face = "plain"),  # x and y axis labels
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  ylim(0, 12)

# --- CvC plot ---
CvC_faceted <- ggplot(newC_C, aes(x = pos/1e6, y = -log10(pval), color = chr)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(values = chr_colors) +
  geom_hline(yintercept = -log10(0.05/nrow(newC_C)), linetype = "dashed", color = "black") +
  facet_wrap(~ chr, scales = "free_x", nrow = 1) +
  labs(y = "-log10(p-value)", x = "Position (Mb)") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "plain"), # chromosome labels
    axis.title = element_text(size = 12, face = "plain"),  # x and y axis labels
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  ylim(0, 12)

# --- Combine plots vertically with AvA on top ---
combined_faceted <- ggarrange(
  AvA_faceted,
  CvC_faceted,
  ncol = 1, nrow = 2,
  labels = c("A", "B"),  
  heights = c(1, 1)
)

# --- Display combined plot ---
combined_faceted


# Save as SVG
ggsave(
  "GLMM_AvA_CvC_Bonferroni_faceted.svg",
  plot = combined_faceted,
  width = 14,
  height = 10
)

# Save PNG
ggsave(
  "GLMM_AvA_CvC_Bonferroni_faceted.png",
  plot = combined_faceted,
  width = 14,
  height = 10,
  units = "in",  
  dpi = 900      
)
