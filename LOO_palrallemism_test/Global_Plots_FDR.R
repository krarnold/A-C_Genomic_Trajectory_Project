setwd("~/Library/CloudStorage/Dropbox/AC Trajectory Genomic Data/LOO_palrallemism_test")

# ===============================================================
# Manhattan plots for A>C and C>A global results
# Significance: FDR < 0.05
# Labeled panels: A. A to C Trajectory / B. C to A Trajectory
# ===============================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# --- Load data ---
AC_global <- read.csv("AC_results_all.csv", header = TRUE)
CA_global <- read.csv("CA_results_all.csv", header = TRUE)

# --- Process data ---
process_for_manhattan <- function(df, label) {
  df %>%
    filter(chr %in% c("X", "2L", "2R", "3L", "3R")) %>%
    mutate(
      fdr = p.adjust(p.value, method = "BH"),
      sig = fdr < 0.05,
      neglog10p = -log10(fdr),
      Treatment = label
    )
}

AC_plot_data <- process_for_manhattan(AC_global, "A>C")
CA_plot_data <- process_for_manhattan(CA_global, "C>A")

# --- Plot function ---
make_manhattan <- function(df, title_label) {
  ggplot(df, aes(x = pos, y = neglog10fdr, color = chr)) +
    
    # ---- 1. Points first ----
  geom_point(size = 0.5, alpha = 0.7) +
    
    # ---- 2. Threshold lines drawn last (on top) ----
  geom_hline(yintercept = thr05, linetype = "dashed",
             color = "black", linewidth = 0.65, inherit.aes = FALSE) +
    geom_hline(yintercept = thr01, linetype = "dashed",
               color = "gray40", linewidth = 0.65, inherit.aes = FALSE) +
    
    # ---- 3. Labels for the lines (also on top) ----
  annotate("text", x = -Inf, y = thr05, label = "0.05",
           hjust = -0.1, vjust = -0.5, size = 3.2, color = "black",
           inherit.aes = FALSE) +
    annotate("text", x = -Inf, y = thr01, label = "0.01",
             hjust = -0.1, vjust = -0.5, size = 3.2, color = "gray40",
             inherit.aes = FALSE) +
    
    scale_color_manual(values = c(
      "2L" = "#4E79A7", "2R" = "#A0CBE8",
      "3L" = "#F28E2B", "3R" = "#FFBE7D",
      "X"  = "#59A14F"
    )) +
    facet_wrap(~ chr, scales = "free_x", nrow = 1) +
    labs(
      title = title_label,
      x = "Genomic position",
      y = expression(-log[10]("FDR-adjusted p-value"))
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.spacing = unit(0.4, "lines"),
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold", hjust = 0),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
}


# --- Build labeled plots ---
p_AC <- make_manhattan(AC_plot_data, "A. A to C Trajectory")
p_CA <- make_manhattan(CA_plot_data, "B. C to A Trajectory")

# --- Combine into two rows ---
combined_plot <- p_AC / p_CA 

# --- Display and save ---
#print(combined_plot)

ggsave("Global_Manhattan_AC_CA_FDR05_labeled.png",
       combined_plot, width = 15, height = 8, dpi = 300)
