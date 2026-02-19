# ===============================
# Master PCA + Mixed-Effect Model Script with Single Timepoint + CSV Export
# and Flexible Faceted PCA Plot
# ===============================

# --- Libraries ---
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(scales)

# ===============================
# 1) Read and clean SNP data
# ===============================

data <- read.csv("SNP_Table/AC_trajectory_modified_snp_table.csv", header = TRUE)
data <- data[, c("chr", "pos", grep("_24_", names(data), value = TRUE))]

# Rename map for samples
rename_map <- c(
  # Trajectory Lines -- "A>C" / "C>A" labels
  CACO1="A>C1", CACO2="A>C2", CACO3="A>C3", CACO4="A>C4", CACO5="A>C5",
  CAO1="A>C6",  CAO2="A>C7",  CAO3="A>C8",  CAO4="A>C9",  CAO5="A>C10",
  NACO1="C>A1", NACO2="C>A2", NACO3="C>A3", NACO4="C>A4", NACO5="C>A5",
  ANCO1="C>A6", ANCO2="C>A7", ANCO3="C>A8", ANCO4="C>A9", ANCO5="C>A10",
  
  # Ancestrals -- A1–A10 / C1–C10
  ACO1="A1", ACO2="A2", ACO3="A3", ACO4="A4", ACO5="A5",
  AO1="A6",  AO2="A7",  AO3="A8",  AO4="A9",  AO5="A10",
  CO1="C1",  CO2="C2",  CO3="C3",  CO4="C4",  CO5="C5",
  NCO1="C6", NCO2="C7", NCO3="C8", NCO4="C9", NCO5="C10"
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

# Calculate SNP frequencies
samples <- gsub("_minor_count$", "", grep("_minor_count$", names(data), value = TRUE))
freq <- data[, c("chr", "pos")]

for (s in samples) {
  k <- data[[paste0(s, "_minor_count")]]
  n <- data[[paste0(s, "_total_coverage")]]
  # Add 1 only to zeros
  k[k == 0] <- k[k == 0] + 1
  freq[[s]] <- k / n
}

# ===============================
# 2) PCA on SNP frequencies
# ===============================
freq_mat <- as.matrix(freq[, -c(1,2)])
freq_mat <- t(freq_mat)
pca_res <- prcomp(freq_mat, center = FALSE, scale. = TRUE)
var_exp <- round((pca_res$sdev)^2 / sum(pca_res$sdev^2) * 100, 1)  # % variance explained

samples <- rownames(pca_res$x)
meta <- data.frame(
  sample = samples,
  year = sub(".*_(\\d+)$", "\\1", samples)
)

n_pcs <- 10
pc_df <- data.frame(
  sample = samples,
  pca_res$x[, 1:n_pcs],
  year = meta$year
)

# ===============================
# 3) Clean metadata & create group variable
# ===============================
pc_df$group <- ifelse(grepl("^A\\d+", pc_df$sample), "A",
                      ifelse(grepl("^C\\d+", pc_df$sample), "C",
                             ifelse(grepl("^A>C", pc_df$sample), "A>C",
                                    ifelse(grepl("^C>A", pc_df$sample), "C>A", NA))))
pc_df$group <- factor(pc_df$group, levels = c("A>C", "C>A", "C", "A"))

# ===============================
# 3b) Create population variable for random effect
# ===============================
# Use replicate number as population
pc_df$population <- sub(".*?(\\d+)_.*", "\\1", pc_df$sample)
pc_df$population <- factor(pc_df$population)

# ===============================
# 4) Fit linear models for PCs with group as predictor (single timepoint)
# ===============================
pc_names <- paste0("PC", 1:n_pcs)
results_list <- list()

for (pc in pc_names) {
  formula <- as.formula(paste0(pc, " ~ group"))
  model <- lm(formula, data = pc_df)
  tidy_res <- broom::tidy(model, conf.int = TRUE)  # use broom for lm
  tidy_res$PC <- pc
  results_list[[pc]] <- tidy_res
}

results_all <- bind_rows(results_list)

# ===============================
# 5) Bonferroni correction and PC ordering
# ===============================
n_tests <- nrow(results_all)
results_clean <- results_all %>%
  mutate(
    p.adj = p.value * n_tests,
    p.adj = ifelse(p.adj > 1, 1, p.adj),
    PC_num = as.numeric(gsub("PC", "", PC)),
    PC = factor(PC, levels = paste0("PC", 1:n_pcs))
  ) %>%
  arrange(PC_num, term) %>%
  select(-PC_num)

sig_results <- results_clean %>% filter(p.adj < 0.005)


# ===============================
# 6) Export tables as CSV
# ===============================
write.csv(results_clean, "Endpoint_AllPCs_Bonferroni.csv", row.names = FALSE)
write.csv(sig_results, "Endpoint_SigPCs_Bonferroni.csv", row.names = FALSE)


# ===============================
# 7) PCA plot with colored groups and ellipses
# ===============================

# Filter only year 24
pc24 <- pc_df[pc_df$year == "24", ]

# Define colors
group_colors <- c(
  "A>C" = "#4DB8FF",   # light blue
  "C"   = "#0066CC",   # dark blue
  "C>A" = "#FF6666",   # light red
  "A"   = "#CC0000"    # dark red
)

# PCA plot
pca_plot <- ggplot(pc24, aes(x = PC1, y = PC2, color = group, fill = group, shape = group)) +
  geom_point(size = 6, alpha = 0.9) +  # bigger points
  stat_ellipse(aes(group = group), type = "norm", level = 0.95, alpha = 0.2, geom = "polygon") +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    color = "Group",
    fill = "Group",
    shape = "Group",
    title = ""
  ) +
  theme_bw(base_size = 24) +  # increase base font size
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.8),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 28, face = "bold"),     # axis titles
    axis.text = element_text(size = 24),                    # axis labels
    legend.title = element_text(size = 26, face = "bold"),  # legend titles
    legend.text = element_text(size = 24),                 # legend labels
    legend.key.size = unit(1.5, "cm"),                     # bigger legend keys
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)
  )

# Display plot
print(pca_plot)


# Save as SVG
ggsave(
  "Convergence_pc_scatter_ellipse.svg",
  plot = pca_plot,
  width = 14,
  height = 10
)

