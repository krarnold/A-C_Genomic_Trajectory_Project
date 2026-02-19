# --- Libraries ---
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(scales)
library(ggtext)   



#Read and clean SNP data

data <- read.csv("SNP_Table/AC_trajectory_modified_snp_table.csv", header = TRUE)

# Remove columns corresponding to 24th timepoint if present
cols_to_remove <- grep("^(ACO|AO|CO|NCO)[0-9]_24", names(data), value = TRUE)
data <- data[ , !names(data) %in% cols_to_remove ]

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

# Calculate SNP frequencies
samples <- gsub("_minor_count$", "", grep("_minor_count$", names(data), value = TRUE))
freq <- data[, c("chr", "pos")]

for (s in samples) {
  k <- data[[paste0(s, "_minor_count")]]
  n <- data[[paste0(s, "_total_coverage")]]
  freq[[s]] <- k / n
}


# PCA on SNP frequencies
freq_mat <- as.matrix(freq[, -c(1,2)])
freq_mat <- t(freq_mat)
pca_res <- prcomp(freq_mat, center = FALSE, scale. = TRUE)
var_exp <- round((pca_res$sdev)^2 / sum(pca_res$sdev^2) * 100, 1)  # % variance explained

samples <- rownames(pca_res$x)
meta <- data.frame(
  sample = samples,
  traj = sub("^([^_]+).*", "\\1", samples),
  year = sub(".*_(\\d+)$", "\\1", samples)
)

n_pcs <- 10
pc_df <- data.frame(
  sample = samples,
  pca_res$x[, 1:n_pcs],
  traj = meta$traj,
  year = meta$year
)


# Clean metadata & create continuous generation variable
pc_df$group <- sub("^(A>C|C>A).*", "\\1", pc_df$sample)
pc_df$year <- sub(".*_(\\d+)$", "\\1", pc_df$sample)
pc_df$group <- factor(pc_df$group, levels = c("A>C", "C>A"))

# Continuous generations for modeling
pc_df <- pc_df %>%
  mutate(
    generation = case_when(
      group == "C>A" & year == "18" ~ 0,
      group == "C>A" & year == "19" ~ 16,
      group == "C>A" & year == "20" ~ 29,
      group == "C>A" & year == "24" ~ 182,
      group == "A>C" & year == "18" ~ 0,
      group == "A>C" & year == "19" ~ 6,
      group == "A>C" & year == "20" ~ 10,
      group == "A>C" & year == "24" ~ 65,
      TRUE ~ NA_real_
    ),
    population = sub(".*([0-9]+)$", "\\1", sample)
  ) %>%
  mutate(across(c(group, population), as.factor))


# Fit mixed-effects models for PCs with generation as continuous variable

pc_names <- paste0("PC", 1:n_pcs)
results_list <- list()

for (pc in pc_names) {
  formula <- as.formula(paste0(pc, " ~ group * generation + (1|population)"))
  model <- lmer(formula, data = pc_df)
  tidy_res <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  tidy_res$PC <- pc
  results_list[[pc]] <- tidy_res
}

results_all <- bind_rows(results_list)


# Bonferroni correction and PC ordering
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


# Export tables as CSV

write.csv(results_clean, "MixedEffectModel_AllPCs_Bonferroni.csv", row.names = FALSE)
write.csv(sig_results, "MixedEffectModel_SignificantEffects_Bonferroni.csv", row.names = FALSE)


# Flexible faceted line plot for PCA

pc_plot_df <- pc_df %>%
  mutate(
    traj = sub("^([AC]>[AC])\\d+_\\d+$", "\\1", sample),
    replicate = sub("^[AC]>[AC](\\d+)_\\d+$", "\\1", sample),
    year = sub(".*_(\\d+)$", "\\1", sample)
  )

pc_plot_df$year <- factor(pc_plot_df$year, levels = c("18", "19", "20", "24"))
pc_plot_df$traj <- factor(pc_plot_df$traj, levels = c("A>C", "C>A"))

# --- Flexible PC selection ---
available_pcs <- grep("^PC[0-9]+$", names(pc_df), value = TRUE)
pcs_to_plot <- c(1, 2)       # default: PC1 & PC2
#pcs_to_plot <- 1:4          
# pcs_to_plot <- c(5, 9)      

pcs_to_plot <- paste0("PC", pcs_to_plot)
pcs_to_plot <- pcs_to_plot[pcs_to_plot %in% available_pcs]  # ensure only available PCs

pc_long <- pc_plot_df %>%
  pivot_longer(cols = all_of(pcs_to_plot),
               names_to = "PC",
               values_to = "value")

# Variance explained for selected PCs
var_exp_plot <- round((pca_res$sdev[as.numeric(sub("PC", "", pcs_to_plot))]^2 / sum(pca_res$sdev^2)) * 100, 1)
names(var_exp_plot) <- pcs_to_plot

# generation mapping
ac_gen <- c("18" = 0,  "19" = 6,  "20" = 10, "24" = 65)
ca_gen <- c("18" = 0,  "19" = 16, "20" = 29, "24" = 182)

pc_plot <- ggplot(
  pc_long,
  aes(x = year, y = value, color = traj,
      group = interaction(traj, replicate))
) +
  geom_line(linewidth = 1.2, alpha = 0.7) +
  geom_point(size = 4) +
  scale_color_manual(values = c("A>C" = "#1f78b4", "C>A" = "#e31a1c")) +
  facet_wrap(
    ~PC, ncol = 2, scales = "free_y",
    labeller = labeller(PC = function(x) {
      paste0(x, " (", var_exp_plot[x], "%)")
    })
  ) +
  scale_x_discrete(
    labels = function(yrs) {
      sapply(yrs, function(y) {
        paste0(
          "<span style='color:#1f78b4'>", ac_gen[[y]], "</span><br>",
          "<span style='color:#e31a1c'>", ca_gen[[y]], "</span>"
        )
      })
    }
  ) +
  labs(
    x = "Generation",
    y = "PC value",
    color = "Trajectory"
  ) +
  theme_classic(base_size = 22) +   # << Increase all default text
  theme(
    strip.text = element_text(size = 24, face = "bold"),       # facet labels
    axis.text.x = ggtext::element_markdown(size = 20),         # colored x labels
    axis.text.y = element_text(size = 20),
    axis.title  = element_text(size = 24, face = "bold"),
    legend.title = element_text(size = 22),
    legend.text  = element_text(size = 20),
    legend.key.size = unit(1.2, "cm"),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(0.3, "cm")
  )


print(pc_plot)


# Save as SVG
ggsave(
  "PC_Plot.svg",
  plot = pc_plot,
  width = 14,
  height = 6
)


library(dplyr)
library(readr)

