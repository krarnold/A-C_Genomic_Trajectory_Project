# ===============================================================
#  Parallelism LOO Results - t-tests + Visualization (median shift)
# ===============================================================

setwd("~/Library/CloudStorage/Dropbox/AC Trajectory Genomic Data/LOO_palrallemism_test")

library(dplyr)
library(broom)
library(readr)
library(ggplot2)

# --- Load LOO median results ---
AC_median <- read.csv("AC_LOO_medians.csv", header = TRUE)
CA_median <- read.csv("CA_LOO_medians.csv", header = TRUE)

# ---------------------------------------------------------------
# 1. Genome-level t-tests
# ---------------------------------------------------------------
AC_genome_t <- AC_median %>%
  filter(chr == "genome") %>%
  t.test(median_delta ~ sig, data = .) %>%
  tidy() %>%
  mutate(treatment = "A>C", level = "genome", chr = "Genome")

CA_genome_t <- CA_median %>%
  filter(chr == "genome") %>%
  t.test(median_delta ~ sig, data = .) %>%
  tidy() %>%
  mutate(treatment = "C>A", level = "genome", chr = "Genome")

# ---------------------------------------------------------------
# 2. Chromosome-level t-tests
# ---------------------------------------------------------------
AC_chr_t <- AC_median %>%
  filter(chr != "genome") %>%
  group_by(chr) %>%
  group_modify(~ tidy(t.test(median_delta ~ sig, data = .x))) %>%
  mutate(treatment = "A>C", level = "chromosome")

CA_chr_t <- CA_median %>%
  filter(chr != "genome") %>%
  group_by(chr) %>%
  group_modify(~ tidy(t.test(median_delta ~ sig, data = .x))) %>%
  mutate(treatment = "C>A", level = "chromosome")

# ---------------------------------------------------------------
# 3. Combine t-test results
# ---------------------------------------------------------------
t_results <- bind_rows(AC_genome_t, CA_genome_t, AC_chr_t, CA_chr_t) %>%
  select(treatment, level, chr, estimate1, estimate2, statistic, p.value, conf.low, conf.high)

t_results$fdr <- p.adjust(t_results$p.value, method = "fdr")

write.csv(t_results, "LOO_ttest_results.csv", row.names = FALSE)

# ---------------------------------------------------------------
# 4. Prepare data for plotting
# ---------------------------------------------------------------
AC_median$Treatment <- "A>C"
CA_median$Treatment <- "C>A"

combined <- bind_rows(AC_median, CA_median) %>%
  mutate(
    chr = ifelse(chr == "genome", "Genome", chr),
    chr = factor(chr, levels = c("Genome", "2L", "2R", "3L", "3R", "X")),
    sig = factor(sig, levels = c("control", "target")),
    category = case_when(
      sig == "control" & Treatment == "A>C" ~ "A>C matched controls",
      sig == "control" & Treatment == "C>A" ~ "C>A matched controls",
      sig == "target"  & Treatment == "A>C" ~ "A>C targets",
      sig == "target"  & Treatment == "C>A" ~ "C>A targets"
    ),
    category = factor(
      category,
      levels = c("A>C targets", "A>C matched controls", "C>A targets", "C>A matched controls")
    )
  )

# Shared y-axis limits
y_range  <- range(combined$median_delta, na.rm = TRUE)
y_pad    <- diff(y_range) * 0.05
y_limits <- c(y_range[1] - y_pad, y_range[2] + y_pad)

# ---------------------------------------------------------------
# Updated colors and shapes (Option 1: desaturated matched colors)
# ---------------------------------------------------------------
col_map <- c(
  "A>C targets"           = "#3B82F6",  # bright blue
  "A>C matched controls"  = "#9EC5E8",  # light blue
  "C>A targets"           = "#EF4444",  # bright red
  "C>A matched controls"  = "#F5A3A3"   # light red
)

shape_map <- c(
  "A>C targets"           = 16,
  "A>C matched controls"  = 16,
  "C>A targets"           = 16,
  "C>A matched controls"  = 16
)

# ---------------------------------------------------------------
# 5. Plot (median shifts; dashed baseline at 0)
# ---------------------------------------------------------------
p <- ggplot(combined, aes(x = Treatment, y = median_delta, color = category, shape = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.5) +
  geom_point(
    position = position_jitter(width = 0.15),
    size = 2.5,
    alpha = 0.9
  ) +
  facet_wrap(~ chr, ncol = 3, scales = "fixed") +
  scale_color_manual(values = col_map, name = "Category") +
  scale_shape_manual(values = shape_map, name = "Category") +
  coord_cartesian(ylim = y_limits) +
  labs(
    x = NULL,
    y = "Median change in allele frequency"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 13),
    panel.grid.major.x = element_blank(),
    panel.spacing      = unit(0.8, "lines"),
    legend.position    = "right",
    legend.title       = element_text(face = "bold"),
    legend.text        = element_text(size = 11),
    plot.title         = element_blank()
  )

ggsave("LOO_median_parallelism.svg", p, width = 10, height = 7, device = "svg")


