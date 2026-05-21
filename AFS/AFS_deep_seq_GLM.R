library(dplyr)

# ---- LOAD DATA ----
data <- read.csv("SNP_Table/AC_trajectory_modified_snp_table.csv", header = TRUE)[,1:6]
deep <- read.csv("SNP_Table/ACO1_Deepseq_Modified.csv", header = T)

# ---- LOAD FDR PERMUTATION SIGNIFICANT SNPs (treatmentC>A:generation interaction) ----
# Significant_SNPs_treatmentC_A_generation.csv contains SNPs that passed the
# permutation-based FDR threshold for the interaction term (treatmentC>A:generation).
sig_glm <- read.csv("BetaBinomial_ME/Significant_SNPs_treatmentC_A_generation.csv", header = TRUE)

# ---- MERGE SIGNIFICANT SNPs WITH SNP TABLE ----
CACO1_sig <- merge(data, sig_glm, by = c("chr","pos"))[,1:6]

# ---- IDENTIFY PSEUDO-FIXED SITES ----
CACO1_sig_fixed <- subset(
  CACO1_sig,
  CACO1_18_minor_count == 0 | CACO1_18_minor_count == CACO1_18_total_coverage
)[,1:2]

cat("Total pseudo-fixed sites:", nrow(CACO1_sig_fixed), "\n")

# ---- MERGE WITH DEEP SEQUENCING DATA ----
deep_pseudofixed <- merge(deep, CACO1_sig_fixed, by = c("chr","pos"))
deep_pseudofixed$freq <- deep_pseudofixed$ACO1_minor_count / deep_pseudofixed$ACO1_total_coverage

cat("Pseudo-fixed sites with deep seq data (unfiltered):", nrow(deep_pseudofixed), "\n")

# ---- COUNT DISTRIBUTION AT PSEUDO-FIXED SITES (diagnostic) ----
cat("\nCount distribution at pseudo-fixed sites (top 10 count values):\n")
count_dist <- table(deep_pseudofixed$ACO1_minor_count)
print(head(count_dist, 10))

# ---- APPLY MINIMUM MINOR ALLELE COUNT THRESHOLD OF 3 ----
# Motivated by:
# (1) Illumina error rate ~0.001 predicts ~0.75 spurious reads/site at 767x depth
# (2) Sharp spike at count=1 and count=2 consistent with error distribution
# (3) Count of 3 corresponds to min detectable freq of ~0.00391 (3/767),
#     modestly above the detection floor of 1/341 = ~0.00293
deep_pseudofixed_filtered <- subset(deep_pseudofixed, ACO1_minor_count >= 3)

cat("\nPseudo-fixed sites AFTER count threshold of 3:", nrow(deep_pseudofixed_filtered), "\n")
cat("Mean MAF (filtered):", round(mean(deep_pseudofixed_filtered$freq, na.rm=TRUE), 6), "\n")
cat("Median MAF (filtered):", round(median(deep_pseudofixed_filtered$freq, na.rm=TRUE), 6), "\n")

# ---- RANDOM SUBSET FOR COMPARISON (same size as filtered set) ----
set.seed(123)
n_sites <- nrow(deep_pseudofixed_filtered)
deep_random <- deep[sample(nrow(deep), n_sites), ]
deep_random$freq <- deep_random$ACO1_minor_count / deep_random$ACO1_total_coverage

cat("\nRandom subset size:", nrow(deep_random), "\n")
cat("Random subset mean MAF:", round(mean(deep_random$freq, na.rm=TRUE), 6), "\n")

# ---- PERMUTATION TEST (filtered set) ----
set.seed(123)
obs_mean <- mean(deep_pseudofixed_filtered$freq, na.rm = TRUE)
k         <- nrow(deep_pseudofixed_filtered)
n_perm    <- 10000

perm_means <- replicate(n_perm, {
  idx <- sample(nrow(deep), k)
  mean(deep$ACO1_minor_count[idx] / deep$ACO1_total_coverage[idx], na.rm = TRUE)
})

p_perm    <- mean(perm_means <= obs_mean)
mean_null <- mean(perm_means)
sd_null   <- sd(perm_means)
z_score   <- (obs_mean - mean_null) / sd_null

cat("\n---- PERMUTATION TEST RESULTS (count >= 3 threshold) ----\n")
perm_summary <- data.frame(
  Statistic = c(
    "Observed mean MAF (pseudo-fixed, filtered)",
    "Median MAF (pseudo-fixed, filtered)",
    "Null mean (random subsets)",
    "Null SD",
    "Z-score",
    "Permutation p-value",
    "Number of permutations",
    "Sites before filter",
    "Sites after filter (count >= 3)",
    "Sites excluded by filter"
  ),
  Value = c(
    round(obs_mean, 6),
    round(median(deep_pseudofixed_filtered$freq, na.rm=TRUE), 6),
    round(mean_null, 6),
    round(sd_null, 6),
    round(z_score, 2),
    formatC(p_perm, format = "e", digits = 2),
    n_perm,
    nrow(deep_pseudofixed),
    nrow(deep_pseudofixed_filtered),
    nrow(deep_pseudofixed) - nrow(deep_pseudofixed_filtered)
  )
)
print(perm_summary, row.names = FALSE)

# ---- FRACTION OF SITES EXCLUDED (sanity check) ----
frac_retained <- nrow(deep_pseudofixed_filtered) / nrow(deep_pseudofixed)
cat("\nFraction of pseudo-fixed sites retained after filter:", round(frac_retained, 4), "\n")
cat("Fraction excluded (likely error-dominated):", round(1 - frac_retained, 4), "\n")

# ---- DETECTION LIMIT CALCULATIONS ----
cat("\n---- DETECTION LIMITS ----\n")
set.seed(42)
x <- sapply(1:10000, function(i) length(unique(sample(1:400, 90,  replace=TRUE))))
y <- sapply(1:10000, function(i) length(unique(sample(1:400, 767, replace=TRUE))))
cat("Mean unique alleles at standard depth (~90x):", round(mean(x), 1), "\n")
cat("Mean unique alleles at deep seq depth (~767x):", round(mean(y), 1), "\n")
cat("Min detectable freq at standard depth:", round(1/floor(mean(x)), 6), "\n")
cat("Min detectable freq at deep seq depth:", round(1/floor(mean(y)), 6), "\n")
cat("Min detectable freq at count >= 3 threshold:", round(3/767, 6), "\n")

# ---- SAVE HISTOGRAM ----
tiff("DeepSeq_PseudoFixed_Histogram.tiff", width = 10, height = 6, units = "in", res = 900, compression = "lzw")

par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.8, mar = c(5,5,4,2))

xlims <- range(c(deep_random$freq, deep_pseudofixed_filtered$freq), na.rm = TRUE)

shared_breaks <- seq(xlims[1], xlims[2], length.out = 101)

h1 <- hist(deep_random$freq, breaks = shared_breaks, plot = FALSE)
h2 <- hist(deep_pseudofixed_filtered$freq, breaks = shared_breaks, plot = FALSE)
ylims <- c(0, max(c(h1$counts, h2$counts)) * 1.05)

hist(deep_random$freq, breaks = shared_breaks, col = rgb(0.6, 0.6, 0.6, 0.6), border = "white",
     xlim = xlims, ylim = ylims, main = "",
     xlab = "Minor Allele Frequency", ylab = "Number of Sites")

hist(deep_pseudofixed_filtered$freq, breaks = shared_breaks, col = rgb(0, 0, 1, 0.5), border = "white", add = TRUE)

legend("topright", legend = c("Random subset", "Pseudo-fixed sites"),
       fill = c(rgb(0.6, 0.6, 0.6, 0.6), rgb(0, 0, 1, 0.5)),
       border = "white", bty = "n", cex = 1.5)

dev.off()

# ---- SAVE DIAGNOSTIC OUTPUTS ----

# Count distribution at pseudo-fixed sites (motivates count >= 3 threshold)
count_dist_df <- as.data.frame(count_dist)
colnames(count_dist_df) <- c("minor_allele_count", "n_sites")
write.csv(count_dist_df, "DeepSeq_PseudoFixed_CountDistribution.csv", row.names = FALSE)

# Permutation test summary
write.csv(perm_summary, "DeepSeq_PseudoFixed_PermutationSummary.csv", row.names = FALSE)