library(dplyr)


#SNP table
data <- read.csv("SNP_Table/AC_trajectory_modified_snp_table.csv", header = TRUE)[,1:6]

#Deep seq table
deep <- read.csv("SNP_Table/ACO1_Deepseq_Modified.csv", header = T)


#A>C GLM Results
glm <- read.csv("BetaBinomial_ME/Raw_GLM_Betabi_Results.csv", header = TRUE)

# Apply Bonferroni correction to glm p-values
glm$bonferroni_p <- p.adjust(glm$p.value, method = "bonferroni")


# Subset SNPs significant at Bonferroni-corrected p < 0.005
sig_glm <- subset(glm, bonferroni_p < 0.005)

CACO1_sig <- merge(data, sig_glm, by = c("chr","pos"))[,1:6]

CACO1_sig_fixed <- subset(
  CACO1_sig,
  CACO1_18_minor_count == 0 | CACO1_18_minor_count == CACO1_18_total_coverage
)[,1:2]

#Subset fixed sites in deep seq and plot histogram
deep_pseudofixed <- merge(deep, CACO1_sig_fixed, by = c("chr","pos"))
deep_pseudofixed$freq <- deep_pseudofixed$ACO1_minor_count/deep_pseudofixed$ACO1_total_coverage

hist(deep_pseudofixed$freq,
     breaks = 50,
     col = "steelblue",
     border = "white",
     main = "Deep sequencing frequencies for pseudo-fixed sites",
     xlab = "Minor allele frequency (deep seq)",
     ylab = "Number of Sites")


#Get a random number of SNPs equal to psudeo fixed and overlap the histograms
# Sample random sites equal to number of pseudo-fixed sites
set.seed(123)
n_sites <- nrow(CACO1_sig_fixed)
deep_random <- deep[sample(nrow(deep), n_sites), ]
deep_random$freq <- deep_random$ACO1_minor_count / deep_random$ACO1_total_coverage

# Set x-axis limits for both datasets
xlims <- range(c(deep_random$freq, deep_pseudofixed$freq), na.rm = TRUE)

# Plot random subset (background, gray)
hist(deep_random$freq,
     breaks = 100,
     col = rgb(0.6, 0.6, 0.6, 0.6),
     border = "white",
     xlim = xlims,
     main = "",
     xlab = "Minor Allele Frequency",
     ylab = "Number of Sites")

# Overlay pseudo-fixed sites (foreground, blue)
hist(deep_pseudofixed$freq,
     breaks = 100,
     col = rgb(0, 0, 1, 0.5),
     border = "white",
     add = TRUE)

# Add legend
legend("topright",
       legend = c("Random subset", "Pseudo-fixed sites"),
       fill = c(rgb(0.6, 0.6, 0.6, 0.6), rgb(0, 0, 1, 0.5)),
       border = "white",
       bty = "n")


#Save Histogram
# Open scalable vector
svg("DeepSeq_PseudoFixed_Histogram.svg", width = 10, height = 6)  # 10x6 inches (16:9 approx)

# Set larger text sizes for axes and labels
par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.8, mar = c(5,5,4,2))

# Determine x-axis limits for consistency
xlims <- range(c(deep_random$freq, deep_pseudofixed$freq), na.rm = TRUE)

# Plot random subset first (background, gray)
hist(deep_random$freq, breaks = 100, col = rgb(0.6, 0.6, 0.6, 0.6), border = "white",
     xlim = xlims, main = "",
     xlab = "Minor Allele Frequency", ylab = "Number of Sites")

# Overlay pseudo-fixed sites (foreground, blue)
hist(deep_pseudofixed$freq, breaks = 100, col = rgb(0, 0, 1, 0.5), border = "white", add = TRUE)

# Add legend with larger text
legend("topright", legend = c("Random subset", "Pseudo-fixed sites"),
       fill = c(rgb(0.6, 0.6, 0.6, 0.6), rgb(0, 0, 1, 0.5)),
       border = "white", bty = "n", cex = 1.5)

# Close the device to write the file
dev.off()




#permutation test
# ---- PERMUTATION TEST: are pseudo-fixed means lower than random? ----
set.seed(123)
obs_mean <- mean(deep_pseudofixed$freq, na.rm = TRUE)
k <- nrow(deep_pseudofixed)
n_perm <- 10000

perm_means <- replicate(n_perm, {
  idx <- sample(nrow(deep), k)
  mean(deep$ACO1_minor_count[idx] / deep$ACO1_total_coverage[idx], na.rm = TRUE)
})

p_perm <- mean(perm_means <= obs_mean)
mean_null <- mean(perm_means)
sd_null <- sd(perm_means)
z_score <- (obs_mean - mean_null) / sd_null


perm_summary <- data.frame(
  Statistic = c("Observed mean (pseudo-fixed)",
                "Null mean (random subsets)",
                "Null SD",
                "Z-score",
                "Permutation p-value",
                "Number of permutations",
                "Number of pseudo-fixed sites"),
  Value = c(round(obs_mean, 6),
            round(mean_null, 6),
            round(sd_null, 6),
            round(z_score, 2),
            formatC(p_perm, format = "e", digits = 2),
            n_perm,
            k)
)

print(perm_summary, row.names = FALSE)




