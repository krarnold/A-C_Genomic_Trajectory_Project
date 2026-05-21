library(dplyr)

# ---- LOAD DATA ----
data    <- read.csv("SNP_Table/AC_trajectory_modified_snp_table.csv", header = TRUE)
sig_int <- read.csv("BetaBinomial_ME/Significant_SNPs_treatmentC_A_generation.csv", header = TRUE)
deep    <- read.csv("SNP_Table/ACO1_Deepseq_Modified.csv", header = TRUE)

# ---- A-TYPE FOUNDER COLUMNS (generation 0, timepoint _18) ----
atype_minor_cols <- grep("^(CACO[1-5]|CAO[1-5])_18_minor_count$",   names(data), value = TRUE)
atype_total_cols <- grep("^(CACO[1-5]|CAO[1-5])_18_total_coverage$", names(data), value = TRUE)

# ---- MEAN FOUNDER FREQUENCY ACROSS 10 A-TYPE POPULATIONS ----
minor_mat          <- as.matrix(data[, atype_minor_cols])
total_mat          <- as.matrix(data[, atype_total_cols])
atype_founder_freq <- rowMeans(minor_mat / total_mat, na.rm = TRUE)

# ---- PSEUDO-FIXED FLAG: minor count = 0 in ALL 10 A-type founders ----
pseudo_fixed_flag <- rowSums(minor_mat == 0) == 10

founders <- data.frame(
  chr                = data$chr,
  pos                = data$pos,
  atype_founder_freq = atype_founder_freq,
  pseudo_fixed_atype = pseudo_fixed_flag
)

# ---- MERGE WITH SIGNIFICANT INTERACTION SNPs ----
sig_founders <- sig_int %>%
  select(chr, pos) %>%
  inner_join(founders, by = c("chr", "pos"))

# ---- PART 1: FOUNDER FREQUENCY SUMMARY ----
cat("=== COMMENT 5 SUMMARY ===\n")
cat("Total significant interaction SNPs:        ", nrow(sig_founders), "\n")
cat("Pseudo-fixed in A-type founders:           ", sum(sig_founders$pseudo_fixed_atype), "\n")
cat("Fraction pseudo-fixed:                     ", round(mean(sig_founders$pseudo_fixed_atype), 4), "\n")
cat("Fraction at detectable frequency:          ", round(mean(!sig_founders$pseudo_fixed_atype), 4), "\n")
cat("Mean founder MAF (non-pseudo-fixed):       ",
    round(mean(sig_founders$atype_founder_freq[!sig_founders$pseudo_fixed_atype], na.rm = TRUE), 4), "\n\n")

# ---- PART 2: DEEP SEQ CROSS-REFERENCE (pseudo-fixed sig SNPs only) ----
deep_filtered      <- deep %>% filter(ACO1_minor_count >= 3)
deep_filtered$freq <- deep_filtered$ACO1_minor_count / deep_filtered$ACO1_total_coverage

pf_snps <- sig_founders %>% filter(pseudo_fixed_atype) %>% select(chr, pos)
pf_deep <- pf_snps %>% inner_join(deep_filtered, by = c("chr", "pos"))

cat("=== DEEP SEQ CROSS-REFERENCE (pseudo-fixed sig SNPs) ===\n")
cat("Pseudo-fixed sig SNPs detectable in deep seq (count >= 3):", nrow(pf_deep), "\n")
cat("Fraction of pseudo-fixed sig SNPs with deep seq coverage: ",
    round(nrow(pf_deep) / nrow(pf_snps), 4), "\n")
cat("Mean deep seq MAF:                                        ",
    round(mean(pf_deep$freq, na.rm = TRUE), 6), "\n")
cat("Median deep seq MAF:                                      ",
    round(median(pf_deep$freq, na.rm = TRUE), 6), "\n")

# ---- PART 3: REPLICATE CONSISTENCY OF PSEUDO-FIXED RESPONDERS ----
# For each pseudo-fixed significant SNP, count how many of the 10 A>C
# replicate populations showed a directional increase by generation 65.
# "Increase" = endpoint minor count / total > 0 (i.e., the allele became
# detectable in that replicate by the final timepoint).

ac_end_minor_cols <- grep("^(CACO[1-5]|CAO[1-5])_24_minor_count$",   names(data), value = TRUE)
ac_end_total_cols <- grep("^(CACO[1-5]|CAO[1-5])_24_total_coverage$", names(data), value = TRUE)

ac_end_freq_mat <- as.matrix(data[, ac_end_minor_cols]) /
  as.matrix(data[, ac_end_total_cols])

# A replicate "responded" if endpoint freq > 0 (allele emerged from zero)
ac_end_detected <- ac_end_freq_mat > 0

endpoint_df <- data.frame(
  chr             = data$chr,
  pos             = data$pos,
  n_reps_detected = rowSums(ac_end_detected, na.rm = TRUE)
)

pf_consistency <- sig_founders %>%
  filter(pseudo_fixed_atype) %>%
  select(chr, pos) %>%
  inner_join(endpoint_df, by = c("chr", "pos"))

cat("=== REPLICATE CONSISTENCY (pseudo-fixed sig SNPs) ===\n")
cat("Pseudo-fixed sig SNPs analyzed:                        ", nrow(pf_consistency), "\n\n")
cat("Distribution across replicates that detected allele:\n")
print(table(pf_consistency$n_reps_detected))
cat("\n")
cat("Fraction detected in >= 5 of 10 replicates:           ",
    round(mean(pf_consistency$n_reps_detected >= 5), 4), "\n")
cat("Fraction detected in >= 8 of 10 replicates:           ",
    round(mean(pf_consistency$n_reps_detected >= 8), 4), "\n")
cat("Fraction detected in all 10 replicates:               ",
    round(mean(pf_consistency$n_reps_detected == 10), 4), "\n")
cat("Mean number of replicates detecting allele:            ",
    round(mean(pf_consistency$n_reps_detected), 2), "\n")

# ---- PART 4: ENDPOINT FREQUENCIES OF PSEUDO-FIXED RESPONDERS AT GEN 65 ----
# For each pseudo-fixed significant SNP, calculate the mean, median, and max
# minor allele frequency reached across the 10 A>C replicate populations
# at the final timepoint (generation 65).

endpoint_freq_df <- data.frame(
  chr             = data$chr,
  pos             = data$pos,
  mean_end_freq   = rowMeans(ac_end_freq_mat, na.rm = TRUE),
  median_end_freq = apply(ac_end_freq_mat, 1, median, na.rm = TRUE),
  max_end_freq    = apply(ac_end_freq_mat, 1, max, na.rm = TRUE)
)

pf_endpoint <- sig_founders %>%
  filter(pseudo_fixed_atype) %>%
  select(chr, pos) %>%
  inner_join(endpoint_freq_df, by = c("chr", "pos"))

cat("=== ENDPOINT FREQUENCIES OF PSEUDO-FIXED RESPONDERS AT GEN 65 ===\n")
cat("N pseudo-fixed responders:               ", nrow(pf_endpoint), "\n")
cat("Mean endpoint MAF (mean across reps):    ", round(mean(pf_endpoint$mean_end_freq,   na.rm = TRUE), 4), "\n")
cat("Median endpoint MAF:                     ", round(median(pf_endpoint$mean_end_freq, na.rm = TRUE), 4), "\n")
cat("Max endpoint MAF observed:               ", round(max(pf_endpoint$max_end_freq,     na.rm = TRUE), 4), "\n")
cat("Distribution of mean endpoint MAF:\n")
print(quantile(pf_endpoint$mean_end_freq,
               probs = c(0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
               na.rm = TRUE))

# ---- SAVE SUMMARY OUTPUT ----

# Part 1: Founder frequency summary
founder_summary <- data.frame(
  Statistic = c(
    "Total significant interaction SNPs",
    "Pseudo-fixed in A-type founders (n)",
    "Fraction pseudo-fixed",
    "Fraction at detectable frequency",
    "Mean founder MAF (non-pseudo-fixed)"
  ),
  Value = c(
    nrow(sig_founders),
    sum(sig_founders$pseudo_fixed_atype),
    round(mean(sig_founders$pseudo_fixed_atype), 4),
    round(mean(!sig_founders$pseudo_fixed_atype), 4),
    round(mean(sig_founders$atype_founder_freq[!sig_founders$pseudo_fixed_atype], na.rm = TRUE), 4)
  )
)
write.csv(founder_summary, "FounderFreq_Summary.csv", row.names = FALSE)

# Part 2: Deep seq cross-reference summary
deepseq_summary <- data.frame(
  Statistic = c(
    "Pseudo-fixed sig SNPs detectable in deep seq (count >= 3)",
    "Fraction of pseudo-fixed sig SNPs with deep seq coverage",
    "Mean deep seq MAF",
    "Median deep seq MAF"
  ),
  Value = c(
    nrow(pf_deep),
    round(nrow(pf_deep) / nrow(pf_snps), 4),
    round(mean(pf_deep$freq, na.rm = TRUE), 6),
    round(median(pf_deep$freq, na.rm = TRUE), 6)
  )
)
write.csv(deepseq_summary, "FounderFreq_DeepSeq_Summary.csv", row.names = FALSE)

# Part 3: Replicate consistency summary
consistency_summary <- data.frame(
  Statistic = c(
    "Pseudo-fixed sig SNPs analyzed",
    "Fraction detected in >= 6 of 10 replicates",
    "Fraction detected in >= 8 of 10 replicates",
    "Fraction detected in all 10 replicates",
    "Mean replicates detecting allele"
  ),
  Value = c(
    nrow(pf_consistency),
    round(mean(pf_consistency$n_reps_detected >= 6), 4),
    round(mean(pf_consistency$n_reps_detected >= 8), 4),
    round(mean(pf_consistency$n_reps_detected == 10), 4),
    round(mean(pf_consistency$n_reps_detected), 2)
  )
)
write.csv(consistency_summary, "FounderFreq_Consistency_Summary.csv", row.names = FALSE)

# Replicate distribution table
rep_dist <- as.data.frame(table(pf_consistency$n_reps_detected))
colnames(rep_dist) <- c("n_replicates_detected", "n_sites")
write.csv(rep_dist, "FounderFreq_Replicate_Distribution.csv", row.names = FALSE)

# Part 4: Endpoint frequency summary
endpoint_summary <- data.frame(
  Statistic = c(
    "N pseudo-fixed responders",
    "Mean endpoint MAF (mean across replicates)",
    "Median endpoint MAF",
    "Max endpoint MAF observed",
    "10th percentile endpoint MAF",
    "25th percentile endpoint MAF",
    "75th percentile endpoint MAF",
    "90th percentile endpoint MAF"
  ),
  Value = c(
    nrow(pf_endpoint),
    round(mean(pf_endpoint$mean_end_freq,                                        na.rm = TRUE), 4),
    round(median(pf_endpoint$mean_end_freq,                                      na.rm = TRUE), 4),
    round(max(pf_endpoint$max_end_freq,                                          na.rm = TRUE), 4),
    round(quantile(pf_endpoint$mean_end_freq, probs = 0.10,                      na.rm = TRUE), 4),
    round(quantile(pf_endpoint$mean_end_freq, probs = 0.25,                      na.rm = TRUE), 4),
    round(quantile(pf_endpoint$mean_end_freq, probs = 0.75,                      na.rm = TRUE), 4),
    round(quantile(pf_endpoint$mean_end_freq, probs = 0.90,                      na.rm = TRUE), 4)
  )
)
write.csv(endpoint_summary, "PseudoFixed_Endpoint_Summary.csv", row.names = FALSE)
write.csv(pf_endpoint,      "PseudoFixed_Endpoint_Frequencies.csv", row.names = FALSE)

message("All summary files saved.")