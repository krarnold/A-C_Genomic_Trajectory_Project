# Load required packages
library(lme4)
library(future.apply)

# Read in data
data <- read.csv("Done/AC_trajectory_modified_snp_table.csv", header = TRUE)


# New A vs Old A Comparison


# Identify relevant columns
cols_A <- grep("^(NACO|ANCO|ACO|AO)[0-9]*_24(_minor_count|_total_coverage)$", names(data), value = TRUE)
A_data <- data[, c("chr", "pos", cols_A)]
rm(data)

# Filter to retain polymorphic sites
sum_maf_A <- rowSums(A_data[, seq(3, ncol(A_data), by = 2)] / A_data[, seq(4, ncol(A_data), by = 2)])
A_data <- cbind(A_data, sum_maf_A)
A_data <- subset(A_data, sum_maf_A > 0 & sum_maf_A < 20)

# Extract count matrices
A_coverage <- A_data[, grep("_total_coverage$", names(A_data), value = TRUE)]
A_minor    <- A_data[, grep("_minor_count$", names(A_data), value = TRUE)]
A_major    <- A_coverage - A_minor

# Metadata
pop_A <- colnames(A_minor)
treatment_A <- c(rep("old", 10), rep("new", 10))

# Robust GLMM test function for A
test_glmer_A <- function(i) {
  A1 <- t(A_minor[i, ])
  A2 <- t(A_major[i, ])
  
  result <- tryCatch({
    fit1 <- suppressMessages(glmer(cbind(A1, A2) ~ treatment_A + (1 | pop_A), family = "binomial",
                                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))))
    fit0 <- suppressMessages(glmer(cbind(A1, A2) ~ (1 | pop_A), family = "binomial",
                                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))))
    LRT <- anova(fit1, fit0)
    LRT$`Pr(>Chisq)`[2]
  }, error = function(e) {
    NA
  })
  
  return(result)
}


# Run A comparison in parallel
plan(multisession, workers = 10)
pval_A <- future_sapply(1:nrow(A_minor), test_glmer_A)
A_glm <- cbind(A_data[, 1:2], pval = pval_A)


#Save output
write.table(A_glm, "GLMM_AvA.txt", quote = FALSE, row.names = FALSE, sep = "\t")


