
# Load required packages
library(lme4)
library(future.apply)

# Read in data
data <- read.csv("SNP_Table/AC_trajectory_modified_snp_table.csv", header = TRUE)


# New C vs Old C Comparison


# Identify relevant columns
cols_C <- grep("^(CACO|CAO|CO|NCO)[0-9]*_24(_minor_count|_total_coverage)$", names(data), value = TRUE)
C_data <- data[, c("chr", "pos", cols_C)]
rm(data)

# Filter to retain polymorphic sites
sum_maf_C <- rowSums(C_data[, seq(3, ncol(C_data), by = 2)] / C_data[, seq(4, ncol(C_data), by = 2)])
C_data <- cbind(C_data, sum_maf_C)
C_data <- subset(C_data, sum_maf_C > 0 & sum_maf_C < 20)

# Extract count matrices
C_coverage <- C_data[, grep("_total_coverage$", names(C_data), value = TRUE)]
C_minor    <- C_data[, grep("_minor_count$", names(C_data), value = TRUE)]
C_major    <- C_coverage - C_minor

# Metadata
pop_C <- colnames(C_minor)
treatment_C <- c(rep("old", 10), rep("new", 10))

# Robust GLMM test function for C
test_glmer_C <- function(i) {
  A1 <- t(C_minor[i, ])
  A2 <- t(C_major[i, ])
  
  result <- tryCatch({
    fit1 <- suppressMessages(glmer(cbind(A1, A2) ~ treatment_C + (1 | pop_C), family = "binomial",
                                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))))
    fit0 <- suppressMessages(glmer(cbind(A1, A2) ~ (1 | pop_C), family = "binomial",
                                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))))
    LRT <- anova(fit1, fit0)
    LRT$`Pr(>Chisq)`[2]
  }, error = function(e) {
    NA
  })
  
  return(result)
}


# Run C comparison in parallel
plan(multisession, workers = 10)
pval_C <- future_sapply(1:nrow(C_minor), test_glmer_C)
C_glm <- cbind(C_data[, 1:2], pval = pval_C)

#Save output
write.table(C_glm, "GLMM_CvC.txt", quote = FALSE, row.names = FALSE, sep = "\t")


