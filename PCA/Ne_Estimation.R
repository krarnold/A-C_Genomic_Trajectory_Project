library(poolSeq)

data <- read.csv("SNP_Table/AC_trajectory_modified_snp_table.csv", header = TRUE)

# Replace minor counts with minor allele frequencies
minor_cols <- grep("_minor_count$", names(data), value = TRUE)
coverage_cols <- gsub("_minor_count$", "_total_coverage", minor_cols)

for (i in seq_along(minor_cols)) {
  data[[minor_cols[i]]] <- data[[minor_cols[i]]] / data[[coverage_cols[i]]]
}
names(data)[names(data) %in% minor_cols] <- gsub("_minor_count$", "_minor_freq", minor_cols)

# Define matched pairs: p0 = timepoint 18, pt = timepoint 24
pair_map <- c(
  CACO1 = "ACO1", CACO2 = "ACO2", CACO3 = "ACO3", CACO4 = "ACO4", CACO5 = "ACO5",
  CAO1  = "AO1",  CAO2  = "AO2",  CAO3  = "AO3",  CAO4  = "AO4",  CAO5  = "AO5"
)

ne_results <- list()

for (p0_name in names(pair_map)) {
  pt_name <- pair_map[[p0_name]]
  
  p0_col   <- paste0(p0_name, "_18_minor_freq")
  cov0_col <- paste0(p0_name, "_18_total_coverage")
  pt_col   <- paste0(pt_name, "_24_minor_freq")
  covt_col <- paste0(pt_name, "_24_total_coverage")
  
  if (all(c(p0_col, cov0_col, pt_col, covt_col) %in% names(data))) {
    ne_result <- estimateNe(
      p0       = data[[p0_col]],
      pt       = data[[pt_col]],
      cov0     = data[[cov0_col]],
      covt     = data[[covt_col]],
      t        = 219,
      ploidy   = 2,
      method   = "P.planII",
      poolSize = c(200, 200)
    )
    ne_results[[p0_name]] <- ne_result
  } else {
    warning(paste("Missing columns for pair:", p0_name, "->", pt_name))
  }
}

Atype_founder_Ne <- data.frame(
  population_t0 = names(ne_results),
  population_t1 = pair_map[names(ne_results)],
  Ne = as.numeric(ne_results)
)

print(Atype_founder_Ne)
write.csv(Atype_founder_Ne, "Atype_founder_Ne.csv", row.names = FALSE)
