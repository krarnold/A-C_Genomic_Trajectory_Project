library(dplyr)
library(stringr)

# ============================================================
# SEQUENCING STATISTICS SUMMARY SCRIPT
# Combines alignment logs (reads, aligned reads) with
# SNP table coverage for standard and deep sequencing
# ============================================================

# ---- PATHS ----
# Set data_dir to the path of your local data directory
data_dir        <- "."
snp_table_path  <- file.path(data_dir, "SNP_Table/arnold25_genomic_trajectory_SNP_table_counts.csv")
deep_raw_path   <- file.path(data_dir, "SNP_Table/ACO1_Deepseq_SNP_20_5000.csv")
output_dir      <- file.path(data_dir, "SNP_Table/AllLogs")

run_dirs <- list(
  Run1    = file.path(data_dir, "SNP_Table/AllLogs/AllLogs/logFiles"),
  Run2    = file.path(data_dir, "SNP_Table/AllLogs/AllLogs/logFilesRun2"),
  Run3    = file.path(data_dir, "SNP_Table/AllLogs/AllLogs/logFilesRun3"),
  DeepSeq = file.path(data_dir, "SNP_Table/AllLogs/AllLogs/ACO1_DeepSeq_Logs/ACO1_DeepSeq_Logs")
)

# ============================================================
# PART 1: PARSE ALIGNMENT LOGS
# ============================================================

parse_novoalign_log <- function(filepath, run_id) {
  lines <- readLines(filepath, warn = FALSE)
  filename <- basename(filepath)

  pop_id <- str_extract(filename, "(?<=arnold_).*(?=_novoalign_log)")
  if (is.na(pop_id)) pop_id <- gsub("\\.txt$", "", filename)

  read_seqs    <- lines[grep("Read Sequences:", lines)]
  unique_align <- lines[grep("Unique Alignment:", lines)]
  multi_mapped <- lines[grep("Multi Mapped:", lines)]
  no_mapping   <- lines[grep("No Mapping Found:", lines)]
  proper_pairs <- lines[grep("Proper Pairs:", lines)]
  homopolymer  <- lines[grep("Homopolymer Filter:", lines)]

  extract_num <- function(x) {
    if (length(x) == 0) return(NA)
    as.numeric(gsub(",", "", str_extract(x[1], "[0-9,]+")))
  }

  extract_pct <- function(x) {
    if (length(x) == 0) return(NA)
    as.numeric(str_extract(x[1], "[0-9]+\\.[0-9]+(?=%)"))
  }

  data.frame(
    run           = run_id,
    filename      = filename,
    population    = pop_id,
    read_seqs     = extract_num(read_seqs),
    unique_align  = extract_num(unique_align),
    pct_unique    = extract_pct(unique_align),
    multi_mapped  = extract_num(multi_mapped),
    pct_multi     = extract_pct(multi_mapped),
    no_mapping    = extract_num(no_mapping),
    pct_no_map    = extract_pct(no_mapping),
    proper_pairs  = extract_num(proper_pairs),
    pct_proper    = extract_pct(proper_pairs),
    homopolymer   = extract_num(homopolymer),
    stringsAsFactors = FALSE
  )
}

all_logs <- bind_rows(lapply(names(run_dirs), function(run_id) {
  dir_path <- run_dirs[[run_id]]
  files <- list.files(dir_path, pattern = "novoalign_log\\.txt$", full.names = TRUE)
  cat("Parsing", length(files), "novoalign log files from", run_id, "\n")
  if (length(files) == 0) {
    cat("  WARNING: No novoalign log files found in", dir_path, "\n")
    return(NULL)
  }
  bind_rows(lapply(files, parse_novoalign_log, run_id = run_id))
}))

cat("\nTotal log files parsed:", nrow(all_logs), "\n")
cat("Parsing failures (NA in read_seqs):", sum(is.na(all_logs$read_seqs)), "\n\n")

# Split standard and deep sequencing logs
standard_logs <- all_logs %>% filter(run %in% c("Run1", "Run2", "Run3"))
deepseq_logs  <- all_logs %>% filter(run == "DeepSeq")

# Aggregate standard sequencing across runs per P-code
standard_log_summary <- standard_logs %>%
  group_by(population) %>%
  summarise(
    n_runs            = n(),
    total_read_seqs   = sum(read_seqs, na.rm = TRUE),
    total_unique      = sum(unique_align, na.rm = TRUE),
    total_multi       = sum(multi_mapped, na.rm = TRUE),
    total_no_map      = sum(no_mapping, na.rm = TRUE),
    total_proper      = sum(proper_pairs, na.rm = TRUE),
    total_homopolymer = sum(homopolymer, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_unique  = round(total_unique / total_read_seqs * 100, 1),
    pct_multi   = round(total_multi  / total_read_seqs * 100, 1),
    pct_no_map  = round(total_no_map / total_read_seqs * 100, 1),
    pct_proper  = round(total_proper / (total_read_seqs / 2) * 100, 1),
    pct_homopol = round(total_homopolymer / total_read_seqs * 100, 1)
  )

# Overall means for standard sequencing logs
standard_log_overall <- standard_log_summary %>%
  summarise(
    mean_read_seqs   = mean(total_read_seqs, na.rm = TRUE),
    mean_unique      = mean(total_unique, na.rm = TRUE),
    mean_pct_unique  = mean(pct_unique, na.rm = TRUE),
    mean_pct_multi   = mean(pct_multi, na.rm = TRUE),
    mean_pct_no_map  = mean(pct_no_map, na.rm = TRUE),
    mean_pct_proper  = mean(pct_proper, na.rm = TRUE),
    mean_pct_homopol = mean(pct_homopol, na.rm = TRUE)
  )

# Deep sequencing log summary
deepseq_log_summary <- deepseq_logs %>%
  summarise(
    n_runs            = n(),
    total_read_seqs   = sum(read_seqs, na.rm = TRUE),
    total_unique      = sum(unique_align, na.rm = TRUE),
    total_multi       = sum(multi_mapped, na.rm = TRUE),
    total_no_map      = sum(no_mapping, na.rm = TRUE),
    total_proper      = sum(proper_pairs, na.rm = TRUE),
    total_homopolymer = sum(homopolymer, na.rm = TRUE)
  ) %>%
  mutate(
    pct_unique  = round(total_unique / total_read_seqs * 100, 1),
    pct_multi   = round(total_multi  / total_read_seqs * 100, 1),
    pct_no_map  = round(total_no_map / total_read_seqs * 100, 1),
    pct_proper  = round(total_proper / (total_read_seqs / 2) * 100, 1),
    pct_homopol = round(total_homopolymer / total_read_seqs * 100, 1)
  )

# ============================================================
# PART 2: SNP COVERAGE FROM RAW SNP TABLES
# ============================================================

cat("Reading standard sequencing SNP table...\n")
raw <- read.csv(snp_table_path, header = TRUE)

mj_cols <- grep("_mj$", names(raw), value = TRUE)

standard_coverage_per_sample <- bind_rows(lapply(mj_cols, function(mj_col) {
  mi_col <- sub("_mj$", "_mi", mj_col)
  sample_name <- sub("_mj$", "", mj_col)
  coverage <- raw[[mj_col]] + raw[[mi_col]]
  data.frame(
    sample     = sample_name,
    mean_cov   = round(mean(coverage, na.rm = TRUE), 1),
    median_cov = round(median(coverage, na.rm = TRUE), 1),
    min_cov    = round(min(coverage, na.rm = TRUE), 1),
    max_cov    = round(max(coverage, na.rm = TRUE), 1)
  )
}))

# Overall SNP coverage summary for standard sequencing
standard_cov_overall <- standard_coverage_per_sample %>%
  summarise(
    mean_snp_cov   = round(mean(mean_cov, na.rm = TRUE), 1),
    median_snp_cov = round(median(mean_cov, na.rm = TRUE), 1),
    min_snp_cov    = round(min(mean_cov, na.rm = TRUE), 1),
    max_snp_cov    = round(max(mean_cov, na.rm = TRUE), 1),
    total_samples  = n()
  )

cat("Reading deep sequencing SNP table...\n")
deep_raw <- read.csv(deep_raw_path, header = TRUE)

# Per sample deep sequencing coverage summary
deep_coverage_per_sample <- data.frame(
  sample     = "ACO1_DeepSeq",
  mean_cov   = round(mean(deep_raw$cov, na.rm = TRUE), 1),
  median_cov = round(median(deep_raw$cov, na.rm = TRUE), 1),
  min_cov    = round(min(deep_raw$cov, na.rm = TRUE), 1),
  max_cov    = round(max(deep_raw$cov, na.rm = TRUE), 1),
  p5_cov     = round(quantile(deep_raw$cov, 0.05, na.rm = TRUE), 1),
  p25_cov    = round(quantile(deep_raw$cov, 0.25, na.rm = TRUE), 1),
  p75_cov    = round(quantile(deep_raw$cov, 0.75, na.rm = TRUE), 1),
  p95_cov    = round(quantile(deep_raw$cov, 0.95, na.rm = TRUE), 1)
)

# Overall deep sequencing coverage summary
deep_cov_overall <- data.frame(
  mean_snp_cov   = round(mean(deep_raw$cov, na.rm = TRUE), 1),
  median_snp_cov = round(median(deep_raw$cov, na.rm = TRUE), 1),
  min_snp_cov    = round(min(deep_raw$cov, na.rm = TRUE), 1),
  max_snp_cov    = round(max(deep_raw$cov, na.rm = TRUE), 1),
  p5_snp_cov     = round(quantile(deep_raw$cov, 0.05, na.rm = TRUE), 1),
  p95_snp_cov    = round(quantile(deep_raw$cov, 0.95, na.rm = TRUE), 1),
  total_sites    = nrow(deep_raw)
)

# ============================================================
# PART 3: PRINT FULL SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("STANDARD SEQUENCING - PER SAMPLE SNP COVERAGE\n")
cat("============================================================\n")
print(as.data.frame(standard_coverage_per_sample), row.names = FALSE)

cat("\n============================================================\n")
cat("STANDARD SEQUENCING - OVERALL MEANS FROM ALIGNMENT LOGS\n")
cat("============================================================\n")
print(t(standard_log_overall))

cat("\n============================================================\n")
cat("STANDARD SEQUENCING - OVERALL SNP COVERAGE\n")
cat("============================================================\n")
print(as.data.frame(standard_cov_overall), row.names = FALSE)

cat("\n============================================================\n")
cat("DEEP SEQUENCING - ALIGNMENT LOG SUMMARY\n")
cat("============================================================\n")
print(t(deepseq_log_summary))

cat("\n============================================================\n")
cat("DEEP SEQUENCING - PER SAMPLE SNP COVERAGE\n")
cat("============================================================\n")
print(as.data.frame(deep_coverage_per_sample), row.names = FALSE)

cat("\n============================================================\n")
cat("DEEP SEQUENCING - OVERALL SNP COVERAGE\n")
cat("============================================================\n")
print(as.data.frame(deep_cov_overall), row.names = FALSE)

cat("\n============================================================\n")
cat("COMBINED SUMMARY TABLE FOR MANUSCRIPT\n")
cat("============================================================\n")
manuscript_table <- data.frame(
  Statistic = c(
    "Total reads (mean per sample)",
    "Uniquely aligned reads (mean per sample)",
    "Unique alignment rate",
    "Proper pair rate",
    "Homopolymer filtered",
    "Mean SNP coverage",
    "Median SNP coverage",
    "SNP coverage range (5th-95th percentile)",
    "SNP coverage range (absolute min-max)",
    "Total sites"
  ),
  Standard_Sequencing = c(
    formatC(round(standard_log_overall$mean_read_seqs), format = "d", big.mark = ","),
    formatC(round(standard_log_overall$mean_unique), format = "d", big.mark = ","),
    paste0(round(standard_log_overall$mean_pct_unique, 1), "%"),
    paste0(round(standard_log_overall$mean_pct_proper, 1), "%"),
    paste0(round(standard_log_overall$mean_pct_homopol, 1), "%"),
    paste0(standard_cov_overall$mean_snp_cov, "x"),
    paste0(standard_cov_overall$median_snp_cov, "x"),
    paste0(standard_cov_overall$min_snp_cov, "x - ",
           standard_cov_overall$max_snp_cov, "x"),
    paste0(standard_cov_overall$min_snp_cov, "x - ",
           standard_cov_overall$max_snp_cov, "x"),
    formatC(nrow(raw), format = "d", big.mark = ",")
  ),
  Deep_Sequencing = c(
    formatC(deepseq_log_summary$total_read_seqs, format = "d", big.mark = ","),
    formatC(deepseq_log_summary$total_unique, format = "d", big.mark = ","),
    paste0(deepseq_log_summary$pct_unique, "%"),
    paste0(deepseq_log_summary$pct_proper, "%"),
    paste0(deepseq_log_summary$pct_homopol, "%"),
    paste0(deep_cov_overall$mean_snp_cov, "x"),
    paste0(deep_cov_overall$median_snp_cov, "x"),
    paste0(deep_cov_overall$p5_snp_cov, "x - ",
           deep_cov_overall$p95_snp_cov, "x"),
    paste0(deep_cov_overall$min_snp_cov, "x - ",
           deep_cov_overall$max_snp_cov, "x"),
    formatC(deep_cov_overall$total_sites, format = "d", big.mark = ",")
  )
)
print(as.data.frame(manuscript_table), row.names = FALSE)

# ============================================================
# PART 4: SAVE ALL OUTPUTS
# ============================================================

write.csv(all_logs,
          file.path(output_dir, "all_novoalign_logs_parsed.csv"),
          row.names = FALSE)

write.csv(standard_log_summary,
          file.path(output_dir, "standard_sequencing_log_summary.csv"),
          row.names = FALSE)

write.csv(standard_coverage_per_sample,
          file.path(output_dir, "standard_snp_coverage_per_sample.csv"),
          row.names = FALSE)

write.csv(deep_coverage_per_sample,
          file.path(output_dir, "deepseq_snp_coverage_per_sample.csv"),
          row.names = FALSE)

write.csv(deep_cov_overall,
          file.path(output_dir, "deepseq_snp_coverage_overall.csv"),
          row.names = FALSE)

write.csv(manuscript_table,
          file.path(output_dir, "manuscript_sequencing_summary_table.csv"),
          row.names = FALSE)

cat("\nAll files saved to:", output_dir, "\n")
