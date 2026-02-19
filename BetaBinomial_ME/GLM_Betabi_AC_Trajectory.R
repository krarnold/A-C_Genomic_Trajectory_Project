#library
library(tidyr)
library(dplyr)
library(stringr)
library(glmmTMB)
library(broom.mixed)
library(furrr)
library(future)
library(purrr)


data <- read.csv("~/Library/CloudStorage/Dropbox/AC Trajectory Genomic Data/SNP_Table/AC_trajectory_modified_snp_table.csv", header = TRUE)

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

#Repalce total coverage with major count
# find all columns ending with "_total_coverage"
total_cols <- grep("_total_coverage$", names(data), value = TRUE)

for (col in total_cols) {
  # find the corresponding minor_count column
  minor_col <- sub("_total_coverage$", "_minor_count", col)
  
  # define the new column name (replace suffix)
  new_colname <- sub("_total_coverage$", "_major_count", col)
  
  # compute major count = total - minor
  data[[col]] <- data[[col]] - data[[minor_col]]
  
  # rename the column
  names(data)[names(data) == col] <- new_colname
}

#change year designations to generations
# rename for A>C samples
names(data) <- sub("A>C([0-9]+)_18", "A>C\\1_0",  names(data))
names(data) <- sub("A>C([0-9]+)_19", "A>C\\1_6",  names(data))
names(data) <- sub("A>C([0-9]+)_20", "A>C\\1_10", names(data))
names(data) <- sub("A>C([0-9]+)_24", "A>C\\1_65", names(data))

# rename for C>A samples
names(data) <- sub("C>A([0-9]+)_18", "C>A\\1_0",   names(data))
names(data) <- sub("C>A([0-9]+)_19", "C>A\\1_16",  names(data))
names(data) <- sub("C>A([0-9]+)_20", "C>A\\1_29",  names(data))
names(data) <- sub("C>A([0-9]+)_24", "C>A\\1_182", names(data))


#make table long format and remove old data
data_long <- data %>%
  pivot_longer(
    cols = -c(chr, pos, Major_Nuc, Minor_Nuc),
    names_to = c("popdir", "replicate", "generation", ".value"),
    names_pattern = "^(A>C|C>A)(\\d+)_([0-9]+)_(minor|major)_count$"
  ) %>%
  mutate(
    population = paste0(popdir, replicate, "_", generation)
  ) %>%
  select(
    population, chr, pos, Major_Nuc, Minor_Nuc,
    Minor_Count = minor,
    Major_Count = major
  )

rm(data)

#add geeneration and treatment columns, clean up population column
data_long <- data_long %>%
  mutate(
    # extract treatment (e.g., "A>C" or "C>A")
    treatment = str_extract(population, "^[AC]?>[AC]?"),
    
    # extract generation (the digits after the final underscore)
    generation = as.numeric(str_extract(population, "(?<=_)\\d+$"))
  ) %>%
  select(population, treatment, generation, chr, pos, Major_Nuc, Minor_Nuc, Minor_Count, Major_Count)


data_long <- data_long %>%
  mutate(
    population = str_replace(population, "_\\d+$", "")
  )


#add a 1 to 0 counts
data_long <- data_long %>%
  mutate(
    Minor_Count = ifelse(Minor_Count == 0, Minor_Count + 1, Minor_Count),
    Major_Count = ifelse(Major_Count == 0, Major_Count + 1, Major_Count)
  )


#run test
# STEP 1: Set up parallel processing (8 cores)
plan(multisession, workers = 10)

# STEP 2: Prepare the data
# Ensure variables are the correct type for modeling
data_long <- data_long %>%
  mutate(
    treatment = factor(treatment),
    population = factor(population),
    generation = as.numeric(generation)
  )

write.csv(data_long, "Long_format_SNP_Data.csv", row.names = FALSE)

# STEP 3: Split data by SNP
# Creates a list where each element is the data for a single SNP (chr + pos) and remove data_long
data_split <- data_long %>%
  group_by(chr, pos) %>%
  group_split()

rm(data_long)

# STEP 4: Define the beta-binomial GLMM function
# Each worker will receive one small SNP subset at a time
fit_glmm <- function(snp_data) {
  if (n_distinct(snp_data$treatment) < 2 ||
      n_distinct(snp_data$generation) < 2) return(NULL)
  
  tryCatch({
    model <- glmmTMB(
      cbind(Minor_Count, Major_Count) ~ treatment * generation + (1 | population),
      data = snp_data,
      family = betabinomial()
    )
    
    # Extract fixed effects with chr/pos info
    broom.mixed::tidy(model, effects = "fixed") %>%
      mutate(chr = unique(snp_data$chr),
             pos = unique(snp_data$pos))
  }, error = function(e) NULL)
}

# STEP 5: Fit models across SNPs in parallel
# Only small SNP-level subsets are sent to workers
results_all <- future_map_dfr(data_split, fit_glmm, .progress = TRUE)

# STEP 6: Clean and organize result table
results_all <- results_all %>%
  select(chr, pos, term, estimate, std.error, statistic, p.value)

#save output
write.csv(results_all, "Raw_GLM_Betabi_Results.csv", row.names = FALSE)