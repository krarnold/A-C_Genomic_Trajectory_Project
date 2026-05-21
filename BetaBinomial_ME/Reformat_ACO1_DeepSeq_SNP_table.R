# Set data_dir to the path of your local data directory
data_dir <- "."
setwd(data_dir)

#load packages
library(dplyr)


#read in data
data <- read.csv("SNP_Table/ACO1_Deepseq_SNP_20_5000.csv", header = T)

#fix chr
data$chr[data$chr == "AE014298.5"] <- "X"
data$chr[data$chr == "AE014134.6"] <- "2L"
data$chr[data$chr == "AE013599.5"] <- "2R"
data$chr[data$chr == "AE014296.5"] <- "3L"
data$chr[data$chr == "AE014297.3"] <- "3R"
data$chr[data$chr == "AE014135.4"] <- "4"
data$chr[data$chr == "CP007106.1"] <- "Y"
data$chr[data$chr == "KJ947872.2"] <- "MT"

#subset to only have major chromsome arms and rename the ref and alt columns
data <- data[data$chr %in% c("X", "2L", "2R", "3L", "3R", "4"), ]
data <- data %>%
  rename(
    Major_Nuc = maj,
    Minor_Nuc = min
  )


# Rename minor and major count columns in data table
names(data) <- gsub("_mi$", "_minor_count", names(data))
names(data) <- gsub("_mj$", "_major_count", names(data))


#make a verison of the table that has chr, post, Pop1_MiC, Pop1_Cov, and so on. MiC is minor count.
# Step 1: Identify columns
minor_cols <- grep("_minor_count$", names(data), value = TRUE)
major_cols <- gsub("minor", "major", minor_cols)

# Step 2: Compute total coverage
coverage_cols <- gsub("_minor_count$", "_total_coverage", minor_cols)
coverage_df <- data[major_cols] + data[minor_cols]
names(coverage_df) <- coverage_cols

# Step 3: Interleave minor and coverage columns
interleaved_df <- bind_cols(data[minor_cols], coverage_df)

# Reorder: minor1, coverage1, minor2, coverage2, ...
ordered_cols <- as.vector(rbind(minor_cols, coverage_cols))

# Step 4: Final table
final_table <- bind_cols(
  data %>% select(chr, pos, Major_Nuc, Minor_Nuc),
  interleaved_df[ordered_cols]
)

write.csv(final_table, file = "ACO1_Deepseq_Modified.csv", row.names = FALSE)




