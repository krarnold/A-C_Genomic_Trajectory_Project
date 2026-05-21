library(dplyr)
library(GenomicRanges)
library(VariantAnnotation)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)

# ── PATHS ─────────────────────────────────────────────────────────────────────
# Set data_dir to the path of your local data directory
data_dir <- "."
# ──────────────────────────────────────────────────────────────────────────────

# STEP 1: Load significant interaction SNPs and build lead SNP table
sig <- read.csv(file.path(data_dir, "Significant_SNPs_treatmentC_A_generation.csv"))

lead_snps_full <- sig %>%
  arrange(chr, pos) %>%
  mutate(window = floor(pos / 10000)) %>%
  group_by(chr, window) %>%
  slice_min(p.value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(chr, pos, p.value, fdr, effect_size)

write.csv(
  lead_snps_full,
  file.path(data_dir, "Lead_SNPs_Interaction.csv"),
  row.names = FALSE
)

message("Total lead SNPs: ", nrow(lead_snps_full))

# STEP 2: Top 100 by p-value
lead_snps <- lead_snps_full %>%
  arrange(fdr, p.value) %>%
  head(100)

# STEP 3: GRanges object
gr <- GRanges(
  seqnames = paste0("chr", lead_snps$chr),
  ranges   = IRanges(start = lead_snps$pos, end = lead_snps$pos),
  pval     = lead_snps$p.value,
  fdr      = lead_snps$fdr,
  eff_size = lead_snps$effect_size
)

# STEP 4: Locate variants
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
loc  <- locateVariants(gr, txdb, AllVariants())

# STEP 5: TxDb returns FlyBase gene IDs (FBgn...) — map to gene symbols
gene_ids <- as.character(na.omit(unique(loc$GENEID)))

gene_map <- as.data.frame(
  AnnotationDbi::select(
    x        = org.Dm.eg.db,
    keys     = gene_ids,
    columns  = c("SYMBOL", "GENENAME"),
    keytype  = "FLYBASE"
  )
)

# STEP 6: Per-SNP annotation summary — join on FLYBASE
loc_df <- as.data.frame(loc) %>%
  mutate(GENEID = as.character(GENEID)) %>%
  left_join(gene_map, by = c("GENEID" = "FLYBASE")) %>%
  group_by(QUERYID) %>%
  summarise(
    location = paste(unique(LOCATION),          collapse = ";"),
    symbol   = paste(unique(na.omit(SYMBOL)),   collapse = ";"),
    flybase  = paste(unique(na.omit(GENEID)),   collapse = ";"),
    genename = paste(unique(na.omit(GENENAME)), collapse = ";"),
    .groups  = "drop"
  )

annotated <- lead_snps %>%
  mutate(QUERYID = row_number()) %>%
  left_join(loc_df, by = "QUERYID") %>%
  dplyr::select(chr, pos, p.value, fdr, effect_size, location, symbol, flybase, genename)

# STEP 7: Save
write.csv(
  annotated,
  file.path(data_dir, "Top100_Lead_SNPs_Annotated.csv"),
  row.names = FALSE
)

head(annotated, 20)
