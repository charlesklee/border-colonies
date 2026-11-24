#!/usr/bin/env Rscript
# ==============================================================================
# 04_decontam_and_filters.R — Plate-Batched Decontam & Symbiont Rescue Rule
# ==============================================================================
# Purpose : Identify reagent contaminants using decontam prevalence (threshold
#           0.5, batched by library plate across 167 biological + 14 reagent
#           blanks [13 NTC + 1 extraction blank S-164]; S-163 ASW excluded from
#           negative-control pool), apply the 3-part biological rescue rule
#           (max tissue RA > 1%, frequency p > 0.1, bio167/ctrl14 ratio > 50),
#           and save the full decontam review table for Table S1B.
# Inputs  : config.R, output/intermediate/seqtab_nochim_stageB.rds,
#           output/intermediate/tax_silva138_2_minBoot50.rds,
#           data/metadata/sample_metadata.csv
# Outputs : output/intermediate/decontam_review.rds
# Supports: Methods §2.5; Table S1B (rescued ASVs); Figure S1 counts
# ==============================================================================

source("config.R")
set.seed(SEED)

suppressPackageStartupMessages({
  library(phyloseq)
  library(decontam)
  library(dplyr)
})

stopifnot(file.exists(PATH_SEQTAB_NOCHIM), file.exists(PATH_TAX_B50), file.exists(PATH_SAMPLE_META))

seqtab <- readRDS(PATH_SEQTAB_NOCHIM) # samples x seqs
tax50_obj <- readRDS(PATH_TAX_B50)
tax50 <- if (is.list(tax50_obj) && "tax" %in% names(tax50_obj)) tax50_obj$tax else tax50_obj
meta  <- read.csv(PATH_SAMPLE_META, stringsAsFactors = FALSE)
rownames(meta) <- meta$SequencingID

# Align samples
common_samps <- intersect(rownames(seqtab), rownames(meta))
stopifnot(length(common_samps) == 182L)
seqtab <- seqtab[common_samps, , drop = FALSE]
meta   <- meta[common_samps, , drop = FALSE]
meta$is.neg <- meta$Sample_Category %in% c("PCR_NTC", "Extraction_blank")

otu_all <- t(seqtab) # seqs x samples
bio167_ids <- meta$SequencingID[meta$Sample_Category %in% c("Ascidian_tissue", "Seawater")]
tis164_ids <- meta$SequencingID[meta$Sample_Category == "Ascidian_tissue"]
ctrl14_ids <- meta$SequencingID[meta$is.neg]
asw_id     <- "S-163"
ids_181    <- c(bio167_ids, ctrl14_ids)

# 1. Plate-batched prevalence decontam (threshold = 0.5, 181 libraries)
phy_181 <- phyloseq(
  otu_table(otu_all[, ids_181, drop = FALSE], taxa_are_rows = TRUE),
  sample_data(meta[ids_181, , drop = FALSE])
)
phy_181 <- prune_taxa(taxa_sums(phy_181) > 0, phy_181)
contam_prev <- isContaminant(
  phy_181, method = "prevalence", neg = "is.neg",
  batch = "Plate", threshold = DECONTAM_PREV_THRESH
)

# 2. Frequency decontam across 167 biological samples using DNA concentration
phy_167 <- phyloseq(
  otu_table(otu_all[, bio167_ids, drop = FALSE], taxa_are_rows = TRUE),
  sample_data(meta[bio167_ids, , drop = FALSE])
)
phy_167 <- prune_taxa(taxa_sums(phy_167) > 0, phy_167)
contam_freq <- isContaminant(
  phy_167, method = "frequency", conc = "DNA_ng_uL", threshold = DECONTAM_PREV_THRESH
)

# 3. Assemble per-ASV metrics & apply 3-part rescue rule
seqs_181 <- rownames(contam_prev)
ra_tis164 <- sweep(otu_all[seqs_181, tis164_ids, drop = FALSE], 2,
                   colSums(otu_all[, tis164_ids, drop = FALSE]), "/")
max_tis_ra <- apply(ra_tis164, 1, max)
reads_bio167 <- rowSums(otu_all[seqs_181, bio167_ids, drop = FALSE])
reads_ctrl14 <- rowSums(otu_all[seqs_181, ctrl14_ids, drop = FALSE])
ratio_bio_ctrl14 <- reads_bio167 / pmax(reads_ctrl14, 0)
freq_p <- contam_freq$p[match(seqs_181, rownames(contam_freq))]

prev_flag <- !is.na(contam_prev$contaminant) & contam_prev$contaminant
rescued <- prev_flag &
           (max_tis_ra > RESCUE_MIN_RA) &
           (!is.na(freq_p) & freq_p > RESCUE_FREQ_P) &
           (ratio_bio_ctrl14 > RESCUE_RATIO)
final_contam <- prev_flag & !rescued

review_df <- data.frame(
  Sequence          = seqs_181,
  Length_bp         = nchar(seqs_181),
  Reads_Bio167      = as.integer(reads_bio167),
  Reads_Ctrl14      = as.integer(reads_ctrl14),
  Reads_S163        = as.integer(otu_all[seqs_181, asw_id]),
  Max_Tissue_RA_Pct = 100 * max_tis_ra,
  Ratio_Bio167_Ctrl14 = ratio_bio_ctrl14,
  Prev_Plate_p      = contam_prev$p,
  Prev_Plate_Flag   = prev_flag,
  Freq_Bio_p        = freq_p,
  Rescued_By_Rule   = rescued,
  Final_Contaminant = final_contam,
  stringsAsFactors  = FALSE
)
tx_sub <- as.data.frame(tax50[seqs_181, , drop = FALSE], stringsAsFactors = FALSE)
review_df <- cbind(review_df, tx_sub)

cat(sprintf("Prev-plate flagged : %d ASVs\n", sum(prev_flag)))
cat(sprintf("Rescued by rule    : %d ASVs\n", sum(rescued)))
cat(sprintf("Final contaminants : %d ASVs removed\n", sum(final_contam)))

saveRDS(list(review_df = review_df, contam_prev = contam_prev, contam_freq = contam_freq),
        PATH_DECONTAM_REV)
print(sessionInfo())
