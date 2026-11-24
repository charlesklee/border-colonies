#!/usr/bin/env Rscript
# ==============================================================================
# 05_build_canonical_object.R — Organelle Filter, Depth Filter & ASV_n Numbering
# ==============================================================================
# Purpose : Remove final contaminants (92 ASVs), remove non-bacterial/organelle
#           ASVs (Kingdom == "Eukaryota" | is.na(Kingdom), Order == "Chloroplast",
#           Family == "Mitochondria"), assign canonical ASV_1..ASV_n names
#           matching data/metadata/asv_crosswalk.csv, and write both
#           ascidian_canonical_allsamples.rds (167 samples, 18,174 ASVs) and
#           ascidian_canonical.rds (>=1,000 reads: 151 samples, 18,131 ASVs).
# Inputs  : config.R, output/intermediate/seqtab_nochim_stageB.rds,
#           output/intermediate/tax_silva138_2_minBoot50.rds,
#           output/intermediate/decontam_review.rds,
#           data/metadata/sample_metadata.csv, data/metadata/asv_crosswalk.csv
# Outputs : output/intermediate/ascidian_canonical_allsamples.rds
#           output/intermediate/ascidian_canonical.rds
# Supports: Methods §2.5; Table S1A–B; Figure S1; all downstream analyses
# ==============================================================================

source("config.R")
set.seed(SEED)

suppressPackageStartupMessages({
  library(phyloseq)
})

seqtab    <- readRDS(PATH_SEQTAB_NOCHIM)
tax50_obj <- readRDS(PATH_TAX_B50)
tax50     <- if (is.list(tax50_obj) && "tax" %in% names(tax50_obj)) tax50_obj$tax else tax50_obj
drev      <- readRDS(PATH_DECONTAM_REV)$review_df
meta      <- read.csv(PATH_SAMPLE_META, stringsAsFactors = FALSE)
rownames(meta) <- meta$SequencingID

bio167_ids <- meta$SequencingID[meta$Sample_Category %in% c("Ascidian_tissue", "Seawater")]
otu_bio <- t(seqtab[bio167_ids, , drop = FALSE]) # seqs x 167 bio samples
otu_bio <- otu_bio[rowSums(otu_bio) > 0, , drop = FALSE]

# 1. Remove final contaminants (92 ASVs total, 87 present in bio167)
rem_seqs <- drev$Sequence[drev$Final_Contaminant]
otu_clean <- otu_bio[!(rownames(otu_bio) %in% rem_seqs), , drop = FALSE]

# 2. Organelle & non-prokaryote filter (minBoot = 50 canonical taxonomy)
tx <- tax50[rownames(otu_clean), , drop = FALSE]
keep_org <- !is.na(tx[, "Kingdom"]) &
            tx[, "Kingdom"] != "Eukaryota" &
            (is.na(tx[, "Order"])  | tx[, "Order"]  != "Chloroplast") &
            (is.na(tx[, "Family"]) | tx[, "Family"] != "Mitochondria")

# Report rescued ASVs surviving organelle filter (expected: 20 / 28)
resc_seqs <- drev$Sequence[drev$Rescued_By_Rule]
n_resc_post_org <- sum(rownames(otu_clean)[keep_org] %in% resc_seqs)
cat(sprintf("Rescued ASVs surviving organelle filter: %d / %d\n", n_resc_post_org, length(resc_seqs)))

otu_clean <- otu_clean[keep_org, , drop = FALSE]
tx_clean  <- tx[keep_org, , drop = FALSE]

# 2b. Abundance prevalence filter: retain ASVs with > 2 reads in >= 1 biological sample
keep_gt2  <- rowSums(otu_clean > MIN_ASV_SAMPLE_GT) >= 1L
otu_clean <- otu_clean[keep_gt2, , drop = FALSE]
tx_clean  <- tx_clean[keep_gt2, , drop = FALSE]

# 3. Assign canonical ASV_n IDs matching data/metadata/asv_crosswalk.csv (order by bio167 reads)
if (file.exists(PATH_ASV_CROSSWALK)) {
  cw <- read.csv(PATH_ASV_CROSSWALK, stringsAsFactors = FALSE)
  m_cw <- match(rownames(otu_clean), cw$ASV_DNA)
  if (!any(is.na(m_cw))) {
    asv_names <- cw$ASV_ID[m_cw]
  } else {
    ord_bio <- order(-rowSums(otu_clean))
    otu_clean <- otu_clean[ord_bio, , drop = FALSE]
    tx_clean  <- tx_clean[ord_bio, , drop = FALSE]
    asv_names <- paste0("ASV_", seq_len(nrow(otu_clean)))
  }
} else {
  ord_bio <- order(-rowSums(otu_clean))
  otu_clean <- otu_clean[ord_bio, , drop = FALSE]
  tx_clean  <- tx_clean[ord_bio, , drop = FALSE]
  asv_names <- paste0("ASV_", seq_len(nrow(otu_clean)))
}
rownames(otu_clean) <- asv_names
rownames(tx_clean)  <- asv_names

# Sort by ASV number
ord_num <- order(as.integer(sub("^ASV_", "", rownames(otu_clean))))
otu_clean <- otu_clean[ord_num, , drop = FALSE]
tx_clean  <- tx_clean[ord_num, , drop = FALSE]

# 4. Build AllSamples canonical phyloseq (167 samples)
phy_all <- phyloseq(
  otu_table(otu_clean, taxa_are_rows = TRUE),
  tax_table(as.matrix(tx_clean)),
  sample_data(meta[bio167_ids, , drop = FALSE])
)
saveRDS(phy_all, PATH_PHY_ALL)

# 5. Depth filter (>= 1,000 reads -> 151 samples; 16 excluded)
keep_samps <- sample_sums(phy_all) >= MIN_READS_SAMPLE
cat(sprintf("Samples excluded (<%d reads): %d\n", MIN_READS_SAMPLE, sum(!keep_samps)))
phy_can <- prune_samples(keep_samps, phy_all)
phy_can <- prune_taxa(taxa_sums(phy_can) > 0, phy_can)
saveRDS(phy_can, PATH_PHY_CANON)

cat(sprintf("Saved %s: %d samples, %d ASVs (%s reads)\n",
            PATH_PHY_ALL, nsamples(phy_all), ntaxa(phy_all),
            format(sum(sample_sums(phy_all)), big.mark = ",")))
cat(sprintf("Saved %s: %d samples, %d ASVs (%s reads)\n",
            PATH_PHY_CANON, nsamples(phy_can), ntaxa(phy_can),
            format(sum(sample_sums(phy_can)), big.mark = ",")))
print(sessionInfo())
