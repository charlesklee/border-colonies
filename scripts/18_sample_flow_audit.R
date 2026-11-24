#!/usr/bin/env Rscript
# ==============================================================================
# 18_sample_flow_audit.R — Sample-Flow Reconciliation & Table S1A (Excluded)
# ==============================================================================
# Purpose : Reconcile all 40 collected ascidian specimens (19 monospecific + 21
#           interaction across 7 pairings), 164 dissected tissues, 3 seawater,
#           15 controls (182 ENA libraries), 16 excluded (<1,000 reads), and
#           151 retained libraries (52 complete non-IB pairs + 6 unpaired +
#           38 IB [33 complete triplets] + 3 seawater). Exports Table S1A.
# Inputs  : config.R, data/metadata/sample_metadata.csv,
#           output/intermediate/ascidian_canonical_allsamples.rds,
#           output/intermediate/ascidian_canonical.rds
# Outputs : output/intermediate/table_S1A_excluded_samples.csv
#           output/intermediate/sample_flow_summary.csv
# Supports: Figure S1 counts (docs/figS1_prompt.md); Table S1A; §2.1–§2.5
# ==============================================================================

source("config.R")
set.seed(SEED)

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
})

meta    <- read.csv(PATH_SAMPLE_META, stringsAsFactors = FALSE)
phy_all <- readRDS(PATH_PHY_ALL)
phy_can <- readRDS(PATH_PHY_CANON)

sd_all <- as(sample_data(phy_all), "data.frame")
sd_all$SequencingID <- rownames(sd_all)
sd_all$Reads <- as.integer(sample_sums(phy_all))
sd_all$Retained <- sd_all$SequencingID %in% sample_names(phy_can)

# Table S1A: 16 excluded samples (<1,000 reads post-decontam)
excl <- sd_all[!sd_all$Retained, c("SequencingID", "Original_Identifier", "Original_sample_name",
                                   "Species", "Zooid_or_Tunic", "Interactive_Border",
                                   "Interaction_or_Monospecies", "Location", "Collection_date", "Reads")]
excl <- excl[order(excl$Interactive_Border, excl$Species, excl$Zooid_or_Tunic), ]
write.csv(excl, file.path(DIR_INTERIM, "table_S1A_excluded_samples.csv"), row.names = FALSE)

# Summary counts for Figure S1
ret <- sd_all[sd_all$Retained, ]
ret$IB <- as.character(ret$Interactive_Border) == "Y"
ret$pair_id <- ifelse(ret$Zooid_or_Tunic == "Seawater", as.character(ret$Original_Identifier),
                      sub("[ZT]$", "", as.character(ret$Original_Identifier)))
ret_nonib <- ret[ret$Zooid_or_Tunic %in% c("Z", "T") & !ret$IB, ]
ptab <- table(ret_nonib$pair_id, ret_nonib$Zooid_or_Tunic)
n_pairs <- sum(ptab[, "Z"] == 1 & ptab[, "T"] == 1)

summ <- data.frame(
  Metric = c("Total_libraries_ENA", "Biological_libraries", "Controls",
             "Excluded_lt1000_reads", "Retained_canonical",
             "Complete_nonIB_pairs", "Unpaired_nonIB_retained",
             "Retained_IB_samples", "Retained_seawater"),
  Count = c(nrow(meta), nsamples(phy_all), nrow(meta) - nsamples(phy_all),
            nrow(excl), nsamples(phy_can),
            n_pairs, nrow(ret_nonib) - 2L * n_pairs,
            sum(ret$IB), sum(ret$Zooid_or_Tunic == "Seawater")),
  stringsAsFactors = FALSE
)
write.csv(summ, file.path(DIR_INTERIM, "sample_flow_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)
print(sessionInfo())
