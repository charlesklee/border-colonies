#!/usr/bin/env Rscript
# ==============================================================================
# 21_supplementary_tables.R — Consolidate All Manuscript Tables (Tables S1–S8)
# ==============================================================================
# Purpose : Programmatically export every supplementary table reported in the
#           manuscript into output/tables/ with standardized filenames:
#             - Table S1A: Excluded samples (<1,000 reads) [from 18_sample_flow_audit.R]
#             - Table S1B: 28 rescued ASVs (with post-organelle flag) [from 04_decontam_and_filters.R]
#             - Table S1C: Control inventory & attribution
#             - Table S1D: Seawater ASV_1:2:3 ratios [from 15_seawater_and_carryover.R]
#             - Table S2 : BLAST full table & summary [from 16_blast.py / 16_blast_local16s.R]
#             - Table S3A–C: Sample-type PERMANOVA, paired tissue tests, PERMDISP
#             - Table S4A–C: Pairwise species, variance partitioning/within-event, per-species tissue
#             - Table S5A–B: Bacterioplankton & Synechococcus enrichment ratios
#             - Table S6A–C: Category RA by species, Alphaproteobacteria %, top-5 classes
#             - Table S7A–C: Dominant ASVs, pooled prevalence, within-specimen co-occurrence
#             - Table S8   : Differential abundance across 3 frameworks
# Inputs  : config.R, output/intermediate/*.rds & *.csv
# Outputs : output/tables/TableS1A_*.csv .. TableS8_*.csv
# Supports: Tables S1–S8
# ==============================================================================

source("config.R")
set.seed(SEED)

# Programmatically assemble Table S1B (28 rescued ASVs + post-organelle status)
if (file.exists(PATH_DECONTAM_REV)) {
  drev <- readRDS(PATH_DECONTAM_REV)$review_df
  resc <- drev[drev$Rescued_By_Rule, , drop = FALSE]
  resc$Retained_Post_Organelle <- !is.na(resc$Kingdom) & resc$Kingdom != "Eukaryota" &
                                  (is.na(resc$Order) | resc$Order != "Chloroplast") &
                                  (is.na(resc$Family) | resc$Family != "Mitochondria")
  if (file.exists(PATH_ASV_CROSSWALK)) {
    cw <- read.csv(PATH_ASV_CROSSWALK, stringsAsFactors = FALSE)
    resc$ASV_ID <- cw$ASV_ID[match(resc$Sequence, cw$ASV_DNA)]
  }
  resc <- resc[order(-resc$Reads_Bio167), ]
  write.csv(resc, file.path(DIR_TABLES, "TableS1B_rescued_ASVs.csv"), row.names = FALSE)
}

# Copy all generated intermediate supplementary tables to standardized TableS# names
copy_map <- c(
  "table_S1A_excluded_samples.csv"           = "TableS1A_excluded_samples.csv",
  "table_S1C_control_attribution.csv"        = "TableS1C_control_attribution.csv",
  "table_S1D_seawater_ratios.csv"            = "TableS1D_seawater_ASV123_ratios.csv",
  "table_S3A_sample_type_permanova.csv"      = "TableS3A_sample_type_permanova.csv",
  "table_S3B_paired_tissue_permanova.csv"    = "TableS3B_paired_tissue_permanova.csv",
  "table_S3C_permdisp.csv"                   = "TableS3C_permdisp.csv",
  "table_S4A_pairwise_species.csv"           = "TableS4A_pairwise_species_permanova.csv",
  "table_S4B_varpart_within_event.csv"       = "TableS4B_variance_partitioning_within_event.csv",
  "table_S4C_per_species_tissue.csv"         = "TableS4C_per_species_tissue_permanova.csv",
  "table_S5A_bacterioplankton_enrichment.csv"= "TableS5A_bacterioplankton_enrichment.csv",
  "table_S5B_synechococcus_by_species.csv"   = "TableS5B_synechococcus_by_species.csv",
  "table_S6A_categories_by_species.csv"      = "TableS6A_categories_by_species.csv",
  "table_S6B_category_alpha.csv"             = "TableS6B_category_alphaproteobacteria.csv",
  "table_S6C_category_classes.csv"           = "TableS6C_category_classes.csv",
  "table_S7A_dominant.csv"                   = "TableS7A_dominant_by_species_tissue.csv",
  "table_S7B_prevalence.csv"                 = "TableS7B_pooled_prevalence.csv",
  "table_S7B_prevalence_summary.csv"         = "TableS7B_prevalence_thresholds.csv",
  "table_S7C_cooccurrence.csv"               = "TableS7C_within_specimen_cooccurrence.csv",
  "table_S8_differential_abundance.csv"      = "TableS8_differential_abundance.csv"
)

for (src in names(copy_map)) {
  sp <- file.path(DIR_INTERIM, src)
  if (file.exists(sp)) {
    file.copy(sp, file.path(DIR_TABLES, copy_map[[src]]), overwrite = TRUE)
    cat(sprintf("  Exported %s\n", copy_map[[src]]))
  }
}

cat("Supplementary table consolidation complete (output/tables/).\n")
print(sessionInfo())
