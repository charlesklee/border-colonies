#!/usr/bin/env Rscript
# ==============================================================================
# 03_taxonomy.R — SILVA v138.2 Taxonomic Assignment (minBoot 50 & 30)
# ==============================================================================
# Purpose : Assign taxonomy against SILVA v138.2 at canonical minBoot = 50 and
#           secondary lower-threshold minBoot = 30 (used solely for the
#           "(bootstrap <50)" best_assignment_minBoot30 labels in Tables S2 & S7A).
# Inputs  : config.R, output/intermediate/seqtab_nochim_stageB.rds,
#           data/reference/silva_nr99_v138.2_toSpecies_trainset.fa.gz
# Outputs : output/intermediate/tax_silva138_2_minBoot50.rds (canonical)
#           output/intermediate/tax_silva138_2_minBoot30.rds (secondary labels)
# Supports: Methods §2.4; Tables S2, S7A; Figure 4A
# ==============================================================================

source("config.R")
set.seed(SEED)

suppressPackageStartupMessages({
  library(dada2)
})

stopifnot(file.exists(PATH_SEQTAB_NOCHIM), file.exists(PATH_SILVA_FASTA))
seqtab_nochim <- readRDS(PATH_SEQTAB_NOCHIM)
seqs <- colnames(seqtab_nochim)
cat(sprintf("=== 03_taxonomy.R: Assigning SILVA v138.2 taxonomy to %d ASVs ===\n", length(seqs)))

# Canonical assignment: minBoot = 50
set.seed(SEED)
tax_b50 <- assignTaxonomy(
  seqs, PATH_SILVA_FASTA,
  minBoot = MIN_BOOT_CANON,
  tryRC = FALSE,
  outputBootstraps = TRUE,
  multithread = N_CORES
)
saveRDS(tax_b50, PATH_TAX_B50)

# Secondary assignment: minBoot = 30 (for best_assignment_minBoot30 column in Tables S2 & S7A)
set.seed(SEED)
tax_b30 <- assignTaxonomy(
  seqs, PATH_SILVA_FASTA,
  minBoot = MIN_BOOT_LOW,
  tryRC = FALSE,
  outputBootstraps = TRUE,
  multithread = N_CORES
)
saveRDS(tax_b30, PATH_TAX_B30)

cat(sprintf("Saved %s (minBoot=%d) and %s (minBoot=%d).\n",
            PATH_TAX_B50, MIN_BOOT_CANON, PATH_TAX_B30, MIN_BOOT_LOW))
print(sessionInfo())
