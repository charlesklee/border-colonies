#!/usr/bin/env Rscript
# ==============================================================================
# 10_sample_type_permanova.R — Sample-Type PERMANOVA & PERMDISP (Table S3A, S3C)
# ==============================================================================
# Purpose : Test global and pairwise community differences among Zooid (Z),
#           Tunic (T), and Seawater (SW) across three sample subsets:
#             (1) Primary 52 complete non-IB pairs + 3 SW (n = 107)
#             (2) Canonical IB-included (n = 151)
#             (3) All-samples sensitivity (n = 167)
#           Computes Bray-Curtis (primary), Jaccard (RA & PA), and Aitchison
#           (CLR on >=5% prevalence ASVs), plus multivariate dispersion
#           (vegan::betadisper + permutest, 9,999 permutations).
# Inputs  : config.R, output/intermediate/ascidian_canonical.rds,
#           output/intermediate/ascidian_canonical_allsamples.rds
# Outputs : output/intermediate/table_S3A_sample_type_permanova.csv
#           output/intermediate/table_S3C_permdisp.csv
# Supports: §3.1; Table S3A, Table S3C
# ==============================================================================

source("config.R")
set.seed(SEED)

suppressPackageStartupMessages({
  library(phyloseq)
  library(microbiome)
  library(vegan)
  library(dplyr)
})

phy_can <- readRDS(PATH_PHY_CANON)
phy_all <- readRDS(PATH_PHY_ALL)

prep_meta <- function(phy) {
  sd <- as(sample_data(phy), "data.frame")
  sd$SequencingID <- rownames(sd)
  sd$Tissue <- factor(ifelse(sd$Zooid_or_Tunic == "Z", "Zooid",
                      ifelse(sd$Zooid_or_Tunic == "T", "Tunic", "Seawater")),
                      levels = c("Zooid", "Tunic", "Seawater"))
  sd$IB <- as.character(sd$Interactive_Border) == "Y"
  sd$pair_id <- ifelse(sd$Zooid_or_Tunic == "Seawater",
                       as.character(sd$Original_Identifier),
                       sub("[ZT]$", "", as.character(sd$Original_Identifier)))
  sd
}

sd_can <- prep_meta(phy_can)
sd_all <- prep_meta(phy_all)

# 52 complete non-IB pairs + 3 SW (n = 107)
sd_nonib <- sd_can[sd_can$Zooid_or_Tunic %in% c("Z", "T") & !sd_can$IB & sd_can$Species %in% SP_LEVELS, ]
ptab <- table(sd_nonib$pair_id, sd_nonib$Zooid_or_Tunic)
comp52 <- sort(rownames(ptab)[ptab[, "Z"] == 1 & ptab[, "T"] == 1])
ids_107 <- c(sd_nonib$SequencingID[sd_nonib$pair_id %in% comp52],
             sd_can$SequencingID[sd_can$Zooid_or_Tunic == "Seawater"])
stopifnot(length(ids_107) == 107L)

compute_dists <- function(phy_sub, sample_order) {
  phy_sub <- prune_taxa(taxa_sums(phy_sub) > 0, phy_sub)
  otu <- as(otu_table(phy_sub), "matrix")
  if (taxa_are_rows(phy_sub)) otu <- t(otu) # samples x taxa
  otu <- otu[sample_order, , drop = FALSE]
  ra <- sweep(otu, 1, rowSums(otu), "/")
  pa <- (otu > 0) * 1.0
  min_s <- ceiling(AITCHISON_MIN_PREV * nrow(otu))
  otu_p5 <- otu[, colSums(otu > 0) >= min_s, drop = FALSE]
  phy_p5 <- phyloseq(otu_table(t(otu_p5), taxa_are_rows = TRUE))
  clr_mat <- t(as(otu_table(microbiome::transform(phy_p5, "clr")), "matrix"))
  list(
    Bray_Curtis = vegdist(ra, method = "bray"),
    Jaccard_RA  = vegdist(ra, method = "jaccard", binary = FALSE),
    Jaccard_PA  = vegdist(pa, method = "jaccard", binary = TRUE),
    Aitchison   = dist(clr_mat, method = "euclidean")
  )
}

subsets <- list(
  Primary_107    = list(phy = prune_samples(ids_107, phy_can), md = sd_can[ids_107, ]),
  Canonical_151  = list(phy = phy_can, md = sd_can),
  AllSamples_167 = list(phy = phy_all, md = sd_all)
)

perm_rows <- list()
disp_rows <- list()

for (sname in names(subsets)) {
  cat(sprintf("\n--- Running %s (n = %d) ---\n", sname, nrow(subsets[[sname]]$md)))
  md <- subsets[[sname]]$md
  dists <- compute_dists(subsets[[sname]]$phy, rownames(md))
  for (dname in names(dists)) {
    dmat <- dists[[dname]]
    set.seed(SEED)
    ag <- adonis2(dmat ~ Tissue, data = md, permutations = N_PERM, parallel = N_CORES)
    perm_rows[[length(perm_rows) + 1L]] <- data.frame(
      Subset = sname, n_samples = nrow(md), Distance = dname, Comparison = "Global (Z vs T vs SW)",
      Df = ag$Df[1], R2 = ag$R2[1], F_stat = ag$F[1], p_raw = ag$`Pr(>F)`[1], p_bonf = ag$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    )
    # Pairwise Z-T, Z-SW, T-SW (Bonferroni x 3)
    pw_list <- list(c("Zooid", "Tunic"), c("Zooid", "Seawater"), c("Tunic", "Seawater"))
    for (pr in pw_list) {
      keep <- md$Tissue %in% pr
      sub_md <- droplevels(md[keep, ])
      sub_d  <- as.dist(as.matrix(dmat)[keep, keep])
      set.seed(SEED)
      ap <- adonis2(sub_d ~ Tissue, data = sub_md, permutations = N_PERM, parallel = N_CORES)
      perm_rows[[length(perm_rows) + 1L]] <- data.frame(
        Subset = sname, n_samples = nrow(sub_md), Distance = dname,
        Comparison = paste(pr, collapse = " vs "),
        Df = ap$Df[1], R2 = ap$R2[1], F_stat = ap$F[1],
        p_raw = ap$`Pr(>F)`[1], p_bonf = min(1, 3 * ap$`Pr(>F)`[1]),
        stringsAsFactors = FALSE
      )
    }
    # PERMDISP (betadisper + permutest)
    bd <- betadisper(dmat, group = md$Tissue)
    set.seed(SEED)
    pt <- permutest(bd, pairwise = TRUE, permutations = N_PERM)
    grp_means <- tapply(bd$distances, md$Tissue, mean)
    disp_rows[[length(disp_rows) + 1L]] <- data.frame(
      Subset = sname, n_samples = nrow(md), Distance = dname,
      Mean_dist_Zooid = unname(grp_means["Zooid"]),
      Mean_dist_Tunic = unname(grp_means["Tunic"]),
      Mean_dist_Seawater = unname(grp_means["Seawater"]),
      F_stat = pt$tab$F[1], p_global = pt$tab$`Pr(>F)`[1],
      p_Z_vs_T = unname(pt$pairwise$permuted["Zooid-Tunic"]),
      p_Z_vs_SW = unname(pt$pairwise$permuted["Zooid-Seawater"]),
      p_T_vs_SW = unname(pt$pairwise$permuted["Tunic-Seawater"]),
      stringsAsFactors = FALSE
    )
  }
}

df_perm <- do.call(rbind, perm_rows)
df_disp <- do.call(rbind, disp_rows)
write.csv(df_perm, file.path(DIR_INTERIM, "table_S3A_sample_type_permanova.csv"), row.names = FALSE)
write.csv(df_disp, file.path(DIR_INTERIM, "table_S3C_permdisp.csv"), row.names = FALSE)
cat("Saved Table S3A and Table S3C intermediate CSVs.\n")
print(sessionInfo())
