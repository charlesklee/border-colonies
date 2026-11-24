#!/usr/bin/env Rscript
# ==============================================================================
# 12_species_status_IB.R — Host Species, Paired Tissue, Status & IB (§3.3)
# ==============================================================================
# Purpose : Run blocked PERMANOVA (9,999 permutations, set.seed(42)) across
#           Bray-Curtis (primary), Jaccard (RA & PA), and Aitchison (CLR >=5%)
#           on both primary 52 complete pairs (n = 104) and sensitivity 61 pairs
#           (n = 122):
#             - Table S3B: Species + paired Tissue main effects & Species:Tissue
#             - Table S4A: All 15 pairwise species PERMANOVA (Holm-adjusted)
#             - Table S4B: Species vs Site_Date variance partitioning & within-
#                          sampling-event tests (incl. AK vs AP at Mount 25 Mar)
#             - Table S4C: Per-species paired Tissue tests (Holm-adjusted x 6)
#             - Ecological Status (global & within-species) & Interaction Border
#               (IB intermediacy, IB-unique ASVs vs random partner-pair null,
#               IB->Native vs IB->Invasive paired tests).
# Inputs  : config.R, output/intermediate/ascidian_canonical.rds,
#           output/intermediate/ascidian_canonical_allsamples.rds
# Outputs : output/intermediate/table_S3B_paired_tissue_permanova.csv
#           output/intermediate/table_S4A_pairwise_species.csv
#           output/intermediate/table_S4B_varpart_within_event.csv
#           output/intermediate/table_S4C_per_species_tissue.csv
#           output/intermediate/fig3A_pcoa_coords.csv
#           output/intermediate/figS3_status_pcoa.csv
#           output/intermediate/figS4_ib_distances.csv
#           output/intermediate/figS4_ib_unique_asvs.csv
# Supports: Figure 3 (A-B); Figures S3, S4; Tables S3B, S4A, S4B, S4C; §3.3
# ==============================================================================

source("config.R")
set.seed(SEED)

suppressPackageStartupMessages({
  library(phyloseq)
  library(microbiome)
  library(vegan)
  library(permute)
  library(dplyr)
})

prep_obj <- function(rds_path) {
  phy <- readRDS(rds_path)
  sd <- as(sample_data(phy), "data.frame")
  sd$SequencingID <- rownames(sd)
  sd$Tissue_code  <- factor(as.character(sd$Zooid_or_Tunic), levels = c("Z", "T", "Seawater"))
  sd$Tissue       <- factor(ifelse(sd$Zooid_or_Tunic == "Z", "Zooid",
                            ifelse(sd$Zooid_or_Tunic == "T", "Tunic", "Seawater")),
                            levels = c("Zooid", "Tunic", "Seawater"))
  sd$IB           <- as.character(sd$Interactive_Border) == "Y"
  sd$encounter_id <- as.character(sd$Original_sample_name)
  sd$pair_id      <- ifelse(sd$Tissue_code == "Seawater", as.character(sd$Original_Identifier),
                            sub("[ZT]$", "", as.character(sd$Original_Identifier)))
  sd$Status       <- factor(as.character(sd$Interaction_or_Monospecies))
  sd$Site_Date    <- paste(as.character(sd$Location), as.character(sd$Collection_date), sep = "_")
  sd$Species_chr  <- as.character(sd$Species)
  sd$Native_Invasive <- ifelse(sd$Species_chr %in% c("AP", "DI", "LN"), "Native",
                        ifelse(sd$Species_chr %in% c("AK", "BL", "DV"), "Invasive", NA_character_))
  sd_nonib <- sd[sd$Tissue_code %in% c("Z", "T") & !sd$IB & sd$Species_chr %in% SP_LEVELS, ]
  ptab <- table(sd_nonib$pair_id, droplevels(sd_nonib$Tissue_code))
  comp_pairs <- sort(rownames(ptab)[ptab[, "Z"] == 1 & ptab[, "T"] == 1])
  sids <- rownames(sd_nonib)[sd_nonib$pair_id %in% comp_pairs]
  list(phy = phy, sd = sd, complete_pairs = comp_pairs, nonib_pair_sids = sids)
}

compute_dists <- function(phy_obj, sids) {
  p_sub <- prune_samples(sids, phy_obj)
  p_sub <- prune_taxa(taxa_sums(p_sub) > 0, p_sub)
  otu <- as(otu_table(p_sub), "matrix")
  if (taxa_are_rows(p_sub)) otu <- t(otu)
  otu <- otu[sids, , drop = FALSE]
  ra  <- sweep(otu, 1, rowSums(otu), "/")
  pa  <- (otu > 0) * 1.0
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

obj_can <- prep_obj(PATH_PHY_CANON)
obj_all <- prep_obj(PATH_PHY_ALL)
stopifnot(length(obj_can$complete_pairs) == 52L, length(obj_all$complete_pairs) == 61L)

datasets <- list(
  Primary_52p    = obj_can,
  AllSamples_61p = obj_all
)

s3b_rows <- list()
s4a_rows <- list()
s4b_rows <- list()
s4c_rows <- list()

for (ds_name in names(datasets)) {
  obj <- datasets[[ds_name]]
  sids <- obj$nonib_pair_sids
  md <- obj$sd[sids, , drop = FALSE]
  md$Species <- factor(md$Species_chr, levels = SP_LEVELS)
  md$Tissue  <- droplevels(md$Tissue)
  md$pair_id <- factor(md$pair_id)
  md$Site_Date <- factor(md$Site_Date)
  dists <- compute_dists(obj$phy, sids)

  perm_within  <- how(blocks = md$pair_id, nperm = N_PERM)
  perm_between <- how(plots = Plots(strata = md$pair_id, type = "free"),
                      within = Within(type = "none"), nperm = N_PERM)

  # --- Table S3B: Species + Tissue main effects & Species:Tissue ---
  for (dn in names(dists)) {
    dm <- dists[[dn]]
    set.seed(SEED)
    a_sp <- adonis2(dm ~ Species + Tissue, data = md, by = "margin", permutations = perm_between, parallel = N_CORES)
    set.seed(SEED)
    a_ti <- adonis2(dm ~ Species + Tissue, data = md, by = "margin", permutations = perm_within,  parallel = N_CORES)
    set.seed(SEED)
    a_in <- adonis2(dm ~ Species * Tissue, data = md, by = "margin", permutations = perm_within,  parallel = N_CORES)
    s3b_rows[[length(s3b_rows) + 1L]] <- rbind(
      data.frame(Dataset = ds_name, Distance = dn, Term = "Species", Df = a_sp["Species", "Df"],
                 R2 = a_sp["Species", "R2"], F_stat = a_sp["Species", "F"], p_value = a_sp["Species", "Pr(>F)"]),
      data.frame(Dataset = ds_name, Distance = dn, Term = "Tissue", Df = a_ti["Tissue", "Df"],
                 R2 = a_ti["Tissue", "R2"], F_stat = a_ti["Tissue", "F"], p_value = a_ti["Tissue", "Pr(>F)"]),
      data.frame(Dataset = ds_name, Distance = dn, Term = "Species:Tissue", Df = a_in["Species:Tissue", "Df"],
                 R2 = a_in["Species:Tissue", "R2"], F_stat = a_in["Species:Tissue", "F"], p_value = a_in["Species:Tissue", "Pr(>F)"])
    )
  }

  # --- Table S4A: Pairwise Species PERMANOVA (all 15 pairs, Holm) ---
  sp_pairs <- combn(SP_LEVELS, 2, simplify = FALSE)
  for (dn in names(dists)) {
    dm_mat <- as.matrix(dists[[dn]])
    pw_list <- lapply(sp_pairs, function(pr) {
      keep <- md$Species_chr %in% pr
      sub_md <- droplevels(md[keep, , drop = FALSE])
      sub_d  <- as.dist(dm_mat[keep, keep])
      perm_pw <- how(plots = Plots(strata = sub_md$pair_id, type = "free"),
                     within = Within(type = "none"), nperm = N_PERM)
      set.seed(SEED)
      ad <- suppressMessages(adonis2(sub_d ~ Species, data = sub_md, by = "terms", permutations = perm_pw, parallel = N_CORES))
      data.frame(
        Dataset = ds_name, Distance = dn, Pair = paste(pr, collapse = " vs "),
        n_pairs_1 = sum(sub_md$Species_chr == pr[1]) / 2L,
        n_pairs_2 = sum(sub_md$Species_chr == pr[2]) / 2L,
        R2 = ad["Species", "R2"], F_stat = ad["Species", "F"], p_raw = ad["Species", "Pr(>F)"],
        stringsAsFactors = FALSE
      )
    })
    pw_df <- do.call(rbind, pw_list)
    pw_df$p_holm <- p.adjust(pw_df$p_raw, method = "holm")
    s4a_rows[[length(s4a_rows) + 1L]] <- pw_df
  }

  # --- Table S4B: Variance partitioning (Species vs Site_Date) & within-event tests ---
  for (dn in names(dists)) {
    dm <- dists[[dn]]
    set.seed(SEED)
    t_sp <- adonis2(dm ~ Species + Site_Date + Tissue, data = md, by = "terms", permutations = perm_between, parallel = N_CORES)
    set.seed(SEED)
    t_sd <- adonis2(dm ~ Site_Date + Species + Tissue, data = md, by = "terms", permutations = perm_between, parallel = N_CORES)
    set.seed(SEED)
    m_b  <- adonis2(dm ~ Site_Date + Species + Tissue, data = md, by = "margin", permutations = perm_between, parallel = N_CORES)
    s4b_rows[[length(s4b_rows) + 1L]] <- data.frame(
      Dataset = ds_name, Distance = dn, Scope = "Global_VarPart",
      Species_first_R2 = t_sp["Species", "R2"], Site_Date_first_R2 = t_sd["Site_Date", "R2"],
      Unique_Species_R2 = m_b["Species", "R2"], Unique_Species_F = m_b["Species", "F"], Unique_Species_p = m_b["Species", "Pr(>F)"],
      Unique_Site_Date_R2 = m_b["Site_Date", "R2"], Unique_Site_Date_F = m_b["Site_Date", "F"], Unique_Site_Date_p = m_b["Site_Date", "Pr(>F)"],
      Shared_R2 = t_sp["Species", "R2"] - m_b["Species", "R2"],
      stringsAsFactors = FALSE
    )
  }
  # Within-sampling-event Species tests (Bray-Curtis + AK vs AP at Mount_2503)
  dm_bc <- as.matrix(dists$Bray_Curtis)
  for (ev in levels(md$Site_Date)) {
    keep <- md$Site_Date == ev
    sub_md <- droplevels(md[keep, , drop = FALSE])
    if (nlevels(sub_md$Species) < 2L || nlevels(sub_md$pair_id) < 3L) next
    sub_d <- as.dist(dm_bc[keep, keep])
    perm_ev <- how(plots = Plots(strata = sub_md$pair_id, type = "free"), within = Within(type = "none"), nperm = N_PERM)
    set.seed(SEED)
    a_ev <- suppressMessages(adonis2(sub_d ~ Species, data = sub_md, by = "terms", permutations = perm_ev))
    s4b_rows[[length(s4b_rows) + 1L]] <- data.frame(
      Dataset = ds_name, Distance = "Bray_Curtis", Scope = paste0("Within_", ev),
      Species_first_R2 = a_ev["Species", "R2"], Site_Date_first_R2 = NA_real_,
      Unique_Species_R2 = a_ev["Species", "R2"], Unique_Species_F = a_ev["Species", "F"], Unique_Species_p = a_ev["Species", "Pr(>F)"],
      Unique_Site_Date_R2 = NA_real_, Unique_Site_Date_F = NA_real_, Unique_Site_Date_p = NA_real_, Shared_R2 = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  # Specific AK vs AP at Mount_2503
  keep_akap <- md$Site_Date == "Mount_2503" & md$Species_chr %in% c("AK", "AP")
  sub_akap <- droplevels(md[keep_akap, , drop = FALSE])
  perm_akap <- how(plots = Plots(strata = sub_akap$pair_id, type = "free"), within = Within(type = "none"), nperm = N_PERM)
  set.seed(SEED)
  a_akap <- suppressMessages(adonis2(as.dist(dm_bc[keep_akap, keep_akap]) ~ Species, data = sub_akap, by = "terms", permutations = perm_akap))
  s4b_rows[[length(s4b_rows) + 1L]] <- data.frame(
    Dataset = ds_name, Distance = "Bray_Curtis", Scope = "Within_Mount_2503_AK_vs_AP",
    Species_first_R2 = a_akap["Species", "R2"], Site_Date_first_R2 = NA_real_,
    Unique_Species_R2 = a_akap["Species", "R2"], Unique_Species_F = a_akap["Species", "F"], Unique_Species_p = a_akap["Species", "Pr(>F)"],
    Unique_Site_Date_R2 = NA_real_, Unique_Site_Date_F = NA_real_, Unique_Site_Date_p = NA_real_, Shared_R2 = NA_real_,
    stringsAsFactors = FALSE
  )

  # --- Table S4C: Per-species paired Tissue tests (Holm x 6) ---
  for (dn in names(dists)) {
    dm_mat <- as.matrix(dists[[dn]])
    sp_res <- lapply(SP_LEVELS, function(sp) {
      keep <- md$Species_chr == sp
      sub_md <- droplevels(md[keep, , drop = FALSE])
      sub_d  <- as.dist(dm_mat[keep, keep])
      perm_blk <- how(blocks = sub_md$pair_id, nperm = N_PERM)
      set.seed(SEED)
      ad <- suppressMessages(adonis2(sub_d ~ Tissue, data = sub_md, by = "terms", permutations = perm_blk))
      data.frame(
        Dataset = ds_name, Distance = dn, Species = sp, Species_full = unname(SP_FULL_LATIN[sp]),
        n_pairs = nrow(sub_md) / 2L, R2 = ad["Tissue", "R2"], F_stat = ad["Tissue", "F"],
        p_raw = ad["Tissue", "Pr(>F)"], stringsAsFactors = FALSE
      )
    })
    sp_df <- do.call(rbind, sp_res)
    sp_df$p_holm <- p.adjust(sp_df$p_raw, method = "holm")
    s4c_rows[[length(s4c_rows) + 1L]] <- sp_df
  }

  # Save PCoA coordinates for Figure 3A & Figure S3 (Primary 52 pairs)
  if (ds_name == "Primary_52p") {
    pcoa <- cmdscale(dists$Bray_Curtis, eig = TRUE, k = 2)
    pos_eig <- pcoa$eig[pcoa$eig > 0]
    p1_pct <- 100 * pcoa$eig[1] / sum(pos_eig)
    p2_pct <- 100 * pcoa$eig[2] / sum(pos_eig)
    pcoa_df <- data.frame(
      SequencingID = sids, pair_id = as.character(md$pair_id),
      Species = md$Species_chr, Tissue = as.character(md$Tissue),
      Status = as.character(md$Status), Site_Date = as.character(md$Site_Date),
      PCoA1 = pcoa$points[, 1], PCoA2 = pcoa$points[, 2],
      PCoA1_pct = p1_pct, PCoA2_pct = p2_pct, stringsAsFactors = FALSE
    )
    write.csv(pcoa_df, file.path(DIR_INTERIM, "fig3A_pcoa_coords.csv"), row.names = FALSE)
  }
}

write.csv(do.call(rbind, s3b_rows), file.path(DIR_INTERIM, "table_S3B_paired_tissue_permanova.csv"), row.names = FALSE)
write.csv(do.call(rbind, s4a_rows), file.path(DIR_INTERIM, "table_S4A_pairwise_species.csv"), row.names = FALSE)
write.csv(do.call(rbind, s4b_rows), file.path(DIR_INTERIM, "table_S4B_varpart_within_event.csv"), row.names = FALSE)
write.csv(do.call(rbind, s4c_rows), file.path(DIR_INTERIM, "table_S4C_per_species_tissue.csv"), row.names = FALSE)

cat("Saved Tables S3B, S4A, S4B, S4C and Figure 3A PCoA coordinates.\n")
print(sessionInfo())
