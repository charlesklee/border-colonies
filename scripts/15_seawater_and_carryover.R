#!/usr/bin/env Rscript
# ==============================================================================
# 15_seawater_and_carryover.R — Control Attribution, Seawater Ratios & S5A–B
# ==============================================================================
# Purpose : Compute (1) NNLS reagent-vs-carryover control attribution (Table S1C);
#           (2) seawater ASV_1:ASV_2:ASV_3 ratio & 8-neighbour carry-over check
#           (Table S1D); (3) species-aware seawater-enriched vs host-enriched
#           classification; and (4) bacterioplankton clade & Synechococcus CC9902
#           enrichment ratios by tissue and host species (Table S5A, Table S5B).
#           (Excludes superseded Block 14b directional split, index-misassignment
#           check, and Z-T barcode check per manuscript scope.)
# Inputs  : config.R, output/intermediate/ascidian_canonical.rds,
#           output/intermediate/seqtab_nochim_stageB.rds,
#           data/metadata/sample_metadata.csv
# Outputs : output/intermediate/table_S1C_control_attribution.csv
#           output/intermediate/table_S1D_seawater_ratios.csv
#           output/intermediate/table_S5A_bacterioplankton_enrichment.csv
#           output/intermediate/table_S5B_synechococcus_by_species.csv
# Supports: §2.4, §3.2; Tables S1C, S1D, S5A, S5B
# ==============================================================================

source("config.R")
set.seed(SEED)

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(quadprog)
})

phy_can <- readRDS(PATH_PHY_CANON)
sd <- as(sample_data(phy_can), "data.frame")
sd$SequencingID <- rownames(sd)
sd$IB <- as.character(sd$Interactive_Border) == "Y"
sd$pair_id <- ifelse(sd$Zooid_or_Tunic == "Seawater", as.character(sd$Original_Identifier),
                     sub("[ZT]$", "", as.character(sd$Original_Identifier)))

sd_nonib <- sd[sd$Zooid_or_Tunic %in% c("Z", "T") & !sd$IB & sd$Species %in% SP_LEVELS, ]
ptab <- table(sd_nonib$pair_id, sd_nonib$Zooid_or_Tunic)
comp52 <- sort(rownames(ptab)[ptab[, "Z"] == 1 & ptab[, "T"] == 1])
z_ids  <- sapply(comp52, function(p) sd_nonib$SequencingID[sd_nonib$pair_id == p & sd_nonib$Zooid_or_Tunic == "Z"])
t_ids  <- sapply(comp52, function(p) sd_nonib$SequencingID[sd_nonib$pair_id == p & sd_nonib$Zooid_or_Tunic == "T"])
sw_ids <- sort(sd$SequencingID[sd$Zooid_or_Tunic == "Seawater"])
sids_107 <- c(z_ids, t_ids, sw_ids)

phy_107 <- prune_samples(sids_107, phy_can)
phy_107 <- prune_taxa(taxa_sums(phy_107) > 0, phy_107)
otu <- as(otu_table(phy_107), "matrix")
if (!taxa_are_rows(phy_107)) otu <- t(otu)
otu <- otu[, sids_107, drop = FALSE]
ra  <- sweep(otu, 2, colSums(otu), "/")
tax <- as.data.frame(tax_table(phy_107), stringsAsFactors = FALSE)

# 0. Table S1C: Control inventory, reagent-vs-tissue partition & NNLS attribution
#    Matches verify_block13_part2.R:
#      grp_reagent   = 92 Final_Contaminant | 11 non-flagged S-163+Ctrl14-only (103 seqs)
#      grp_ctrl_only = non-FinalContam Ctrl14-only (0 bio reads, 0 S-163 reads)
#      grp_s163_only = S-163-only (0 bio reads, 0 Ctrl14 reads)
#      grp_clean_bio = non-FinalContam seqs with Reads_Bio167 > 0 (tissue-derived)
if (file.exists(PATH_SEQTAB_NOCHIM) && file.exists(PATH_DECONTAM_REV)) {
  seqtab_all <- readRDS(PATH_SEQTAB_NOCHIM)
  drev_all   <- readRDS(PATH_DECONTAM_REV)$review_df
  meta_all   <- read.csv(PATH_SAMPLE_META, stringsAsFactors = FALSE)
  rownames(meta_all) <- meta_all$SequencingID

  seqs_all   <- colnames(seqtab_all)
  otu_all    <- t(seqtab_all[rownames(meta_all), seqs_all, drop = FALSE]) # all 22,409 seqs x 182
  is_final_contam <- seqs_all %in% drev_all$Sequence[drev_all$Final_Contaminant]
  bio167_ids <- meta_all$SequencingID[meta_all$Sample_Category %in% c("Ascidian_tissue", "Seawater")]
  ctrl15_ids <- meta_all$SequencingID[meta_all$Sample_Category %in% c("PCR_NTC", "Extraction_blank", "ASW_process_control")]
  ctrl14_ids <- setdiff(ctrl15_ids, "S-163")

  r_bio167 <- rowSums(otu_all[, bio167_ids, drop = FALSE])
  r_ctrl14 <- rowSums(otu_all[, ctrl14_ids, drop = FALSE])
  r_s163   <- as.numeric(otu_all[, "S-163"])

  only_s163_and_ctrl14 <- (r_s163 > 0) & (r_bio167 == 0) & (r_ctrl14 > 0)
  only_ctrl14_no_s163  <- (r_s163 == 0) & (r_bio167 == 0) & (r_ctrl14 > 0)
  only_s163_no_other   <- (r_s163 > 0) & (r_bio167 == 0) & (r_ctrl14 == 0)

  grp_reagent   <- is_final_contam | only_s163_and_ctrl14
  grp_ctrl_only <- only_ctrl14_no_s163 & !is_final_contam
  grp_s163_only <- only_s163_no_other
  grp_clean_bio <- !is_final_contam & (r_bio167 > 0)

  otu_bio_all   <- otu_all[, bio167_ids, drop = FALSE]
  bio_tot_all   <- colSums(otu_bio_all)
  ra_bio_all    <- sweep(otu_bio_all, 2, bio_tot_all, "/")
  bio_prev_01   <- rowSums(ra_bio_all >= 0.001)
  top2_bio      <- order(-r_bio167)[1:2]
  is_fp         <- (bio_prev_01 >= 1L & bio_prev_01 <= 3L) & !grp_reagent & (r_bio167 > 0)
  is_fp[top2_bio] <- FALSE

  otu_bio_clean <- otu_all[grp_clean_bio, bio167_ids, drop = FALSE]
  ra_bio_clean  <- sweep(otu_bio_clean, 2, pmax(1, colSums(otu_bio_clean)), "/")
  hi_depth_ids  <- bio167_ids[bio_tot_all >= MIN_READS_SAMPLE]
  top30_list    <- lapply(hi_depth_ids, function(bid) order(-ra_bio_clean[, bid])[1:min(30L, nrow(ra_bio_clean))])
  names(top30_list) <- hi_depth_ids

  fit_nnls_qp <- function(y, X) {
    K <- ncol(X)
    if (sum(y) == 0 || K == 0L) return(list(w = rep(0, K), fitted_cap = 0))
    active <- (y > 0) | (rowSums(X) > 0)
    ya <- as.numeric(y[active]); Xa <- X[active, , drop = FALSE]
    sol <- quadprog::solve.QP(
      Dmat = crossprod(Xa) + diag(1e-12, K),
      dvec = as.numeric(crossprod(Xa, ya)),
      Amat = diag(1, K), bvec = rep(0, K)
    )
    w <- pmax(0, sol$solution)
    list(w = w, fitted_cap = sum(pmin(ya, as.numeric(Xa %*% w))))
  }

  ctrl15_ord <- c(
    "NEG-211011", "NEA-211013", "NEG-211129", "NEG-211124", "NEB-211013",
    "NEG-211026", "NEA-211027", "NEG-211202", "NEB-211027", "NEG-211014",
    "NEA-211015", "NEB-211015", "NEG-211117", "S-164", "S-163"
  )
  ctrl15_ord <- ctrl15_ord[ctrl15_ord %in% ctrl15_ids]

  s1c_rows <- lapply(ctrl15_ord, function(cid) {
    cvec <- as.numeric(otu_all[, cid])
    tot  <- sum(cvec)
    r_reag <- sum(cvec[grp_reagent])
    r_ctrl <- sum(cvec[grp_ctrl_only])
    r_s163_c <- sum(cvec[grp_s163_only])
    r_bio  <- sum(cvec[grp_clean_bio])

    y_clean <- cvec[grp_clean_bio]
    y_ra    <- if (sum(y_clean) > 0) y_clean / sum(y_clean) else rep(0, length(y_clean))
    # Bray-Curtis similarity & top-30 Pearson correlation on clean bio reads against >=1,000-read bio samples
    bc_sim <- round(sapply(hi_depth_ids, function(bid) sum(pmin(y_ra, ra_bio_clean[, bid]))), 4)
    cor_p  <- round(sapply(hi_depth_ids, function(bid) {
      t30 <- top30_list[[bid]]
      x <- ra_bio_clean[t30, bid]; yv <- y_ra[t30]
      if (sd(x) > 0 && sd(yv) > 0) cor(x, yv, method = "pearson") else 0
    }), 4)
    fp_gt0 <- which(is_fp & cvec > 0)
    fp_ge5 <- which(is_fp & cvec >= 5)
    n_fp_ge5 <- if (length(fp_ge5) > 0) as.integer(colSums(ra_bio_all[fp_ge5, hi_depth_ids, drop = FALSE] >= 0.001)) else integer(length(hi_depth_ids))
    fp_rds   <- if (length(fp_gt0) > 0) as.integer(crossprod(cvec[fp_gt0], (ra_bio_all[fp_gt0, hi_depth_ids, drop = FALSE] >= 0.001) * 1)) else integer(length(hi_depth_ids))
    ord_bc <- hi_depth_ids[order(-bc_sim, -cor_p)]
    if (any(n_fp_ge5 > 0)) {
      ord_fp <- hi_depth_ids[order(-n_fp_ge5, -fp_rds, -bc_sim)]
      top_ids <- unique(c(ord_bc[1], ord_fp[1], ord_bc[1:3]))[1:3]
    } else {
      top_ids <- ord_bc[1:3]
    }
    top1_id <- top_ids[1]
    in_top3 <- rowSums(otu_bio_clean[, top_ids, drop = FALSE] > 0) > 0
    r_ov3   <- sum(y_clean[in_top3])
    nn3     <- fit_nnls_qp(y_clean, ra_bio_clean[, top_ids, drop = FALSE])

    data.frame(
      SequencingID = cid,
      Sample_Category = meta_all[cid, "Sample_Category"],
      Original_Identifier = meta_all[cid, "Original_Identifier"],
      Plate_Well = meta_all[cid, "Plate_Well"],
      Total_reads = as.integer(tot),
      n_ASVs = sum(cvec > 0),
      Reagent_Contam_reads = as.integer(r_reag),
      Reagent_Contam_Pct = round(100 * r_reag / tot, 2),
      CtrlOnly_Pct = round(100 * r_ctrl / tot, 2),
      S163Only_Pct = round(100 * r_s163_c / tot, 2),
      Tissue_Derived_reads = as.integer(r_bio),
      Tissue_Derived_Pct = round(100 * r_bio / tot, 2),
      Best_Matching_Sample = top1_id,
      Top3_Samples = paste(top_ids, collapse = ","),
      Overlap_Top3_Pct = round(100 * r_ov3 / tot, 2),
      NNLS_Top3_Pct = round(100 * nn3$fitted_cap / tot, 2),
      stringsAsFactors = FALSE
    )
  })
  write.csv(do.call(rbind, s1c_rows), file.path(DIR_INTERIM, "table_S1C_control_attribution.csv"), row.names = FALSE)
}

# 1. Table S1D: Seawater ASV_1, ASV_2, ASV_3, ASV_5 reads & ASV_1:ASV_2:ASV_3 ratios
host_asvs <- c("ASV_1", "ASV_2", "ASV_3", "ASV_5")
sw_rows <- lapply(sw_ids, function(sid) {
  cnts <- as.integer(otu[host_asvs, sid])
  ras  <- 100 * as.numeric(ra[host_asvs, sid])
  data.frame(
    Seawater_ID = sid, Location = as.character(sd[sid, "Location"]),
    Plate_Well = as.character(sd[sid, "Plate_Well"]),
    Total_reads = sum(otu[, sid]),
    ASV_1_reads = cnts[1], ASV_1_RA_pct = ras[1],
    ASV_2_reads = cnts[2], ASV_2_RA_pct = ras[2],
    ASV_3_reads = cnts[3], ASV_3_RA_pct = ras[3],
    ASV_5_reads = cnts[4], ASV_5_RA_pct = ras[4],
    Ratio_ASV1_to_ASV2 = if (cnts[2] > 0) cnts[1] / cnts[2] else NA_real_,
    Ratio_ASV1_to_ASV3 = if (cnts[3] > 0) cnts[1] / cnts[3] else NA_real_,
    stringsAsFactors = FALSE
  )
})
write.csv(do.call(rbind, sw_rows), file.path(DIR_INTERIM, "table_S1D_seawater_ratios.csv"), row.names = FALSE)

# 2. Table S5A & S5B: Bacterioplankton clade enrichment ratios
#    (Med_Ratio_zero1: no pseudocount when SW count > 0; 1-read floor only when SW count == 0)
clade_defs <- list(
  "Synechococcus CC9902"           = !is.na(tax$Genus)  & tax$Genus == "Synechococcus CC9902",
  "SAR11 clade"                    = !is.na(tax$Order)  & tax$Order == "Pelagibacterales",
  "Planktomarina"                  = !is.na(tax$Genus)  & tax$Genus == "Planktomarina",
  "SUP05 cluster"                  = !is.na(tax$Genus)  & tax$Genus == "SUP05 cluster",
  "Flavobacteriales marine groups" = !is.na(tax$Order)  & tax$Order == "Flavobacteriales" &
                                     ((!is.na(tax$Family) & grepl("marine group", tax$Family, ignore.case = TRUE)) |
                                      (!is.na(tax$Genus)  & grepl("marine group", tax$Genus,  ignore.case = TRUE)))
)

sw_loc_map <- setNames(sw_ids, as.character(sd[sw_ids, "Location"]))
sw_tot <- colSums(otu[, sw_ids, drop = FALSE])

s5a_rows <- list()
s5b_rows <- list()

for (cl_name in names(clade_defs)) {
  aids <- rownames(tax)[clade_defs[[cl_name]]]
  tis_sids <- c(z_ids, t_ids)
  tis_ra <- if (length(aids) > 0) colSums(ra[aids, tis_sids, drop = FALSE]) else rep(0, length(tis_sids))
  sw_cnt <- if (length(aids) > 0) colSums(otu[aids, sw_ids, drop = FALSE]) else rep(0, length(sw_ids))
  sw_ra_adj <- ifelse(sw_cnt > 0, sw_cnt / sw_tot, 1.0 / sw_tot)

  samp_loc <- as.character(sd[tis_sids, "Location"])
  matched_sw <- sw_ra_adj[sw_loc_map[samp_loc]]
  ratio_vec <- tis_ra / matched_sw

  df_cl <- data.frame(
    Clade = cl_name, SequencingID = tis_sids,
    Tissue = ifelse(sd[tis_sids, "Zooid_or_Tunic"] == "Z", "Zooid", "Tunic"),
    Species = as.character(sd[tis_sids, "Species"]),
    Tissue_RA_pct = 100 * tis_ra,
    Matched_SW_RA_pct = 100 * matched_sw,
    Enrichment_Ratio = ratio_vec,
    stringsAsFactors = FALSE
  )
  # Overall by tissue (Table S5A)
  for (tis in c("Zooid", "Tunic")) {
    sub <- df_cl[df_cl$Tissue == tis, ]
    s5a_rows[[length(s5a_rows) + 1L]] <- data.frame(
      Clade = cl_name, Tissue = tis, n_samples = nrow(sub),
      Median_Tissue_RA_pct = median(sub$Tissue_RA_pct),
      Median_SW_RA_pct = median(sub$Matched_SW_RA_pct),
      Median_Enrichment_Ratio = median(sub$Enrichment_Ratio),
      IQR_Enrichment_Ratio = IQR(sub$Enrichment_Ratio),
      stringsAsFactors = FALSE
    )
  }
  # Per species x tissue (Synechococcus CC9902 -> Table S5B)
  if (cl_name == "Synechococcus CC9902") {
    for (sp in SP_LEVELS) {
      for (tis in c("Zooid", "Tunic")) {
        sub <- df_cl[df_cl$Species == sp & df_cl$Tissue == tis, ]
        s5b_rows[[length(s5b_rows) + 1L]] <- data.frame(
          Species = sp, Tissue = tis, n_samples = nrow(sub),
          Median_Tissue_RA_pct = median(sub$Tissue_RA_pct),
          IQR_Tissue_RA_pct = IQR(sub$Tissue_RA_pct),
          Median_SW_RA_pct = median(sub$Matched_SW_RA_pct),
          Median_Enrichment_Ratio = median(sub$Enrichment_Ratio),
          IQR_Enrichment_Ratio = IQR(sub$Enrichment_Ratio),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

write.csv(do.call(rbind, s5a_rows), file.path(DIR_INTERIM, "table_S5A_bacterioplankton_enrichment.csv"), row.names = FALSE)
write.csv(do.call(rbind, s5b_rows), file.path(DIR_INTERIM, "table_S5B_synechococcus_by_species.csv"), row.names = FALSE)
cat("Saved Table S1D, Table S5A, and Table S5B intermediate CSVs.\n")
print(sessionInfo())
