#!/usr/bin/env Rscript
# ==============================================================================
# 11_categories_and_null.R — ASV Categories, 999-Swap Null & Table 1 / S6A
# ==============================================================================
# Purpose : Classify ASVs into 7 detection categories across 52 complete non-IB
#           Zooid-Tunic pairs + 3 Seawater samples (n = 107, presence >= 2 reads),
#           run the 999-permutation within-pair Z/T label-swap null (set.seed(42)),
#           compute Table 1 (counts, per-sample median [IQR], Mean_RA_tissue_pct),
#           Table S6A (per-species category medians), and the 96 Z-T-SW breakdown.
# Inputs  : config.R, output/intermediate/ascidian_canonical.rds
# Outputs : output/tables/Table1_asv_categories.csv
#           output/intermediate/table_S6A_categories_by_species.csv
#           output/intermediate/fig2A_category_ra.csv
#           output/intermediate/fig2A_swap_null.csv
#           output/intermediate/fig2B_nmds_coords.csv
#           output/intermediate/ztsw_96_summary.csv
# Supports: Table 1; Figure 2 (A-C); Table S6A; §3.2
# ==============================================================================

source("config.R")
set.seed(SEED)

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(parallel)
})

phy_all <- readRDS(PATH_PHY_CANON)
sd <- as(sample_data(phy_all), "data.frame")
sd$SequencingID <- rownames(sd)
sd$IB <- as.character(sd$Interactive_Border) == "Y"
sd$pair_id <- ifelse(sd$Zooid_or_Tunic == "Seawater",
                     as.character(sd$Original_Identifier),
                     sub("[ZT]$", "", as.character(sd$Original_Identifier)))

sd_nonib <- sd[sd$Zooid_or_Tunic %in% c("Z", "T") & !sd$IB & sd$Species %in% SP_LEVELS, ]
ptab <- table(sd_nonib$pair_id, sd_nonib$Zooid_or_Tunic)
comp52 <- sort(rownames(ptab)[ptab[, "Z"] == 1 & ptab[, "T"] == 1])
stopifnot(length(comp52) == 52L)

z_ids  <- sapply(comp52, function(p) sd_nonib$SequencingID[sd_nonib$pair_id == p & sd_nonib$Zooid_or_Tunic == "Z"])
t_ids  <- sapply(comp52, function(p) sd_nonib$SequencingID[sd_nonib$pair_id == p & sd_nonib$Zooid_or_Tunic == "T"])
sw_ids <- sort(sd$SequencingID[sd$Zooid_or_Tunic == "Seawater"])
sids_107 <- c(z_ids, t_ids, sw_ids)

phy_107 <- prune_samples(sids_107, phy_all)
phy_107 <- prune_taxa(taxa_sums(phy_107) > 0, phy_107)
otu <- as(otu_table(phy_107), "matrix")
if (!taxa_are_rows(phy_107)) otu <- t(otu)
otu <- otu[, sids_107, drop = FALSE]
ra  <- sweep(otu, 2, colSums(otu), "/")

pa  <- otu >= PRESENCE_MIN_COUNTS
inZ  <- rowSums(pa[, z_ids, drop = FALSE]) > 0
inT  <- rowSums(pa[, t_ids, drop = FALSE]) > 0
inSW <- rowSums(pa[, sw_ids, drop = FALSE]) > 0

cat_levels <- c("Zooid-only", "Tunic-only", "Zooid–Tunic", "Zooid–Tunic–Seawater",
                "Seawater-only", "Zooid–Seawater", "Tunic–Seawater")
asv_cat <- rep("Below_threshold", nrow(otu))
names(asv_cat) <- rownames(otu)
asv_cat[ inZ & !inT & !inSW] <- "Zooid-only"
asv_cat[!inZ &  inT & !inSW] <- "Tunic-only"
asv_cat[ inZ &  inT & !inSW] <- "Zooid–Tunic"
asv_cat[ inZ &  inT &  inSW] <- "Zooid–Tunic–Seawater"
asv_cat[!inZ & !inT &  inSW] <- "Seawater-only"
asv_cat[ inZ & !inT &  inSW] <- "Zooid–Seawater"
asv_cat[!inZ &  inT &  inSW] <- "Tunic–Seawater"

# Per-sample summed RA by category
per_samp_list <- lapply(cat_levels, function(ct) {
  idx <- which(asv_cat == ct)
  vals <- if (length(idx) > 0) 100 * colSums(ra[idx, , drop = FALSE]) else rep(0, ncol(ra))
  data.frame(
    SequencingID = sids_107,
    Sample_Type  = c(rep("Zooid", 52), rep("Tunic", 52), rep("Seawater", 3)),
    Species      = as.character(sd[sids_107, "Species"]),
    ASV_Category = ct,
    RA_pct       = as.numeric(vals),
    stringsAsFactors = FALSE
  )
})
per_samp_df <- do.call(rbind, per_samp_list)
write.csv(per_samp_df, file.path(DIR_INTERIM, "fig2A_category_ra.csv"), row.names = FALSE)

# Table 1 summary (including Mean_RA_tissue_pct across 104 tissue samples)
t1_rows <- lapply(cat_levels, function(ct) {
  idx <- which(asv_cat == ct)
  sub <- per_samp_df[per_samp_df$ASV_Category == ct, ]
  vz  <- sub$RA_pct[sub$Sample_Type == "Zooid"]
  vt  <- sub$RA_pct[sub$Sample_Type == "Tunic"]
  vsw <- sub$RA_pct[sub$Sample_Type == "Seawater"]
  vtis <- c(vz, vt)
  data.frame(
    ASV_Category = ct,
    n_ASVs = length(idx),
    Summed_reads_107 = if (length(idx) > 0) sum(otu[idx, ]) else 0L,
    Mean_RA_tissue_pct = mean(vtis),
    Zooid_median_pct = median(vz), Zooid_IQR_pct = IQR(vz),
    Tunic_median_pct = median(vt), Tunic_IQR_pct = IQR(vt),
    Seawater_median_pct = median(vsw), Seawater_IQR_pct = IQR(vsw),
    stringsAsFactors = FALSE
  )
})
t1_df <- do.call(rbind, t1_rows)
write.csv(t1_df, file.path(DIR_TABLES, "Table1_asv_categories.csv"), row.names = FALSE)

# 999-swap within-pair Z/T label-swap null (exact set.seed(42) matching make_fig2_v3.R)
set.seed(SEED)
swaps <- matrix(runif(52L * N_SWAP_NULL) < 0.5, nrow = 52L, ncol = N_SWAP_NULL)
z_pos <- seq_len(52L); t_pos <- 52L + seq_len(52L)
pa_num <- pa * 1.0

eval_swap <- function(fl) {
  cz <- ifelse(fl, t_pos, z_pos)
  ct <- ifelse(fl, z_pos, t_pos)
  iz <- rowSums(pa_num[, cz, drop = FALSE]) > 0
  it <- rowSums(pa_num[, ct, drop = FALSE]) > 0
  zo <- iz & !it & !inSW
  to <- !iz & it & !inSW
  c(
    n_Z_only = sum(zo),
    n_T_only = sum(to),
    med_Z_only_in_Z = median(100 * colSums(ra[zo, cz, drop = FALSE])),
    med_T_only_in_T = median(100 * colSums(ra[to, ct, drop = FALSE]))
  )
}
obs_null <- eval_swap(rep(FALSE, 52L))
null_mat <- do.call(rbind, mclapply(seq_len(N_SWAP_NULL), function(b) eval_swap(swaps[, b]), mc.cores = N_CORES))

null_summary <- do.call(rbind, lapply(names(obs_null), function(m) {
  ov <- unname(obs_null[m]); nv <- null_mat[, m]
  p_lo <- (1 + sum(nv <= ov)) / (N_SWAP_NULL + 1)
  p_hi <- (1 + sum(nv >= ov)) / (N_SWAP_NULL + 1)
  data.frame(
    Metric = m, Observed = ov,
    Null_median = median(nv),
    Null_q025 = unname(quantile(nv, 0.025)),
    Null_q975 = unname(quantile(nv, 0.975)),
    p_two_sided = min(1, 2 * min(p_lo, p_hi)),
    stringsAsFactors = FALSE
  )
}))
write.csv(null_summary, file.path(DIR_INTERIM, "fig2A_swap_null.csv"), row.names = FALSE)

# Table S6A: per-species category contributions (Zooid and Tunic)
s6a_df <- per_samp_df %>%
  filter(Sample_Type %in% c("Zooid", "Tunic")) %>%
  group_by(Species, Sample_Type, ASV_Category) %>%
  summarise(
    n_Samples  = n(),
    Median_pct = median(RA_pct),
    Q25_pct    = as.numeric(quantile(RA_pct, 0.25)),
    Q75_pct    = as.numeric(quantile(RA_pct, 0.75)),
    IQR_pct    = IQR(RA_pct),
    .groups    = "drop"
  )
write.csv(s6a_df, file.path(DIR_INTERIM, "table_S6A_categories_by_species.csv"), row.names = FALSE)

# Bray-Curtis NMDS for Figure 2B (k = 2, trymax = 500, maxit = 999, set.seed(42))
set.seed(SEED)
d_bc <- vegdist(t(ra), method = "bray")
nmds <- metaMDS(d_bc, k = 2, trymax = NMDS_TRYMAX, maxit = NMDS_MAXIT, trace = FALSE)
nmds_df <- data.frame(
  SequencingID = sids_107,
  Sample_Type  = c(rep("Zooid", 52), rep("Tunic", 52), rep("Seawater", 3)),
  Species      = as.character(sd[sids_107, "Species"]),
  NMDS1        = scores(nmds, display = "sites")[, 1],
  NMDS2        = scores(nmds, display = "sites")[, 2],
  Stress       = nmds$stress,
  Tries        = nmds$tries,
  stringsAsFactors = FALSE
)
write.csv(nmds_df, file.path(DIR_INTERIM, "fig2B_nmds_coords.csv"), row.names = FALSE)

cat("Table 1 category counts:\n")
print(t1_df[, c("ASV_Category", "n_ASVs", "Mean_RA_tissue_pct")], row.names = FALSE)
print(sessionInfo())
