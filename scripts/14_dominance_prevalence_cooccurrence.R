#!/usr/bin/env Rscript
# ==============================================================================
# 14_dominance_prevalence_cooccurrence.R — Tables S6B–C & S7A–C (§3.5–§3.6)
# ==============================================================================
# Purpose : Compute (1) top-3 dominant ASVs per species x tissue with minBoot 50
#           and best_assignment_minBoot30 labels & single-sample flags (Table S7A);
#           (2) pooled tissue top-10, prevalence thresholds, and strict cross-
#           species core check (Table S7B); (3) within-specimen co-occurrence
#           via CMH + 9,999 within-species permutations (Table S7C); and
#           (4) category Alphaproteobacteria share & top-5 class RA (Table S6B–C).
# Inputs  : config.R, output/intermediate/ascidian_canonical.rds,
#           data/metadata/asv_crosswalk.csv
# Outputs : output/intermediate/table_S7A_dominant.csv
#           output/intermediate/table_S7B_prevalence.csv
#           output/intermediate/table_S7C_cooccurrence.csv
#           output/intermediate/table_S6B_category_alpha.csv
#           output/intermediate/table_S6C_category_classes.csv
# Supports: §3.5–§3.6; Table S6B, Table S6C, Table S7A, Table S7B, Table S7C
# ==============================================================================

source("config.R")
set.seed(SEED)

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(parallel)
})

phy_all <- readRDS(PATH_PHY_CANON)
cw <- read.csv(PATH_ASV_CROSSWALK, stringsAsFactors = FALSE)

sd <- as(sample_data(phy_all), "data.frame")
sd$SequencingID <- rownames(sd)
sd$IB <- as.character(sd$Interactive_Border) == "Y"
sd$pair_id <- ifelse(sd$Zooid_or_Tunic == "Seawater",
                     as.character(sd$Original_Identifier),
                     sub("[ZT]$", "", as.character(sd$Original_Identifier)))

sd_nonib <- sd[sd$Zooid_or_Tunic %in% c("Z", "T") & !sd$IB & sd$Species %in% SP_LEVELS, ]
ptab <- table(sd_nonib$pair_id, sd_nonib$Zooid_or_Tunic)
comp52 <- sort(rownames(ptab)[ptab[, "Z"] == 1 & ptab[, "T"] == 1])
z_ids  <- sapply(comp52, function(p) sd_nonib$SequencingID[sd_nonib$pair_id == p & sd_nonib$Zooid_or_Tunic == "Z"])
t_ids  <- sapply(comp52, function(p) sd_nonib$SequencingID[sd_nonib$pair_id == p & sd_nonib$Zooid_or_Tunic == "T"])
pair_sp <- sapply(comp52, function(p) as.character(sd_nonib$Species[sd_nonib$pair_id == p][1]))
sw_ids <- sort(sd$SequencingID[sd$Zooid_or_Tunic == "Seawater"])
sids_104 <- c(z_ids, t_ids)
sids_107 <- c(z_ids, t_ids, sw_ids)

phy_107 <- prune_samples(sids_107, phy_all)
phy_107 <- prune_taxa(taxa_sums(phy_107) > 0, phy_107)
otu_107 <- as(otu_table(phy_107), "matrix")
if (!taxa_are_rows(phy_107)) otu_107 <- t(otu_107)
otu_107 <- otu_107[, sids_107, drop = FALSE]
ra_107  <- sweep(otu_107, 2, colSums(otu_107), "/")

ranks_6 <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
tax_b50 <- as.data.frame(tax_table(phy_107), stringsAsFactors = FALSE)[, ranks_6, drop = FALSE]
for (r in ranks_6) {
  v <- as.character(tax_b50[[r]])
  v[is.na(v) | trimws(v) == "" | v == "NA"] <- NA_character_
  tax_b50[[r]] <- v
}
rownames(tax_b50) <- rownames(otu_107)

# Secondary minBoot=30 labels from crosswalk (best_assignment_minBoot30)
tax_b30 <- tax_b50
m_cw <- match(rownames(tax_b30), cw$ASV_ID)
hit_cw <- which(!is.na(m_cw) & !is.na(cw$best_assignment_minBoot30[m_cw]) & nzchar(cw$best_assignment_minBoot30[m_cw]))
if (length(hit_cw) > 0) {
  sp_toks <- strsplit(cw$best_assignment_minBoot30[m_cw[hit_cw]], ";", fixed = TRUE)
  for (k in seq_along(hit_cw)) {
    toks <- trimws(sp_toks[[k]])
    toks[toks %in% c("", "NA")] <- NA_character_
    for (j in seq_along(ranks_6)) if (j <= length(toks)) tax_b30[hit_cw[k], ranks_6[j]] <- toks[j]
  }
}

fmt_plain <- function(v50, v30) {
  a <- ifelse(is.na(v50), "NA", v50)
  b <- ifelse(is.na(v30), "NA", v30)
  if (a == "NA" && b != "NA") sprintf("%s (bootstrap <50)", b)
  else if (a != b && b != "NA") sprintf("%s [%s (bootstrap <50)]", a, b)
  else a
}

fmt_taxon_label <- function(aid) {
  if (aid == "ASV_13") return("UBA10353 marine group")
  for (rk in c("Genus", "Family", "Order", "Class", "Phylum")) {
    lbl <- fmt_plain(tax_b50[aid, rk], tax_b30[aid, rk])
    if (lbl != "NA") return(lbl)
  }
  "Unclassified Bacteria"
}

# Category membership (presence >= 2 reads)
pa_107 <- otu_107 >= PRESENCE_MIN_COUNTS
inZ <- rowSums(pa_107[, z_ids, drop = FALSE]) > 0
inT <- rowSums(pa_107[, t_ids, drop = FALSE]) > 0
inSW <- rowSums(pa_107[, sw_ids, drop = FALSE]) > 0
asv_cat <- rep("Below_threshold", nrow(otu_107)); names(asv_cat) <- rownames(otu_107)
asv_cat[ inZ & !inT & !inSW] <- "Zooid-only"
asv_cat[!inZ &  inT & !inSW] <- "Tunic-only"
asv_cat[ inZ &  inT & !inSW] <- "Zooid–Tunic"
asv_cat[ inZ &  inT &  inSW] <- "Zooid–Tunic–Seawater"
asv_cat[!inZ & !inT &  inSW] <- "Seawater-only"
asv_cat[ inZ & !inT &  inSW] <- "Zooid–Seawater"
asv_cat[!inZ &  inT &  inSW] <- "Tunic–Seawater"

# 1. Table S7A: Dominant top-3 ASVs per species x tissue
dom_rows <- list()
for (sp in SP_LEVELS) {
  z_sp <- z_ids[pair_sp == sp]; t_sp <- t_ids[pair_sp == sp]
  med_Z <- apply(ra_107[, z_sp, drop = FALSE], 1, median); mean_Z <- rowMeans(ra_107[, z_sp, drop = FALSE])
  med_T <- apply(ra_107[, t_sp, drop = FALSE], 1, median); mean_T <- rowMeans(ra_107[, t_sp, drop = FALSE])
  ord_Z <- order(-med_Z, -mean_Z); ord_T <- order(-med_T, -mean_T)
  rk_Z <- integer(length(med_Z)); rk_Z[ord_Z] <- seq_along(ord_Z); names(rk_Z) <- names(med_Z)
  rk_T <- integer(length(med_T)); rk_T[ord_T] <- seq_along(ord_T); names(rk_T) <- names(med_T)
  for (tis in c("Zooid", "Tunic")) {
    su <- if (tis == "Zooid") z_sp else t_sp
    ou <- if (tis == "Zooid") ord_Z else ord_T
    or <- if (tis == "Zooid") rk_T else rk_Z
    om <- if (tis == "Zooid") med_T else med_Z
    for (k in 1:3) {
      aid <- rownames(ra_107)[ou[k]]
      vr <- as.numeric(ra_107[aid, su]); vc <- as.integer(otu_107[aid, su])
      ms <- if (sum(vr) > 0) max(vr) / sum(vr) else 0
      dom_rows[[length(dom_rows) + 1L]] <- data.frame(
        Species = sp, Species_full = unname(SP_FULL_LATIN[sp]), Tissue = tis,
        n_samples = length(su), Rank = k, ASV = aid,
        Taxon_label = fmt_taxon_label(aid),
        Class = ifelse(is.na(tax_b50[aid, "Class"]), "NA", tax_b50[aid, "Class"]),
        Order = fmt_plain(tax_b50[aid, "Order"], tax_b30[aid, "Order"]),
        Family = fmt_plain(tax_b50[aid, "Family"], tax_b30[aid, "Family"]),
        Genus = fmt_plain(tax_b50[aid, "Genus"], tax_b30[aid, "Genus"]),
        Median_RA_pct = 100 * median(vr),
        Q25_RA_pct = 100 * as.numeric(quantile(vr, 0.25)),
        Q75_RA_pct = 100 * as.numeric(quantile(vr, 0.75)),
        IQR_RA_pct = 100 * IQR(vr),
        Mean_RA_pct = 100 * mean(vr),
        Prev_ge2_n = sum(vc >= PRESENCE_MIN_COUNTS),
        Prev_ge2_pct = 100 * mean(vc >= PRESENCE_MIN_COUNTS),
        Rank_in_other_tissue = as.integer(or[aid]),
        Median_RA_other_tissue_pct = 100 * as.numeric(om[aid]),
        Max_single_sample_share_pct = 100 * ms,
        Single_sample_driven_flag = ms > 0.50,
        ASV_Category = unname(asv_cat[aid]),
        stringsAsFactors = FALSE
      )
    }
  }
}
write.csv(do.call(rbind, dom_rows), file.path(DIR_INTERIM, "table_S7A_dominant.csv"), row.names = FALSE)

# 2. Table S7B: Pooled tissue top-10 by median RA & by prevalence + threshold counts
pool_rows <- list()
thresh_rows <- list()
for (tis in c("Zooid", "Tunic")) {
  su <- if (tis == "Zooid") z_ids else t_ids
  sr <- ra_107[, su, drop = FALSE]; sc <- otu_107[, su, drop = FALSE]
  med_v <- apply(sr, 1, median); iqr_v <- apply(sr, 1, IQR); mean_v <- rowMeans(sr)
  pn <- rowSums(sc >= PRESENCE_MIN_COUNTS); pp <- 100 * pn / length(su)
  # Strict cross-species core check (>=90% prevalence in EVERY species)
  sp_prev_mat <- sapply(SP_LEVELS, function(sp) {
    s_sp <- su[pair_sp == sp]
    rowMeans(otu_107[, s_sp, drop = FALSE] >= PRESENCE_MIN_COUNTS) >= 0.90
  })
  n_strict_core <- sum(rowSums(sp_prev_mat) == length(SP_LEVELS))
  thresh_rows[[length(thresh_rows) + 1L]] <- data.frame(
    Tissue = tis, n_samples = length(su),
    n_ASVs_ge50pct = sum(pp >= 50),
    n_ASVs_ge75pct = sum(pp >= 75),
    n_ASVs_ge90pct = sum(pp >= 90),
    n_ASVs_strict_core_all6sp_ge90pct = n_strict_core,
    stringsAsFactors = FALSE
  )
  for (rk_type in c("Median_RA", "Prevalence")) {
    idx <- if (rk_type == "Median_RA") order(-med_v, -mean_v)[1:10] else order(-pn, -med_v, -mean_v)[1:10]
    aids <- rownames(sr)[idx]
    pool_rows[[length(pool_rows) + 1L]] <- data.frame(
      Tissue = tis, Ranking = rk_type, Rank = 1:10, ASV = aids,
      Taxon_label = sapply(aids, fmt_taxon_label),
      Family = sapply(aids, function(a) fmt_plain(tax_b50[a, "Family"], tax_b30[a, "Family"])),
      Genus  = sapply(aids, function(a) fmt_plain(tax_b50[a, "Genus"], tax_b30[a, "Genus"])),
      Median_RA_pct = 100 * med_v[aids], IQR_RA_pct = 100 * iqr_v[aids],
      Prev_ge2_n = pn[aids], Prev_ge2_pct = pp[aids],
      stringsAsFactors = FALSE
    )
  }
}
write.csv(do.call(rbind, pool_rows), file.path(DIR_INTERIM, "table_S7B_prevalence.csv"), row.names = FALSE)
thresh_df <- do.call(rbind, thresh_rows)
write.csv(thresh_df, file.path(DIR_INTERIM, "table_S7B_prevalence_summary.csv"), row.names = FALSE)
cat("Table S7B prevalence threshold counts (>=50% / >=75% / >=90%):\n")
print(thresh_df, row.names = FALSE)

# 3. Table S7C: Within-specimen co-occurrence (CMH + 9,999 within-species permutations, >=2 counts)
run_cmh_safe <- function(zp, tp, sp) {
  arr <- table(Z = factor(zp, levels = c(FALSE, TRUE)), T = factor(tp, levels = c(FALSE, TRUE)), Species = sp)
  if (sum(apply(arr, 3, function(m) all(rowSums(m) > 0) && all(colSums(m) > 0))) == 0) {
    return(list(n_info = 0L, p = NA_real_, or = NA_real_))
  }
  res <- tryCatch(mantelhaen.test(arr), error = function(e) e)
  if (inherits(res, "error") || is.nan(res$p.value)) return(list(n_info = 0L, p = NA_real_, or = NA_real_))
  list(n_info = sum(apply(arr, 3, function(m) all(rowSums(m) > 0) && all(colSums(m) > 0))),
       p = res$p.value, or = unname(res$estimate))
}

set.seed(SEED)
n_pairs <- length(comp52)
sp_vec  <- factor(pair_sp, levels = SP_LEVELS)
sp_idx  <- split(seq_len(n_pairs), sp_vec, drop = TRUE)
Z_pa <- (otu_107[, z_ids, drop = FALSE] >= PRESENCE_MIN_COUNTS) * 1.0
T_pa <- (otu_107[, t_ids, drop = FALSE] >= PRESENCE_MIN_COUNTS) * 1.0
pair_prev <- rowSums((Z_pa + T_pa) > 0)
keep <- pair_prev >= COOCCUR_MIN_PAIRS
Z_k <- Z_pa[keep, , drop = FALSE]; T_k <- T_pa[keep, , drop = FALSE]
n_test <- nrow(Z_k)

perm_mat <- matrix(seq_len(n_pairs), nrow = n_pairs, ncol = N_PERM)
for (b in seq_len(N_PERM)) {
  p_idx <- seq_len(n_pairs)
  for (idx_s in sp_idx) if (length(idx_s) > 1) p_idx[idx_s] <- sample(idx_s)
  perm_mat[, b] <- p_idx
}
obs_nb <- rowSums(Z_k * T_k)
chunk_idx <- split(seq_len(N_PERM), cut(seq_len(N_PERM), breaks = N_CORES, labels = FALSE))
perm_chunks <- mclapply(chunk_idx, function(b_vec) {
  ge <- integer(n_test); sm <- numeric(n_test)
  for (b in b_vec) {
    nb <- rowSums(Z_k * T_k[, perm_mat[, b], drop = FALSE])
    ge <- ge + (nb >= obs_nb); sm <- sm + nb
  }
  list(ge = ge, sm = sm)
}, mc.cores = N_CORES)
perm_p <- (1.0 + Reduce("+", lapply(perm_chunks, `[[`, "ge"))) / (N_PERM + 1.0)
perm_q <- p.adjust(perm_p, method = "BH")
null_mean_nb <- Reduce("+", lapply(perm_chunks, `[[`, "sm")) / N_PERM

cmh_mat <- do.call(rbind, mclapply(seq_len(n_test), function(i) {
  cm <- run_cmh_safe(Z_k[i, ] > 0, T_k[i, ] > 0, sp_vec)
  c(n_info = cm$n_info, cmh_p = cm$p, cmh_or = cm$or)
}, mc.cores = N_CORES))

aids <- rownames(Z_k)
cooc_df <- data.frame(
  ASV = aids,
  Order  = sapply(aids, function(a) fmt_plain(tax_b50[a, "Order"], tax_b30[a, "Order"])),
  Family = sapply(aids, function(a) fmt_plain(tax_b50[a, "Family"], tax_b30[a, "Family"])),
  Genus  = sapply(aids, function(a) fmt_plain(tax_b50[a, "Genus"], tax_b30[a, "Genus"])),
  n_pairs_present = as.integer(pair_prev[keep]),
  n_both   = as.integer(obs_nb),
  n_Z_only = as.integer(rowSums(Z_k > 0 & T_k == 0)),
  n_T_only = as.integer(rowSums(Z_k == 0 & T_k > 0)),
  null_mean_both = null_mean_nb,
  cmh_n_info = as.integer(cmh_mat[, "n_info"]),
  cmh_OR = cmh_mat[, "cmh_or"],
  cmh_p  = cmh_mat[, "cmh_p"],
  cmh_q  = p.adjust(cmh_mat[, "cmh_p"], method = "BH"),
  perm_p = perm_p,
  perm_q = perm_q,
  stringsAsFactors = FALSE
)
cooc_df$sig_cmh  <- !is.na(cooc_df$cmh_q) & cooc_df$cmh_q < ALPHA_SIG
cooc_df$sig_perm <- cooc_df$perm_q < ALPHA_SIG
cooc_df$sig_both <- cooc_df$sig_cmh & cooc_df$sig_perm
cooc_df <- cooc_df[order(!cooc_df$sig_both, !cooc_df$sig_perm, cooc_df$perm_q, cooc_df$cmh_q), ]
write.csv(cooc_df[cooc_df$sig_cmh | cooc_df$sig_perm, ],
          file.path(DIR_INTERIM, "table_S7C_cooccurrence.csv"), row.names = FALSE)

# 4. Table S6B & S6C: Category Alphaproteobacteria share & top-5 classes
cats_4 <- c("Zooid–Tunic", "Zooid–Tunic–Seawater", "Zooid-only", "Tunic-only")
cls_vec <- ifelse(is.na(tax_b50$Class), "Unclassified", tax_b50$Class)
names(cls_vec) <- rownames(tax_b50)
cls_ra_104 <- rowsum(ra_107[, sids_104, drop = FALSE], group = cls_vec, reorder = TRUE)
top5_cls <- names(sort(rowMeans(cls_ra_104), decreasing = TRUE))[1:5]
cls_grp <- ifelse(cls_vec %in% top5_cls, cls_vec, "Other")
names(cls_grp) <- names(cls_vec)

alpha_rows <- list(); class_rows <- list()
for (ct in cats_4) {
  ac <- names(asv_cat)[asv_cat == ct]
  na <- sum(cls_vec[ac] == "Alphaproteobacteria")
  rt <- sum(otu_107[ac, sids_104, drop = FALSE])
  rta <- sum(otu_107[ac[cls_vec[ac] == "Alphaproteobacteria"], sids_104, drop = FALSE])
  r107 <- sum(otu_107[ac, sids_107, drop = FALSE])
  r107a <- sum(otu_107[ac[cls_vec[ac] == "Alphaproteobacteria"], sids_107, drop = FALSE])
  alpha_rows[[length(alpha_rows) + 1L]] <- data.frame(
    Category = ct, n_ASVs = length(ac), n_Alpha_ASVs = na, Pct_Alpha_ASVs = 100 * na / length(ac),
    Reads_tissue_104 = rt, Reads_Alpha_tissue_104 = rta, Pct_Alpha_reads_tissue_104 = 100 * rta / rt,
    Reads_all_107 = r107, Pct_Alpha_reads_all_107 = 100 * r107a / r107, stringsAsFactors = FALSE
  )
  for (tis in c("Zooid", "Tunic")) {
    su <- if (tis == "Zooid") z_ids else t_ids
    tot_c <- colSums(ra_107[ac, su, drop = FALSE])
    for (cl in c(top5_cls, "Other")) {
      ca <- ac[cls_grp[ac] == cl]
      vc <- if (length(ca) > 0) colSums(ra_107[ca, su, drop = FALSE]) else rep(0, length(su))
      class_rows[[length(class_rows) + 1L]] <- data.frame(
        Category = ct, Tissue = tis, n_samples = length(su), Class = cl,
        Median_RA_pct = 100 * median(vc),
        Q25_RA_pct = 100 * as.numeric(quantile(vc, 0.25)),
        Q75_RA_pct = 100 * as.numeric(quantile(vc, 0.75)),
        IQR_RA_pct = 100 * IQR(vc),
        Mean_RA_pct = 100 * mean(vc),
        Category_Total_Median_RA_pct = 100 * median(tot_c),
        stringsAsFactors = FALSE
      )
    }
  }
}
write.csv(do.call(rbind, alpha_rows), file.path(DIR_INTERIM, "table_S6B_category_alpha.csv"), row.names = FALSE)
write.csv(do.call(rbind, class_rows), file.path(DIR_INTERIM, "table_S6C_category_classes.csv"), row.names = FALSE)
cat("Saved Tables S6B, S6C, S7A, S7B, S7C intermediate CSVs.\n")
print(sessionInfo())
