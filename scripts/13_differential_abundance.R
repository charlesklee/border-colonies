#!/usr/bin/env Rscript
# ==============================================================================
# 13_differential_abundance.R — Differential Abundance (Z vs T; Table S8, Fig 4)
# ==============================================================================
# Purpose : Run three complementary DA frameworks across 52 complete non-IB
#           Zooid-Tunic pairs (>=10% prevalence per level: Order, Family, Genus,
#           ASV) — ANCOM-BC2 (species + tissue + (1|pair_id), BH q<0.05 + ss=0),
#           paired ALDEx2 (128 MC, denom='all', wi.eBH & we.eBH), and MaAsLin2
#           (TSS_LOG, TSS_AST, CLR_NONE, tissue-only BH q<0.05) — plus ASV_24,
#           Cyanobiaceae, Vibrionaceae, and ASV_8 (in DV only) within-pair log2
#           ratios and Figure 4A/4B source data.
# Inputs  : config.R, output/intermediate/ascidian_canonical.rds,
#           data/metadata/asv_crosswalk.csv
# Outputs : output/intermediate/table_S8_differential_abundance.csv
#           output/intermediate/fig4A_top15_families.csv
#           output/intermediate/fig4B_within_pair_log2.csv
# Supports: Figure 4 (A-B); Table S8; §3.4
# ==============================================================================

source("config.R")
set.seed(SEED)

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(ALDEx2)
  library(Maaslin2)
  library(dplyr)
})

phy_all <- readRDS(PATH_PHY_CANON)
sd_all  <- as(sample_data(phy_all), "data.frame")
sd_all$SequencingID <- rownames(sd_all)
sd_all$IB <- as.character(sd_all$Interactive_Border) == "Y"
sd_all$pair_id <- ifelse(sd_all$Zooid_or_Tunic == "Seawater", as.character(sd_all$Original_Identifier),
                         sub("[ZT]$", "", as.character(sd_all$Original_Identifier)))

sd_nonib <- sd_all[sd_all$Zooid_or_Tunic %in% c("Z", "T") & !sd_all$IB & sd_all$Species %in% SP_LEVELS, ]
ptab <- table(sd_nonib$pair_id, sd_nonib$Zooid_or_Tunic)
comp52 <- sort(rownames(ptab)[ptab[, "Z"] == 1 & ptab[, "T"] == 1])
stopifnot(length(comp52) == 52L)

t_ids <- sapply(comp52, function(p) sd_nonib$SequencingID[sd_nonib$pair_id == p & sd_nonib$Zooid_or_Tunic == "T"])
z_ids <- sapply(comp52, function(p) sd_nonib$SequencingID[sd_nonib$pair_id == p & sd_nonib$Zooid_or_Tunic == "Z"])
pair_sp <- sapply(comp52, function(p) as.character(sd_nonib$Species[sd_nonib$pair_id == p][1]))
samp_order <- c(t_ids, z_ids)

phy_52 <- prune_samples(samp_order, phy_all)
phy_52 <- prune_taxa(taxa_sums(phy_52) > 0, phy_52)
otu_cnt <- as(otu_table(phy_52), "matrix")
if (!taxa_are_rows(phy_52)) otu_cnt <- t(otu_cnt)
otu_cnt <- otu_cnt[, samp_order, drop = FALSE]
otu_ra  <- sweep(otu_cnt, 2, colSums(otu_cnt), "/")
p_pseudo <- 0.5 * min(otu_ra[otu_ra > 0])

fill_tax <- function(tt) {
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  out <- as.data.frame(tt[, ranks, drop = FALSE], stringsAsFactors = FALSE)
  for (i in seq_along(ranks)) {
    r <- ranks[i]; vals <- as.character(out[[r]])
    na_idx <- is.na(vals) | trimws(vals) == ""
    if (i == 1) vals[na_idx] <- "Unclassified"
    else {
      pv <- as.character(out[[ranks[i - 1]]])
      vals[na_idx] <- ifelse(grepl("^Unclassified_", pv[na_idx]), pv[na_idx], paste0("Unclassified_", pv[na_idx]))
    }
    out[[r]] <- vals
  }
  out
}
agg_rank <- function(m, tt, rk) {
  grp <- if (rk == "Order") paste(tt$Phylum, tt$Class, tt$Order, sep = " | ")
         else if (rk == "Family") paste(tt$Phylum, tt$Class, tt$Order, tt$Family, sep = " | ")
         else paste(tt$Phylum, tt$Class, tt$Order, tt$Family, tt$Genus, sep = " | ")
  rowsum(m, group = grp, reorder = TRUE)
}

tax_b50 <- fill_tax(as.data.frame(tax_table(phy_52), stringsAsFactors = FALSE))
rownames(tax_b50) <- rownames(otu_cnt)

meta_ord <- data.frame(
  Species = factor(c(pair_sp, pair_sp), levels = SP_LEVELS),
  Tissue  = factor(c(rep("T", 52), rep("Z", 52)), levels = c("T", "Z")),
  pair_id = factor(c(comp52, comp52)),
  row.names = samp_order
)

run_da_level <- function(cmat, lv_name, tt_df = NULL) {
  set.seed(SEED)
  cmat <- cmat[rowSums(cmat) > 0, , drop = FALSE]
  cmat_10 <- cmat[rowSums(cmat > 0) >= ceiling(DA_MIN_PREV * ncol(cmat)), , drop = FALSE]

  ps_in <- phyloseq(otu_table(cmat, taxa_are_rows = TRUE), sample_data(meta_ord))
  anc <- suppressMessages(suppressWarnings(ANCOMBC::ancombc2(
    data = ps_in, assay_name = "counts", tax_level = NULL,
    fix_formula = "Species + Tissue", rand_formula = "(1 | pair_id)",
    p_adj_method = "BH", pseudo = 0, pseudo_sens = TRUE,
    prv_cut = DA_MIN_PREV, lib_cut = 0, s0_perc = 0.05, group = "Tissue",
    struc_zero = TRUE, neg_lb = FALSE, alpha = ALPHA_SIG, n_cl = N_CL_DA, verbose = FALSE
  )))
  res_anc <- anc$res
  res_anc$passed_ss <- !is.null(anc$pseudo_sens_tab) && "TissueZ" %in% colnames(anc$pseudo_sens_tab) &
                       (anc$pseudo_sens_tab$TissueZ[match(res_anc$taxon, anc$pseudo_sens_tab$taxon)] == 0)

  set.seed(SEED)
  ald_p <- suppressMessages(ALDEx2::aldex(
    cmat_10, conditions = as.character(meta_ord$Tissue),
    mc.samples = N_ALDEX_MC, test = "t", effect = TRUE, denom = "all", paired.test = TRUE
  ))
  ald_p$taxon <- rownames(ald_p)

  safe_nm <- sprintf("F%05d", seq_len(nrow(cmat)))
  nm_map  <- setNames(rownames(cmat), safe_nm)
  cmat_m2 <- cmat; rownames(cmat_m2) <- safe_nm
  df_in   <- as.data.frame(t(cmat_m2))
  run_m2 <- function(nrm, trn) {
    td <- tempfile(pattern = paste0("m2_", nrm, "_"))
    dir.create(td, recursive = TRUE, showWarnings = FALSE)
    fit <- suppressMessages(suppressWarnings({
      capture.output({
        m2 <- Maaslin2::Maaslin2(
          input_data = df_in, input_metadata = meta_ord, output = td,
          min_abundance = 0, min_prevalence = DA_MIN_PREV, min_variance = 0,
          normalization = nrm, transform = trn, analysis_method = "LM",
          max_significance = 0.25, random_effects = "pair_id",
          fixed_effects = c("Species", "Tissue"), correction = "BH",
          standardize = TRUE, cores = N_CL_DA, plot_heatmap = FALSE, plot_scatter = FALSE,
          reference = "Species,AK;Tissue,T"
        )
      }); m2
    }))
    unlink(td, recursive = TRUE)
    rt <- fit$results[fit$results$metadata == "Tissue" & fit$results$value == "Z", , drop = FALSE]
    rt$Taxon <- unname(nm_map[rt$feature])
    rt$q_tis <- p.adjust(rt$pval, method = "BH")
    rt
  }
  m2_log <- run_m2("TSS", "LOG")
  m2_ast <- run_m2("TSS", "AST")
  m2_clr <- run_m2("CLR", "NONE")

  log_m <- log(cmat + 0.5); clr_m <- sweep(log_m, 2, colMeans(log_m), "-")
  diff_clr <- clr_m[, z_ids, drop = FALSE] - clr_m[, t_ids, drop = FALSE]
  all_tx <- unique(c(res_anc$taxon, ald_p$taxon, m2_log$Taxon, m2_ast$Taxon, m2_clr$Taxon))

  rows <- lapply(all_tx, function(tx) {
    ia <- match(tx, res_anc$taxon); il <- match(tx, ald_p$taxon)
    ig <- match(tx, m2_log$Taxon);  is <- match(tx, m2_ast$Taxon); ic <- match(tx, m2_clr$Taxon)
    anc_q <- if (!is.na(ia)) res_anc$q_TissueZ[ia] else NA_real_
    anc_ss <- if (!is.na(ia)) isTRUE(res_anc$passed_ss[ia]) else FALSE
    wi_q <- if (!is.na(il)) ald_p$wi.eBH[il] else NA_real_
    we_q <- if (!is.na(il)) ald_p$we.eBH[il] else NA_real_
    lq <- if (!is.na(ig)) m2_log$q_tis[ig] else NA_real_
    aq <- if (!is.na(is)) m2_ast$q_tis[is] else NA_real_
    cq <- if (!is.na(ic)) m2_clr$q_tis[ic] else NA_real_
    s_anc <- !is.na(anc_q) && anc_q < ALPHA_SIG
    n_t <- as.integer(s_anc) + as.integer(s_anc && anc_ss) +
           as.integer(!is.na(wi_q) && wi_q < ALPHA_SIG) + as.integer(!is.na(we_q) && we_q < ALPHA_SIG) +
           as.integer(!is.na(lq) && lq < ALPHA_SIG) + as.integer(!is.na(aq) && aq < ALPHA_SIG) +
           as.integer(!is.na(cq) && cq < ALPHA_SIG)
    if (n_t == 0L) return(NULL)
    n_fw <- as.integer(s_anc) + as.integer((!is.na(wi_q) && wi_q < ALPHA_SIG) || (!is.na(we_q) && we_q < ALPHA_SIG)) +
            as.integer((!is.na(lq) && lq < ALPHA_SIG) || (!is.na(aq) && aq < ALPHA_SIG) || (!is.na(cq) && cq < ALPHA_SIG))
    dir_e <- if (!is.na(ia)) ifelse(res_anc$lfc_TissueZ[ia] > 0, "Z-enriched", "T-enriched")
             else ifelse(mean(diff_clr[tx, ]) > 0, "Z-enriched", "T-enriched")
    sp_d <- sapply(SP_LEVELS, function(sp) ifelse(mean(diff_clr[tx, pair_sp == sp]) > 0, "Z", "T"))
    tx_lbl <- if (lv_name == "ASV" && !is.null(tt_df)) sprintf("%s (%s | %s | %s)", tx, tt_df[tx, "Order"], tt_df[tx, "Family"], tt_df[tx, "Genus"]) else tx
    data.frame(
      Level = lv_name, Taxon_ID = tx, Taxon = tx_lbl, Direction = dir_e,
      N_tests_sig = n_t, N_frameworks_sig = n_fw,
      ANCOMBC2_LFC = if (!is.na(ia)) res_anc$lfc_TissueZ[ia] else NA_real_,
      ANCOMBC2_q = anc_q, ANCOMBC2_passed_ss = anc_ss,
      ALDEx2_effect = if (!is.na(il)) ald_p$effect[il] else NA_real_,
      ALDEx2_wi_q = wi_q, ALDEx2_we_q = we_q,
      MaAsLin2_CLR_coef = if (!is.na(ic)) m2_clr$coef[ic] else NA_real_,
      MaAsLin2_CLR_q = cq, MaAsLin2_TSS_LOG_q = lq, MaAsLin2_TSS_AST_q = aq,
      Sp_Consistency = sprintf("%d/6", sum(sp_d == ifelse(dir_e == "Z-enriched", "Z", "T"))),
      Per_Sp_Dir = paste(sp_d, collapse = "/"), stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(-out$N_frameworks_sig, -out$N_tests_sig, out$ANCOMBC2_q), ]
}

da_s8 <- rbind(
  run_da_level(agg_rank(otu_cnt, tax_b50, "Order"),  "Order"),
  run_da_level(agg_rank(otu_cnt, tax_b50, "Family"), "Family"),
  run_da_level(agg_rank(otu_cnt, tax_b50, "Genus"),  "Genus"),
  run_da_level(otu_cnt, "ASV", tax_b50)
)
write.csv(da_s8, file.path(DIR_INTERIM, "table_S8_differential_abundance.csv"), row.names = FALSE)

# Figure 4A & 4B source data
fam_ra <- rowsum(otu_ra, group = tax_b50$Family, reorder = TRUE)
top15  <- names(sort(rowMeans(fam_ra), decreasing = TRUE))[1:15]
f4a_rows <- list()
for (sp in SP_LEVELS) {
  for (tis in c("Z", "T")) {
    su <- if (tis == "Z") z_ids[pair_sp == sp] else t_ids[pair_sp == sp]
    mv <- rowMeans(fam_ra[, su, drop = FALSE])
    vals <- c(mv[top15], Other = sum(mv[setdiff(names(mv), top15)]))
    f4a_rows[[length(f4a_rows) + 1L]] <- data.frame(
      Species = sp, Tissue = ifelse(tis == "Z", "Zooid", "Tunic"), Tissue_abbr = tis,
      Family_raw = names(vals), Family = gsub("^Unclassified_", "Uncl. ", names(vals)),
      Mean_RA_pct = 100 * as.numeric(vals), stringsAsFactors = FALSE
    )
  }
}
write.csv(do.call(rbind, f4a_rows), file.path(DIR_INTERIM, "fig4A_top15_families.csv"), row.names = FALSE)

mk_l2 <- function(vz, vt, tx, sps = SP_LEVELS) {
  k <- pair_sp %in% sps
  data.frame(
    Taxon = tx, pair_id = comp52[k], Species = pair_sp[k],
    Z_RA = as.numeric(vz[k]), T_RA = as.numeric(vt[k]),
    log2_ratio = log2((as.numeric(vz[k]) + p_pseudo) / (as.numeric(vt[k]) + p_pseudo)),
    Pseudocount_p = p_pseudo, stringsAsFactors = FALSE
  )
}
f4b_df <- rbind(
  mk_l2(fam_ra["Cyanobiaceae", z_ids], fam_ra["Cyanobiaceae", t_ids], "Cyanobiaceae"),
  mk_l2(fam_ra["Vibrionaceae", z_ids], fam_ra["Vibrionaceae", t_ids], "Vibrionaceae"),
  mk_l2(otu_ra["ASV_8", z_ids], otu_ra["ASV_8", t_ids], "ASV_8", sps = "DV")
)
write.csv(f4b_df, file.path(DIR_INTERIM, "fig4B_within_pair_log2.csv"), row.names = FALSE)
cat("Saved Table S8 and Figure 4A/4B intermediate CSVs.\n")
print(sessionInfo())
