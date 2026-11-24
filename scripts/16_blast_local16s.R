#!/usr/bin/env Rscript
# ==============================================================================
# 16_blast_local16s.R — Local RefSeq 16S Type-Strain & COTS27 Alignment (Table S2)
# ==============================================================================
# Purpose : Align the 39 target ASVs against official NCBI RefSeq Targeted Loci
#           16S rRNA bacterial + archaeal type strains (data/reference/) and
#           align ASV_8 and ASV_17 against COTS27 (LC490103.1) using Biostrings
#           local Smith-Waterman alignment with exact BLASTn scoring (+2/-3,
#           gap open 5, gap extend 2). Updates block19_blast_cache.json for
#           scripts/16_blast.py.
# Inputs  : config.R, output/intermediate/block19_target_ASVs.csv,
#           data/reference/bacteria.16SrRNA.fna.gz,
#           data/reference/archaea.16SrRNA.fna.gz,
#           data/reference/cots27_LC490103.fasta
# Outputs : output/intermediate/block19_blast_cache.json (16S_type & cots27)
# Supports: Table S2; §3.5 / §4 (ASV_8, ASV_17, COTS27 confirmation)
# ==============================================================================

source("config.R")
set.seed(SEED)

suppressPackageStartupMessages({
  library(phyloseq)
  library(Biostrings)
  library(jsonlite)
  library(parallel)
})

targets_csv <- file.path(DIR_INTERIM, "block19_target_ASVs.csv")
cache_file  <- file.path(DIR_INTERIM, "block19_blast_cache.json")

# Auto-generate block19_target_ASVs.csv if absent so the pipeline runs end-to-end
if (!file.exists(targets_csv)) {
  cat("Generating block19_target_ASVs.csv from canonical phyloseq, crosswalk, Table S7A & Table S8...\n")
  phy_can <- readRDS(PATH_PHY_CANON)
  cw      <- read.csv(PATH_ASV_CROSSWALK, stringsAsFactors = FALSE)
  dom_tb  <- read.csv(file.path(DIR_INTERIM, "table_S7A_dominant.csv"), stringsAsFactors = FALSE)
  da_tb   <- read.csv(file.path(DIR_INTERIM, "table_S8_differential_abundance.csv"), stringsAsFactors = FALSE)

  sd <- as(sample_data(phy_can), "data.frame"); sd$SequencingID <- rownames(sd)
  sd$IB <- as.character(sd$Interactive_Border) == "Y"
  sd$pair_id <- ifelse(sd$Zooid_or_Tunic == "Seawater", as.character(sd$Original_Identifier),
                       sub("[ZT]$", "", as.character(sd$Original_Identifier)))
  sd_nonib <- sd[sd$Zooid_or_Tunic %in% c("Z", "T") & !sd$IB & sd$Species %in% SP_LEVELS, ]
  ptab <- table(sd_nonib$pair_id, sd_nonib$Zooid_or_Tunic)
  comp52 <- sort(rownames(ptab)[ptab[, "Z"] == 1 & ptab[, "T"] == 1])
  z_ids <- sapply(comp52, function(p) sd_nonib$SequencingID[sd_nonib$pair_id == p & sd_nonib$Zooid_or_Tunic == "Z"])
  t_ids <- sapply(comp52, function(p) sd_nonib$SequencingID[sd_nonib$pair_id == p & sd_nonib$Zooid_or_Tunic == "T"])
  sids_104 <- c(z_ids, t_ids)

  otu_all <- as(otu_table(phy_can), "matrix"); if (!taxa_are_rows(phy_can)) otu_all <- t(otu_all)
  otu_104 <- otu_all[, sids_104, drop = FALSE]
  ra_Z    <- sweep(otu_104[, z_ids, drop = FALSE], 2, colSums(otu_104[, z_ids, drop = FALSE]), "/")
  tax_b50 <- as.data.frame(tax_table(phy_can), stringsAsFactors = FALSE)

  set_1a <- c("ASV_1", "ASV_2", "ASV_3", "ASV_4", "ASV_8", "ASV_17",
              "ASV_21", "ASV_22", "ASV_24", "ASV_270", "ASV_278", "ASV_414", "ASV_565")
  set_1b <- unique(dom_tb$ASV)
  da_asvs <- unique(da_tb$Taxon_ID[da_tb$Level == "ASV"])
  syn_ids <- rownames(tax_b50)[!is.na(tax_b50$Genus) & tax_b50$Genus == "Synechococcus CC9902" & rowMeans(ra_Z) > 0.001]
  vib_ids <- rownames(tax_b50)[!is.na(tax_b50$Family) & tax_b50$Family == "Vibrionaceae" & rowMeans(ra_Z) > 0.001]
  all_targets <- unique(c(set_1a, set_1b, da_asvs, syn_ids, vib_ids))
  all_targets <- all_targets[order(as.integer(sub("^ASV_", "", all_targets)))]

  fmt_silva_b30 <- function(t50, t30) {
    p50 <- trimws(unlist(strsplit(t50, ";", fixed = TRUE)))
    p30 <- trimws(unlist(strsplit(t30, ";", fixed = TRUE)))
    out <- character(6)
    for (j in 1:6) {
      v50 <- if (j <= length(p50) && p50[j] != "NA" && nzchar(p50[j])) p50[j] else "NA"
      v30 <- if (j <= length(p30) && p30[j] != "NA" && nzchar(p30[j])) p30[j] else "NA"
      if (v50 == "NA" && v30 != "NA") out[j] <- sprintf("%s (bootstrap <50)", v30)
      else if (v50 != v30 && v30 != "NA") out[j] <- sprintf("%s [%s (bootstrap <50)]", v50, v30)
      else out[j] <- v50
    }
    paste(out[2:6], collapse = "; ")
  }
  grp_104 <- paste0(as.character(sd[sids_104, "Species"]), "_", ifelse(sd[sids_104, "Zooid_or_Tunic"] == "Z", "Zooid", "Tunic"))
  t_rows <- lapply(all_targets, function(aid) {
    idx <- match(aid, cw$ASV_ID)
    v_cnt <- if (aid %in% rownames(otu_104)) otu_104[aid, ] else rep(0, length(sids_104))
    gs <- sort(tapply(v_cnt, grp_104, sum), decreasing = TRUE); gs <- gs[gs > 0]
    dom_s <- if (length(gs) >= 2) sprintf("%s (%s), %s (%s)", names(gs)[1], format(gs[1], big.mark=","), names(gs)[2], format(gs[2], big.mark=","))
             else if (length(gs) == 1) sprintf("%s (%s)", names(gs)[1], format(gs[1], big.mark=",")) else "IB/SW only"
    data.frame(
      ASV = aid, Inclusion_Reason = "Auto_generated",
      Old_IDs_from_31table = "", Old_IDs_crosswalk = cw$Original_ASV_IDs_Trimmed_Match[idx],
      Length_bp = cw$Length_bp[idx], Total_reads = cw$Reads_Canonical_151[idx],
      Dominant_species_tissue = dom_s,
      SILVA_138.2_taxonomy = fmt_silva_b30(cw$Taxonomy_minBoot50[idx], cw$best_assignment_minBoot30[idx]),
      Sequence = cw$ASV_DNA[idx], stringsAsFactors = FALSE
    )
  })
  write.csv(do.call(rbind, t_rows), targets_csv, row.names = FALSE)
}

stopifnot(file.exists(targets_csv), file.exists(PATH_REFSEQ_BAC), file.exists(PATH_REFSEQ_ARC), file.exists(PATH_COTS27_FASTA))

targets <- read.csv(targets_csv, stringsAsFactors = FALSE)
q_seqs  <- DNAStringSet(setNames(targets$Sequence, targets$ASV))
ref_all <- c(readDNAStringSet(PATH_REFSEQ_BAC), readDNAStringSet(PATH_REFSEQ_ARC))
db_len  <- sum(width(ref_all))
cots27  <- readDNAStringSet(PATH_COTS27_FASTA)

submat <- nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = FALSE)
LAMBDA <- 1.28; K_PAR <- 0.46

calc_hit <- function(q_str, s_str, s_header, db_size_bp) {
  qlen <- nchar(q_str)
  pa <- pairwiseAlignment(DNAString(q_str), DNAString(s_str), type = "local",
                          substitutionMatrix = submat, gapOpening = 5, gapExtension = 2)
  raw_sc <- score(pa)
  bit_sc <- (LAMBDA * raw_sc - log(K_PAR)) / log(2)
  eval_v <- K_PAR * qlen * db_size_bp * exp(-LAMBDA * raw_sc)
  p_ch   <- strsplit(as.character(alignedPattern(pa)), "", fixed = TRUE)[[1]]
  s_ch   <- strsplit(as.character(alignedSubject(pa)), "", fixed = TRUE)[[1]]
  al_len <- length(p_ch)
  pct_id <- if (al_len > 0) 100 * sum(p_ch == s_ch & p_ch != "-") / al_len else 0
  qcov   <- 100 * width(pattern(pa)) / qlen
  acc    <- sub("\\s.*$", "", s_header)
  desc   <- sub("^\\S+\\s+", "", s_header)
  sci    <- trimws(sub("\\s+(strain|16S|partial|complete|ribosomal).*$", "", desc, ignore.case = TRUE))
  list(Accession = acc, Description = desc, Scientific_name = sci,
       Pct_identity = round(pct_id, 2), Query_cover = round(qcov, 1),
       Align_len = as.integer(al_len), E_value = signif(eval_v, 3), Bit_score = round(bit_sc, 1))
}

k <- 7L
q_kmers <- lapply(as.character(q_seqs), function(s) unique(substring(s, 1:(nchar(s) - k + 1), k:nchar(s))))
names(q_kmers) <- names(q_seqs)
all_kmers <- unique(unlist(q_kmers)); all_kmers <- all_kmers[grepl("^[ACGT]+$", all_kmers)]
pdict <- PDict(DNAStringSet(all_kmers))
hits_by_ref <- vwhichPDict(pdict, ref_all)
k2i <- setNames(seq_along(all_kmers), all_kmers)
q_ki <- lapply(q_kmers, function(km) unname(na.omit(k2i[km])))
refs_per_k <- split(rep(seq_along(ref_all), lengths(hits_by_ref)), unlist(hits_by_ref, use.names = FALSE))

res_16s <- mclapply(names(q_seqs), function(aid) {
  cands <- unlist(refs_per_k[as.character(q_ki[[aid]])], use.names = FALSE)
  if (length(cands) == 0) return(list())
  tab <- tabulate(cands, nbins = length(ref_all))
  top40 <- order(tab, decreasing = TRUE)[1:40]; top40 <- top40[tab[top40] > 0]
  al <- lapply(top40, function(ri) calc_hit(as.character(q_seqs[[aid]]), as.character(ref_all[[ri]]), names(ref_all)[ri], db_len))
  al[order(sapply(al, function(x) -x$Bit_score))[1:min(3L, length(al))]]
}, mc.cores = N_CORES)
names(res_16s) <- names(q_seqs)

res_cots <- lapply(c("ASV_8", "ASV_17"), function(aid) {
  list(calc_hit(as.character(q_seqs[[aid]]), as.character(cots27[[1]]), names(cots27)[1], nchar(as.character(cots27[[1]]))))
})
names(res_cots) <- c("ASV_8", "ASV_17")

cache <- if (file.exists(cache_file)) fromJSON(cache_file, simplifyVector = FALSE) else list(core_nt = list(), `16S_type` = list(), nuccore_meta = list(), cots27 = list())
cache[["16S_type"]] <- res_16s
cache[["cots27"]]   <- res_cots
write_json(cache, cache_file, auto_unbox = TRUE, pretty = TRUE)
cat("Updated local RefSeq 16S and COTS27 alignments.\n")
print(sessionInfo())
