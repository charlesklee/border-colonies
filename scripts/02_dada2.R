#!/usr/bin/env Rscript
# ==============================================================================
# 02_dada2.R — DADA2 Quality Filtering, Denoising, Merging & Chimera Removal
# ==============================================================================
# Purpose : Run DADA2 1.26.0 on cutadapt 5.2 primer-trimmed paired-end FASTQs
#           (182 libraries: 167 biological + 13 NTC + 1 extraction blank + 1 ASW).
# Inputs  : config.R, data/fastq_trimmed/*_R1_001.fastq.gz & *_R2_001.fastq.gz
# Outputs : output/intermediate/seqtab_nochim_stageB.rds
#           output/intermediate/dada2_read_tracking.csv
# Supports: Methods §2.4 (DADA2 denoising, chimera removal, read counts)
# ==============================================================================

source("config.R")
set.seed(SEED)

suppressPackageStartupMessages({
  library(dada2)
})

cat("=== 02_dada2.R: DADA2 Denoising on Primer-Trimmed Reads ===\n")
fnFs <- sort(list.files(DIR_FASTQ_TRIM, pattern = "_R1_001\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(DIR_FASTQ_TRIM, pattern = "_R2_001\\.fastq\\.gz$", full.names = TRUE))
stopifnot(length(fnFs) > 0L, length(fnFs) == length(fnRs))

sample_names <- sub("_S[0-9]+_L001_R1_001\\.fastq\\.gz$", "", basename(fnFs))
filt_dir <- file.path(DIR_INTERIM, "dada2_filtered")
dir.create(filt_dir, showWarnings = FALSE, recursive = TRUE)

filtFs <- file.path(filt_dir, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sample_names, "_R_filt.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

# Truncation lengths adjusted for 515F-Y (19 nt) and 926R (20 nt) primer removal:
# c(235 - 19, 235 - 20) = c(216, 215)
out_filt <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen = c(216L, 215L),
  maxN = 0L,
  maxEE = c(2, 2),
  truncQ = 2L,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = N_CORES
)

keep_idx <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[keep_idx]
filtRs <- filtRs[keep_idx]

set.seed(SEED)
errF <- learnErrors(filtFs, multithread = N_CORES)
set.seed(SEED)
errR <- learnErrors(filtRs, multithread = N_CORES)

dadaFs  <- dada(filtFs, err = errF, multithread = N_CORES)
dadaRs  <- dada(filtRs, err = errR, multithread = N_CORES)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

seqtab <- makeSequenceTable(mergers)
seqtab_nochim <- removeBimeraDenovo(
  seqtab, method = "consensus", multithread = N_CORES, verbose = TRUE
)

saveRDS(seqtab_nochim, PATH_SEQTAB_NOCHIM)

getN <- function(x) sum(getUniques(x))
track <- data.frame(
  sample_id   = names(filtFs),
  trimmed_in  = out_filt[keep_idx, "reads.in"],
  filtered    = out_filt[keep_idx, "reads.out"],
  denoisedF   = sapply(dadaFs, getN),
  denoisedR   = sapply(dadaRs, getN),
  merged      = sapply(mergers, getN),
  nonchim     = rowSums(seqtab_nochim),
  stringsAsFactors = FALSE
)
write.csv(track, file.path(DIR_INTERIM, "dada2_read_tracking.csv"), row.names = FALSE)
cat(sprintf("Saved %s (%d samples x %d ASVs, %s nonchim reads).\n",
            PATH_SEQTAB_NOCHIM, nrow(seqtab_nochim), ncol(seqtab_nochim),
            format(sum(seqtab_nochim), big.mark = ",")))
print(sessionInfo())
