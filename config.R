# ==============================================================================
# config.R — Central Configuration for Colonial Ascidian Microbiome Pipeline
# ==============================================================================
# All paths are relative to the repository root (submission_code/).
# Modify ONLY this file to adjust threads, NCBI email placeholder, or paths.
# ==============================================================================

# Threading & Reproducibility
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
SEED        <- 42L
N_CORES     <- min(4L, max(1L, parallel::detectCores(logical = TRUE)))
N_CL_DA     <- min(3L, N_CORES)  # ANCOM-BC2 / MaAsLin2 worker cap

# Permutation & Monte Carlo Counts (Methods §2.6 — do not alter)
N_PERM          <- 9999L   # PERMANOVA, betadisper, within-species co-occurrence
N_SWAP_NULL     <- 999L    # Figure 2A / Table 1 within-pair Z/T label-swap null
N_ALDEX_MC      <- 128L    # Paired ALDEx2 Monte Carlo Dirichlet instances
NMDS_TRYMAX     <- 500L    # Figure 2B Bray-Curtis NMDS trymax
NMDS_MAXIT      <- 999L    # Figure 2B Bray-Curtis NMDS maxit

# Bioinformatic & QC Thresholds (Methods §2.4–§2.5)
MIN_BOOT_CANON       <- 50L     # Canonical SILVA v138.2 bootstrap threshold
MIN_BOOT_LOW         <- 30L     # Secondary lower-threshold label (Tables S2 & S7A)
DECONTAM_PREV_THRESH <- 0.50    # Prevalence threshold, batched by library plate
RESCUE_MIN_RA        <- 0.01    # >1% RA in >=1 tissue sample
RESCUE_FREQ_P        <- 0.10    # Frequency-method p > 0.10
RESCUE_RATIO         <- 50.0    # >=50x biological:control read ratio
MIN_READS_SAMPLE     <- 1000L   # Minimum post-decontam reads per biological sample
MIN_ASV_SAMPLE_GT    <- 2L      # Retain ASVs with >2 reads in >=1 biological sample (sum(x > 2) >= 1)
PRESENCE_MIN_COUNTS  <- 2L      # Presence threshold for categories, prevalence, co-occurrence
AITCHISON_MIN_PREV   <- 0.05    # >=5% sample prevalence for CLR / Aitchison sensitivity
DA_MIN_PREV          <- 0.10    # >=10% sample prevalence for differential abundance
COOCCUR_MIN_PAIRS    <- 5L      # >=5 pairs for within-specimen co-occurrence tests
ALPHA_SIG            <- 0.05    # Significance threshold (FDR / Holm / two-sided)

# Relative Directory Paths
DIR_DATA       <- "data"
DIR_META       <- file.path(DIR_DATA, "metadata")
DIR_REF        <- file.path(DIR_DATA, "reference")
DIR_FASTQ_RAW  <- file.path(DIR_DATA, "fastq_raw")
DIR_FASTQ_TRIM <- file.path(DIR_DATA, "fastq_trimmed")
DIR_OUT        <- "output"
DIR_INTERIM    <- file.path(DIR_OUT, "intermediate")
DIR_TABLES     <- file.path(DIR_OUT, "tables")
DIR_FIGURES    <- file.path(DIR_OUT, "figures")

for (d in c(DIR_OUT, DIR_INTERIM, DIR_TABLES, DIR_FIGURES)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Key File Paths
PATH_SAMPLE_META   <- file.path(DIR_META, "sample_metadata.csv")
PATH_ENA_RUN_MAP   <- file.path(DIR_META, "ena_run_map.csv")
PATH_ASV_CROSSWALK <- file.path(DIR_META, "asv_crosswalk.csv.gz")
PATH_SILVA_FASTA   <- file.path(DIR_REF,  "silva_nr99_v138.2_toSpecies_trainset.fa.gz")
PATH_REFSEQ_BAC    <- file.path(DIR_REF,  "bacteria.16SrRNA.fna.gz")
PATH_REFSEQ_ARC    <- file.path(DIR_REF,  "archaea.16SrRNA.fna.gz")
PATH_COTS27_FASTA  <- file.path(DIR_REF,  "cots27_LC490103.fasta")

# Canonical Intermediate Objects
PATH_SEQTAB_NOCHIM <- file.path(DIR_INTERIM, "seqtab_nochim_stageB.rds")
PATH_TAX_B50       <- file.path(DIR_INTERIM, "tax_silva138_2_minBoot50.rds")
PATH_TAX_B30       <- file.path(DIR_INTERIM, "tax_silva138_2_minBoot30.rds")
PATH_DECONTAM_REV  <- file.path(DIR_INTERIM, "decontam_review.rds")
PATH_PHY_CANON     <- file.path(DIR_INTERIM, "ascidian_canonical.rds")
PATH_PHY_ALL       <- file.path(DIR_INTERIM, "ascidian_canonical_allsamples.rds")

# Placeholders for External APIs (fill before running BLAST / ENA download)
ENA_STUDY_ACCESSION <- "PRJEBXXXXXX"
NCBI_EMAIL          <- "your.email@institution.edu"

# Fixed Host Species Metadata & Colour-Blind-Safe Palettes (Figs 1, 2, 3, 4)
SP_LEVELS <- c("AK", "AP", "BL", "DI", "DV", "LN")

SP_FULL_LATIN <- c(
  "AK" = "Aplidium kottae",
  "AP" = "Aplidium phortax",
  "BL" = "Botrylloides leachii",
  "DI" = "Didemnum incanum",
  "DV" = "Didemnum vexillum",
  "LN" = "Lissoclinum notti"
)

SP_LEGEND_EXPR <- c(
  "AK" = expression(italic("Aplidium kottae")~"(AK)"),
  "AP" = expression(italic("Aplidium phortax")~"(AP)"),
  "BL" = expression(italic("Botrylloides leachii")~"(BL)"),
  "DI" = expression(italic("Didemnum incanum")~"(DI)"),
  "DV" = expression(italic("Didemnum vexillum")~"(DV)"),
  "LN" = expression(italic("Lissoclinum notti")~"(LN)")
)

# Fixed Species Palette (Figures 3, 4, S1, S3)
COLS_SPECIES <- c(
  "AK" = "#D55E00", # Vermilion (Invasive)
  "AP" = "#56B4E9", # Sky Blue (Native)
  "BL" = "#CC79A7", # Reddish Purple (Invasive)
  "DI" = "#009E73", # Bluish Green (Native)
  "DV" = "#0072B2", # Blue (Invasive)
  "LN" = "#E69F00"  # Orange (Native)
)

# Fixed Sample-Type Palette (Figure 2)
COLS_SAMPLE_TYPE <- c(
  "Zooid"    = "#D55E00",
  "Tunic"    = "#009E73",
  "Seawater" = "#0072B2"
)

# Fixed Site Marker Palette (Figure 1B)
COLS_SITES <- c(
  "Mount"  = "#1976d2",
  "Pier A" = "#388e3c",
  "Bridge" = "#d32f2f"
)
