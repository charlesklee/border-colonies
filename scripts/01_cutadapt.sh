#!/usr/bin/env bash
# ==============================================================================
# 01_cutadapt.sh — Primer Trimming of 515F-Y / 926R 16S V4–V5 Paired-End Reads
# ==============================================================================
# Purpose : Trim 515F-Y (19 nt) and 926R (20 nt) primers using linked adapters
#           (anchored 5' required + optional 3' read-through) with cutadapt 5.2
#           (Python 3.11.2). Exact reported DADA2 counts require cutadapt 5.2
#           (cutadapt 4.2 differs by 1 edge read out of 61,587 in S-163).
# Inputs  : data/fastq_raw/*_R1_001.fastq(.gz) & *_R2_001.fastq(.gz) (182 pairs)
# Outputs : data/fastq_trimmed/*_R1_001.fastq.gz, *_R2_001.fastq.gz, logs
# Supports: Methods §2.4 (Bioinformatic preprocessing)
# ==============================================================================
set -euo pipefail

RAW_DIR="${1:-data/fastq_raw}"
OUT_DIR="${2:-data/fastq_trimmed}"
LOG_DIR="output/intermediate/cutadapt_logs"
THREADS="${THREADS:-4}"
CUTADAPT_BIN="${CUTADAPT_BIN:-cutadapt}"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

# Linked adapters (515F-Y / 926R):
#   515F-Y (19 nt): GTGYCAGCMGCCGCGGTAA ... rc(926R) (20 nt): AAACTYAAAKRAATTGRCGG
#   926R   (20 nt): CCGYCAATTYMTTTRAGTTT ... rc(515F-Y) (19 nt): TTACCGCGGCKGCTGRCAC
FWD_ADAPTER="^GTGYCAGCMGCCGCGGTAA...AAACTYAAAKRAATTGRCGG;optional"
REV_ADAPTER="^CCGYCAATTYMTTTRAGTTT...TTACCGCGGCKGCTGRCAC;optional"

echo "Running $(${CUTADAPT_BIN} --version) with ${THREADS} threads..."
shopt -s nullglob
R1_FILES=("${RAW_DIR}"/*_R1_001.fastq "${RAW_DIR}"/*_R1_001.fastq.gz)
if [[ ${#R1_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No *_R1_001.fastq(.gz) files found in ${RAW_DIR}"
  exit 1
fi

for r1 in "${R1_FILES[@]}"; do
  base="$(basename "${r1}")"
  r2="${r1/_R1_001.fastq/_R2_001.fastq}"
  out_r1="${OUT_DIR}/${base%.gz}.gz"
  out_r2="${OUT_DIR}/$(basename "${r2%.gz}").gz"
  sample_id="${base%%_L001_R1_001.fastq*}"

  "${CUTADAPT_BIN}" \
    -j "${THREADS}" \
    -g "${FWD_ADAPTER}" \
    -G "${REV_ADAPTER}" \
    --discard-untrimmed \
    --minimum-length 50 \
    -e 0.15 \
    -o "${out_r1}" \
    -p "${out_r2}" \
    "${r1}" "${r2}" \
    > "${LOG_DIR}/${sample_id}.cutadapt.log"
done

echo "Cutadapt trimming complete for ${#R1_FILES[@]} libraries."
