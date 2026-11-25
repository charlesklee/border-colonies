#!/usr/bin/env bash
# ==============================================================================
# 00_download_ena.sh — Download Raw Paired-End FASTQ Files from ENA
# ==============================================================================
# Purpose : Query the ENA Portal API for study accession PRJEB130768 (ERP209816)
#           and download all 182 paired-end libraries listed in
#           data/metadata/ena_run_map.csv into data/fastq_raw/.
# Inputs  : data/metadata/ena_run_map.csv
# Outputs : data/fastq_raw/*_R1_001.fastq.gz, data/fastq_raw/*_R2_001.fastq.gz
# Supports: Methods §2.3 (Data availability & raw sequencing reads)
# ==============================================================================
set -euo pipefail

STUDY_ACC="${1:-PRJEB130768}"
MAP_CSV="data/metadata/ena_run_map.csv"
OUT_DIR="data/fastq_raw"

mkdir -p "${OUT_DIR}"
echo "Fetching ENA filereport for study ${STUDY_ACC}..."
REPORT="${OUT_DIR}/ena_filereport_${STUDY_ACC}.tsv"
curl -fsSL \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${STUDY_ACC}&result=read_run&fields=run_accession,sample_alias,fastq_ftp,fastq_md5" \
  -o "${REPORT}"

echo "Downloading paired-end FASTQ files listed in ${REPORT}..."
tail -n +2 "${REPORT}" | while IFS=$'\t' read -r run_acc alias ftp_urls md5s; do
  IFS=';' read -ra URL_ARR <<< "${ftp_urls}"
  for u in "${URL_ARR[@]}"; do
    fname="$(basename "${u}")"
    if [[ ! -f "${OUT_DIR}/${fname}" ]]; then
      echo "  Downloading ${fname} (${run_acc})..."
      curl -fsSL "https://${u}" -o "${OUT_DIR}/${fname}"
    fi
  done
done

echo "Download complete. Verify downloaded filenames against ${MAP_CSV}."
