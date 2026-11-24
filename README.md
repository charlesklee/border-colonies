# border-colonies

**Microbiome analysis pipeline for colonial ascidian zooid–tunic compartments and interspecific interaction borders in Tauranga Harbour, New Zealand.**

> **Manuscript**: [Draft Title] Tissue-specific microbiome in colonial ascidians: are zooids more than filter-feeding vessels?
> **Authors**: Shir Livne Fadlon, Mafalda S. Baptista, Michele Prinsep, Christopher N. Battershill, Charles K. Lee
> **Citation**: *[Journal, Year, DOI Placeholder]*  
> **Code Archive (Zenodo)**: `https://doi.org/10.5281/zenodo.XXXXXXX` *(placeholder)*  
> **Raw Sequencing Data (ENA)**: Study Accession `PRJEBXXXXXX` *(placeholder)*  
> **License**: Code under [MIT License](LICENSE); metadata, tables, and documentation under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

---

## 1. Overview & Data Availability

`border-colonies` is a self-contained, end-to-end R / Bash / Python pipeline that reproduces every bioinformatic processing step, statistical test, figure, and supplementary table in the manuscript directly from the **182 paired-end 16S rRNA (V4–V5, 515F-Y / 926R) MiSeq libraries** deposited at the European Nucleotide Archive (`PRJEBXXXXXX`):

- **167 biological libraries**: `164` ascidian tissue samples across **6 colonial ascidian species** — native (*Aplidium phortax* [`AP`], *Didemnum incanum* [`DI`], *Lissoclinum notti* [`LN`]) and invasive (*Aplidium kottae* [`AK`], *Botrylloides leachii* [`BL`], *Didemnum vexillum* [`DV`]) — from `40` field specimens (`19` monospecific colonies × 2 tissues + `21` interaction encounters × 6 tissues, including interaction borders), plus `3` site-matched ambient seawater samples (`Bridge`, `Pier A`, `Mount`).
- **15 control libraries**: `13` PCR no-template controls (`PCR_NTC`), `1` extraction blank (`S-164`), and `1` artificial seawater procedural control (`S-163`, flagged as `ASW_process_control` in `data/metadata/ena_run_map.csv` and excluded from the `decontam` negative-control pool).
- **Non-deposited run files**: No `Undetermined`, PhiX, or index-only FASTQ files were present in the raw sequencing run directory; run-level control files are not deposited.

All paths, random seeds (`set.seed(42)`), permutation counts, and QC thresholds are centralized in [`config.R`](config.R).

---

## 2. Repository Layout

```text
border-colonies/
├── README.md
├── LICENSE
├── config.R                          # Centralized paths, seeds, thresholds, palettes
├── sessionInfo.txt                   # Exact R 4.2.2 / Bioconductor package versions
├── INVENTORY.md                      # Provenance mapping from development blocks to submission scripts
├── data/
│   ├── metadata/
│   │   ├── sample_metadata.csv       # Full metadata for all 182 sequenced libraries (incl. Plate/Well, DNA ng/uL)
│   │   ├── ena_run_map.csv           # Sample ID <-> ENA run accession map (182 libraries)
│   │   └── asv_crosswalk.csv.gz      # Canonical ASV_n IDs, DNA sequences, minBoot 50 & minBoot 30 taxonomy
│   └── reference/
│       └── README.md                 # Download instructions for SILVA 138.2, RefSeq 16S, and COTS27
├── docs/
│   └── figS1_prompt.md               # Prompt used to generate Figure S1 artwork from audited counts
├── scripts/
│   ├── 00_download_ena.sh            # Download 182 paired-end FASTQ libraries from ENA
│   ├── 01_cutadapt.sh                # Primer trimming (515F-Y / 926R, cutadapt 5.2)
│   ├── 02_dada2.R                    # DADA2 denoising, merging, and chimera removal
│   ├── 03_taxonomy.R                 # SILVA v138.2 assignment (minBoot 50 canonical + minBoot 30 secondary)
│   ├── 04_decontam_and_filters.R     # Plate-batched prevalence decontam + biological rescue rule
│   ├── 05_build_canonical_object.R   # Organelle & depth filters, canonical ASV_n numbering
│   ├── 10_sample_type_permanova.R    # Sample-type PERMANOVA & PERMDISP (Tables S3A, S3C)
│   ├── 11_categories_and_null.R      # Zooid/Tunic/Seawater categories, 999-swap null, NMDS (Table 1, Fig 2, S6A)
│   ├── 12_species_status_IB.R        # Host species, paired tissue, varpart, status & IB (Fig 3, S3, S4, S3B, S4A–C)
│   ├── 13_differential_abundance.R   # Paired ANCOM-BC2, ALDEx2, MaAsLin2 (Fig 4, Table S8)
│   ├── 14_dominance_prevalence_cooccurrence.R # Dominant ASVs, prevalence, CMH/perm co-occurrence (S6B–C, S7A–C)
│   ├── 15_seawater_and_carryover.R   # Control attribution, seawater ASV_1:2:3 ratios, bacterioplankton (S1C–D, S5A–B)
│   ├── 16_blast_local16s.R           # Local Biostrings alignment vs RefSeq 16S type strains & COTS27 (Table S2)
│   ├── 16_blast.py                   # NCBI REST megablast vs core_nt + BioSample isolation source (Table S2)
│   ├── 18_sample_flow_audit.R        # Metadata-to-phyloseq sample-flow audit (Table S1A, Figure S1 counts)
│   ├── 20_figures.R                  # Publication figures (Figures 1A–B, 2, 3, 4, S3, S4)
│   └── 21_supplementary_tables.R     # Consolidates Tables S1–S8 into output/tables/
└── output/                           # Created at runtime (intermediate/, tables/, figures/)
```

---

## 3. Software Requirements & Installation

Environment used for all reported manuscript results (Debian GNU/Linux 12, `x86_64`; full dump in [`sessionInfo.txt`](sessionInfo.txt)):

- **R**: `4.2.2 Patched (2022-11-10 r83330)`
- **Core R / Bioconductor packages**: `dada2 1.26.0`, `decontam 1.18.0`, `phyloseq 1.42.0`, `Biostrings 2.66.0`, `microbiome 1.20.0`, `vegan 2.7-2`, `permute 0.9-10`, `quadprog 1.5-8`, `ANCOMBC 2.0.2`, `ALDEx2 1.30.0`, `Maaslin2 1.12.0`, `ggplot2 4.0.2`, `patchwork 1.3.2`, `ragg 1.2.5`, `dplyr 1.2.0`, `tidyr 1.3.2`, `sf`, `maps`, `mapdata`, `jsonlite 2.0.0`.
- **CLI / Python**: `cutadapt 5.2` (Python 3.11.2; exact version used for reported DADA2 counts — `cutadapt 4.2` differs by 1 edge read in `S-163`), `Python 3.11+` (standard library only for `16_blast.py`), optional `BLAST+ 2.x`.
- **Hardware**: `8–16 GB RAM`, `<= 4` CPU threads (`N_CORES <- 4L` in `config.R`).

### Reference Databases

Before running `scripts/03_taxonomy.R` or `scripts/16_blast_local16s.R`, download the external reference databases following [`data/reference/README.md`](data/reference/README.md):
1. **SILVA v138.2** DADA2 training set (`silva_nr99_v138.2_toSpecies_trainset.fa.gz`).
2. **NCBI RefSeq Targeted Loci 16S rRNA** type-strain FASTAs (`bacteria.16SrRNA.fna.gz`, `archaea.16SrRNA.fna.gz`).
3. **COTS27** 16S rRNA sequence (`LC490103.1`, Wada et al. 2020).

---

## 4. Quick Start & Pipeline Run Order

Run all scripts from the repository root (`border-colonies/`):

```bash
# 1. Download raw reads from ENA & trim primers
bash scripts/00_download_ena.sh PRJEBXXXXXX
bash scripts/01_cutadapt.sh

# 2. DADA2 denoising, taxonomy, decontamination & canonical phyloseq objects
Rscript scripts/02_dada2.R
Rscript scripts/03_taxonomy.R
Rscript scripts/04_decontam_and_filters.R
Rscript scripts/05_build_canonical_object.R

# 3. Statistical analyses, BLAST & sample-flow audit
Rscript scripts/10_sample_type_permanova.R
Rscript scripts/11_categories_and_null.R
Rscript scripts/12_species_status_IB.R
Rscript scripts/13_differential_abundance.R
Rscript scripts/14_dominance_prevalence_cooccurrence.R
Rscript scripts/15_seawater_and_carryover.R
Rscript scripts/16_blast_local16s.R
NCBI_EMAIL="your.email@institution.edu" python3 scripts/16_blast.py
Rscript scripts/18_sample_flow_audit.R

# 4. Generate publication figures & consolidate Tables S1–S8
Rscript scripts/20_figures.R
Rscript scripts/21_supplementary_tables.R
```

### Runtime & Manuscript Mapping

| Order | Script | Approx. Runtime | Manuscript Figures / Tables Produced |
| :---: | :--- | :--- | :--- |
| 0 | `scripts/00_download_ena.sh` | 10–20 min | Downloads all 182 paired-end FASTQ libraries (`data/metadata/ena_run_map.csv`) |
| 1 | `scripts/01_cutadapt.sh` | ~5 min | §2.4 primer trimming (`515F-Y` / `926R`, `cutadapt 5.2`) |
| 2 | `scripts/02_dada2.R` | ~45–90 min | §2.4 DADA2 denoising, merging, chimera removal (`seqtab_nochim_stageB.rds`) |
| 3 | `scripts/03_taxonomy.R` | ~15 min | §2.4 SILVA v138.2 (`minBoot = 50` canonical; `minBoot = 30` secondary column) |
| 4 | `scripts/04_decontam_and_filters.R` | ~1 min | §2.5 plate-batched `decontam` (prev 0.5), 28 rescued (20 post-organelle), 92 removed; **Table S1B** |
| 5 | `scripts/05_build_canonical_object.R` | <1 min | `ascidian_canonical.rds` (`151 × 18,131`), `ascidian_canonical_allsamples.rds` (`167 × 18,174`), `ASV_n` numbering |
| 6 | `scripts/10_sample_type_permanova.R` | ~5 min | §3.1; **Table S3A** (sample-type PERMANOVA: 107, 151, 167), **Table S3C** (PERMDISP) |
| 7 | `scripts/11_categories_and_null.R` | ~2 min | §3.2; **Table 1**, **Figure 2** data (999-swap null, NMDS), **Table S6A** |
| 8 | `scripts/12_species_status_IB.R` | ~15 min | §3.3; **Figure 3** data, **Figures S3–S4** data, **Table S3B**, **Table S4A–C** |
| 9 | `scripts/13_differential_abundance.R` | ~15 min | §3.4; **Figure 4** data, **Table S8** (ANCOM-BC2, paired ALDEx2, MaAsLin2) |
| 10 | `scripts/14_dominance_prevalence_cooccurrence.R` | ~2 min | §3.5–§3.6; **Table S6B–C**, **Table S7A–C** |
| 11 | `scripts/15_seawater_and_carryover.R` | ~2 min | §2.4, §3.2; **Table S1C–D**, **Table S5A–B** |
| 12 | `scripts/16_blast_local16s.R` & `scripts/16_blast.py` | ~10–30 min | §3.5, §4; **Table S2** (`31_ASV_FULL_BLAST_TABLE_v6.csv`, `blast_summary_v6.csv`) |
| 13 | `scripts/18_sample_flow_audit.R` | <10 s | **Table S1A**; **Figure S1** counts (`block18_counts.md`) |
| 14 | `scripts/20_figures.R` | ~1 min | **Figures 1A–B, 2, 3, 4, S3, S4** |
| 15 | `scripts/21_supplementary_tables.R` | <10 s | Consolidates **Tables S1–S8** into `output/tables/` |

### Manuscript Items Not Produced by Code

- **Figure 1C**: Google Earth satellite imagery of Tauranga Harbour sampling sites.
- **Figure 1D**: Field photographs of colonial ascidians.
- **Figure S1**: Counts produced by `scripts/18_sample_flow_audit.R` (`block18_counts.md`); final artwork generated with Gemini 4 Neon from these counts (see [`docs/figS1_prompt.md`](docs/figS1_prompt.md), [`docs/figS1_sample_flow.svg`](docs/figS1_sample_flow.svg), and manuscript LLM statement).
- **Figure S2**: Dissection schematics adapted from Shir Livne Fadlon's 2022 MSc thesis.

---

## 5. Permutation Counts & Random Seeds (Methods §2.6)

| Analysis | Script | Permutations / Iterations | Seed |
| :--- | :--- | :---: | :---: |
| Sample-type PERMANOVA & PERMDISP (Tables S3A, S3C) | `10_sample_type_permanova.R` | `9,999` | `set.seed(42)` |
| Within-pair Z/T label-swap null (Figure 2A, Table 1) | `11_categories_and_null.R` | `999` | `set.seed(42)` |
| Bray–Curtis NMDS (Figure 2B) | `11_categories_and_null.R` | `trymax = 500, maxit = 999` | `set.seed(42)` |
| Species, paired Tissue, Pairwise Species, VarPart, Within-event, Status PERMANOVA (Tables S3B, S4A–C) | `12_species_status_IB.R` | `9,999` | `set.seed(42)` |
| Paired ALDEx2 Dirichlet Monte Carlo (Table S8) | `13_differential_abundance.R` | `128` MC instances | `set.seed(42)` |
| Within-specimen co-occurrence permutation (Table S7C) | `14_dominance_prevalence_cooccurrence.R` | `9,999` | `set.seed(42)` |

---

## 6. BLAST Notes (`Table S2`)

- **Reported results**: `core_nt` queries were executed via the official NCBI BLAST REST API (`scripts/16_blast.py`, megablast, `HITLIST_SIZE = 50`, queried on **2026-11-23**); because NCBI `core_nt` is updated continuously, re-running at a later date may return newly deposited accessions. RefSeq 16S type-strain (`28,048` bacterial + archaeal sequences) and COTS27 (`LC490103.1`) alignments were executed locally via `Biostrings 2.66.0` (`scripts/16_blast_local16s.R`, exact BLASTn `+2/-3`, gap `5/2` scoring).
- **Optional equivalent BLAST+ 2.x CLI commands (not used for reported results)**:
  ```bash
  blastn -task megablast -db core_nt -query block19_targets.fasta -max_target_seqs 50 -outfmt 6 -num_threads 4
  blastn -task blastn    -db 16S_ribosomal_RNA -query block19_targets.fasta -max_target_seqs 10 -outfmt 6 -num_threads 4
  blastn -task blastn    -query asv8_17.fasta -subject data/reference/cots27_LC490103.fasta -outfmt 6
  ```

---

## 7. License

- **Code** (`scripts/`, `config.R`): [MIT License](LICENSE) — permissive reuse with attribution.
- **Data, metadata, tables & documentation** (`data/`, `docs/`, `README.md`): [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

## 8. Citation

If you use this pipeline or dataset, please cite the accompanying manuscript and the archived Zenodo release:

```bibtex
@article{author2026bordercolonies,
  title   = {[Manuscript Title Placeholder]},
  author  = {[Author List Placeholder]},
  journal = {[Journal Placeholder]},
  year    = {2026},
  doi     = {[Manuscript DOI Placeholder]}
}
```
