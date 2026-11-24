# Reference Databases (`data/reference/`)

Large external reference databases are not bundled in this repository. Download them into `data/reference/` before running `scripts/03_taxonomy.R` and `scripts/16_blast.py` / `scripts/16_blast_local16s.R`.

## 1. SILVA v138.2 16S rRNA Training Set (DADA2 format)
Used by `scripts/03_taxonomy.R` (`minBoot = 50` canonical; `minBoot = 30` secondary column for Tables S2 & S7A):
```bash
# Zenodo record for SILVA v138.2 DADA2 training sets (McLaren / Callahan)
curl -L -o data/reference/silva_nr99_v138.2_toSpecies_trainset.fa.gz \
  "https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toSpecies_trainset.fa.gz?download=1"
```

## 2. NCBI RefSeq Targeted Loci 16S rRNA Type-Strain FASTAs
Used by `scripts/16_blast_local16s.R` (local Smith–Waterman BLASTn alignment via `Biostrings`, producing the `16S_type` hits in Table S2):
```bash
curl -fsSL -o data/reference/bacteria.16SrRNA.fna.gz \
  "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"

curl -fsSL -o data/reference/archaea.16SrRNA.fna.gz \
  "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz"
```

## 3. COTS27 Reference 16S rRNA Sequence (`LC490103.1`, Wada et al. 2020)
Used by `scripts/16_blast_local16s.R` to confirm pairwise identity of `ASV_8` and `ASV_17` against the crown-of-thorns starfish spirochete symbiont COTS27:
```bash
curl -fsSL -o data/reference/cots27_LC490103.fasta \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC490103.1&rettype=fasta&retmode=text"
```
