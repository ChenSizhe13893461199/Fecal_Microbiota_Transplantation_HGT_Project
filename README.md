# Fecal_Microbiota_Transplantation_HGT_Project

This pipeline identifies putative horizontal gene transfer (HGT) events from longitudinal metagenomic samples (pre‑FMT, donor, post‑FMT). It combines stringent sequence alignment filtering, taxonomic assignment, gene prediction, functional annotation, and curation steps that remove false positives potentially caused by assembly edge effects and other possible factors.

## Requirements
### System Dependencies
Python 3.9+ with packages: pandas, openpyxl, biopython

BLAST+ (makeblastdb, blastn)

Prodigal (gene prediction)

eggNOG‑mapper (emapper.py) and its diamond dependency

Kraken2 and Taxonkit (species assignment)

seqtk (optional, for faster FASTA extraction)

Standard Unix tools: awk, sed, sort, paste, mkdir

Conda Environments (recommended)
bash
# Create environment for main pipeline
conda create -n hgt_pipeline python=3.9 pandas openpyxl biopython blast prodigal
conda activate hgt_pipeline

# For eggnog‑mapper (separate environment recommended)
conda create -n eggnog -c bioconda eggnog-mapper
conda activate eggnog

# For Kraken2 (separate environment)
conda create -n kraken2 -c bioconda kraken2 taxonkit
conda activate kraken2
