# Fecal_Microbiota_Transplantation_HGT_Project

This pipeline identifies putative horizontal gene transfer (HGT) events from longitudinal metagenomic samples (pre‑FMT, donor, post‑FMT). It combines sequence alignment filtering, taxonomic assignment, gene prediction, functional annotation, and curation steps that remove false positives potentially caused by assembly edge effects and other possible factors.

## Requirements
### System Dependencies
- Python 3.9+ with packages: pandas, openpyxl, biopython

- BLAST+ (makeblastdb, blastn)

- Prodigal (gene prediction)

- eggNOG‑mapper (emapper.py) and its diamond dependency

- Kraken2 and Taxonkit (species assignment)

- seqtk (optional, for faster FASTA extraction)

- Standard Unix tools: awk, sed, sort, paste, mkdir

## Core Framework of the Pipeline

![](framework.png)

- To systematically identify HGT events driven by FMT, we developed HGTector, a computational framework that integrates metagenomic assembly, homology search, phylogenetic assignment, and functional annotation (as shown above). 

- Briefly, raw metagenomic reads from donor and recipient (pre‑FMT and post‑FMT) samples were subjected to de-contamination and de novo metagenomic assembly, followed by rigorous quality control and decontamination to remove potential contaminants.
![](framework1.png)

### Notes: 
- Please appropriately assign the corresponding path of directory containing assembled FMT contig files (e.g. pre-FMT recipient, post-FMT recipient, and Donor contig files), FMT metadata reference, taxonomic database, and gene annotation reference database in the corresponding file (e.g. HGT_main_implementing.sh, qc_annotation.py, convert.sh, and etc.). Full details are available in annotations of each computational script;

- For user who attempt to process long-reads sequencing data (e.g. PacBio type) by HGTector, please revise "./root.sh" in line 253 of HGT_main_implementing.sh as "./root1.sh";

- Please feel free to contact us via Chen2422679942@163.com and we are willing to provide necessary assistance for implementing HGTector.

### A Prelimnary Step-by-step Example
- In your initial path, you have 3 different directories named as __"FMT_HGT"__, __"tax"__, and __"reference_data_base"__, as shown below:

![](1.png)

Metagenomic Assembled Files (pre-FMT recipient, post-FMT recipient, and Donor) are deposited in __"FMT_HGT"__, with phylogenetic annotation files contained in __"tax"__ and __"reference_data_base"__, respectively. In __"FMT_HGT"__, we have 3 metagenomic assembled files (as shown below), with contigs < 5000 bp filtered.

![](2.png)

