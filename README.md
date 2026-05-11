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

## Core Framework of the Pipeline

![](framework.png)

To systematically detect HGT events driven by FMT, we developed HGTector, a computational framework that integrates metagenomic assembly, homology search, phylogenetic assignment, and functional annotation (as shown above). Briefly, raw metagenomic reads from donor and recipient (pre‑FMT and post‑FMT) samples were subjected to de-contamination and de novo metagenomic assembly, followed by rigorous quality control and decontamination to remove potential contaminants, short contigs (＜5000 bp) or host-derived reads. By using thresholds of (1) ≥99.0% identity rate and (2) ≥90% length coverage, all qualified post-FMT recipient contigs were categorized into four categories: donor‑only contigs, recipient‑only contigs, shared contigs, and suspected HGT contigs.

For inferring HGT from donor microbes to recipient microbes, all suspected HGT contigs were subsequently screened to search for certain genomic regions (≥ 500 bp, ≤ 50% of full contig length) that are not exist in the corresponding contig of pre-FMT recipient sample but share high homology rate (≥99.0%, e-value ≤10-10) with the certain region of contigs in donor sample. Additionally, all suspected homologous regions must not overlap with any region within the contig from pre‑FMT recipient, with flanking and non-HGT regions sufficiently aligned to pre-FMT contig. Only contigs satisfying all criteria were retained as candidate harboring HGT events.

![](framework1.png)

The flanking loci of each HGT regions were further checked to ensure that the homologous segments were not false-positives due to insufficient coverage of sequencing or genome assembly. To further confirm the factuality of those detected HGT regions, the recipient- and donor-source context were extracted and independently aligned against standard reference microbial genome database. Only those cases with HGT regions and non-HGT loci phylogenetically assigned to divergent species were retained to implement final round checking, with filtered cases as bona fide HGT events. Simultaneously, we implemented similar approaches for inferring HGT from recipient microbes to the engrafted donor microbes.

## Usage and Examples
