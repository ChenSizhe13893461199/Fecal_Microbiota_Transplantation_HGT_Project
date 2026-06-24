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

- For user who attempt to process long-reads sequencing data (e.g. PacBio type) by HGTector, please revise __"./root.sh"__ in line 253 of HGT_main_implementing.sh as "./root1.sh";

- Please feel free to contact us via Chen2422679942@163.com (__Dr. CHEN Sizhe__) and we are willing to provide necessary assistance for implementing HGTector.

### A Preliminary Step-by-step Guideline and Example

- The content below is an example for appropriately using HGTector. Assuming in your initial path, you have 3 different directories named as __"FMT_HGT"__, __"tax"__, and __"reference_data_base"__, as shown below:

![](workflow1.png)

- Metagenomic Assembled Files (pre-FMT recipient, post-FMT recipient, and Donor) are deposited in __"FMT_HGT"__, with phylogenetic annotation files contained in __"tax"__ and __"reference_data_base"__, respectively. In __"FMT_HGT"__, we have 3 metagenomic assembled files (as shown below), with contigs < 5000 bp filtered.

![](workflow2.png)

- In linux system, the necessary components in __"tax"__ can be downloaded via:

- wget -c https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
- tar -zxvf taxdump.tar.gz

- The necessary components in __"reference_data_base"__ can be prepared by following the guideline at https://benlangmead.github.io/aws-indexes/k2

- After all basic files prepared as aforementioned, we can start to use HGTector directly.

- First, please make sure that all scripts in this GitHub inventory have been appropriately put in your current user directory path, like the following:

![](workflow3.png)

- Here, the __"FMT_list.xlsx"__ contains FMT order information for HGTector to understand the longitudinal information of each FMT sample (__preFMT__, __postFMT__, and __Donor__ refers to prefix of __preFMT_filter.fasta__, __postFMT_filter.fast__, and __Donor_filter.fasta__, respectively). The content of it is similar to the format shown below:

![](workflow4.png)

- In next step, users can open the script of __HGT_main_implementing.sh__ for more details and information.

__HGT_main_implementing.sh__ is the script serves as the entry point and workflow manager for manipulating the entire HGTector.

- In this script, lines 1-6 shows the corresponding directory for __"FMT_HGT"__, __"FMT_list.xlsx"__, and intermediate processing directory path (as shown below). Please define the appropriate path for any user-customized conditions.

![](workflow5.png)

- __HGT_main_implementing.sh__ call modules or scripts of blastn, blastn_process.py, classifer.py, HGT.py, prodigal, emapper.py, qc_annotation.py, convert.sh, add_species_script.sh, root.sh, fill_species.sh, add_gc_nonhgt.py, filterchecking.py, quality_control_processing.py, etc. Each of the script has unique functions and full annotations of these scripts have been available in the corresponding script and the appropriate loaction in __HGT_main_implementing.sh__.

- Before utilization of HGTector, please pay attention to lines 170-173 in __HGT_main_implementing.sh__, and assign a user-customized parameter n to --min_recipient_cov. We have provided full information and reference regarding the selection of n in annotation part. In any kind of condition, n ≥ 0.5 (n < 1) should be maintained, with exact reasons available in the corresponding annottaion parts in scripts.

![](workflow6.png)

- Then, please open scripts of __convert.sh__, and __add_species_script.sh__ to indicate the exact directory path for the aforementioned __"tax"__, and __"reference_data_base"__, so that these scripts can smoothly function as wished.

- Lastly, implement HGTector by inputting commands of __"nohup ./HGT_main_implementing.sh &"__, and the HGTector pipeline will automatically start.

  #### Result and Explanation
- After implementing HGTector, users will obtain a .xlsx table and seperate .txt report for each FMT under directory path "HGT1_filtered". The .xlsx table looks like the following format:
  
![](workflow7.png)

The exact meanings of each column have been provided below: 

- Recipient_Base:
The base identifier of the recipient contig (from the post‑FMT sample) that acquired the HGT region. This column links the event to the specific contig in the recipient’s assembly.

- Donor_Base:
The base identifier of the donor contig (from the donor sample, e.g., D10H_1) that is the source of the transferred DNA. This column indicates which donor sequence contributed the HGT fragment.

- Pre_Recipient:
The pre‑FMT contig that best matches the recipient contig acquiring HGT region, extracted from the alignment between post‑FMT and pre‑FMT assemblies. This information is used in the judgment calculation to decide whether the HGT is a likely false positive.

- Judgment:
A binary value (0 or 1) that flags whether the HGT event is likely false-positive. In downstream parts of the HGTector, the flanking loci of each HGT regions were further checked to ensure that the homologous segments were not false-positives due to insufficient coverage of sequencing or genome assembly. Judgment = 1 means the event passes the edge‑based filter (only those with value of 1 are reserved in .xlsx table).

- Region_Length:
The length (in base pairs) of the HGT region on the recipient contig. It is calculated from the start and end coordinates of the transferred segment. This value describes the size of the potential horizontally acquired fragment.

- Homologous_Rate:
The percentage sequence identity between the recipient (HGT-acquiring) HGT region and its matched donor (HGT-delivery) contig region. It is usually very high (≥99.0%). This column quantifies the similarity of the transferred DNA.

- Gene_Count:
The number of protein‑coding genes predicted within the HGT region. It is aggregated from the eggNOG‑mapper annotations. This count reflects the functional complexity of the horizontally transferred segment.

- Gene_Description:
The functional annotation of all genes within the HGT region. This description comes from the eggNOG‑mapper output. It provides biological insight into the potential role of the transferred genetic material.

- Recipient_Species:
The taxonomic species assigned to the recipient contig (post‑FMT), based on the extracted contig sequence. This column indicates the species that harbours the HGT event in the post‑transplant sample.

- Donor_Species:
The taxonomic species assigned to the donor contig. It identifies the source organism from which the HGT fragment originated.

In addition, the .txt (report for each FMT) looks like the following format:

![](workflow8.png)

- The .txt result shown above indicate a high-confidence horizontal gene transfer event between two bacterial species: *Roseburia intestinalis* (HGt-acquiring) and *Anaerostipes hadrus* (HGT-delivery). 

- The entire HGT region spans approximately 31 kb (length 31010, __length__ column) and contains nine predicted genes, with HGT region sharing __99.291%__ identity (__rate__ column) (The column of __GC_HGT__ and __GC_nonHGT__ represent the GC content (percentage) of the HGT region and flanking non-HGT region on the HGT-acquiring contig.

- It differs by 44.95%-42.37%=2.58% GC content difference, suggesting that the HGT regions and background (non‑HGT) regions have distinct GC compositions). In addition, The recipient and donor species are phylogenetically different, satisfying the key filtering criterion for a HGT event.

- ### Installation

- git clone https://github.com/ChenSizhe13893461199/Fecal_Microbiota_Transplantation_HGT_Project.git
- cd Fecal_Microbiota_Transplantation_HGT_Project
- conda env create -f envs/environment.yml
- conda activate HGTector
