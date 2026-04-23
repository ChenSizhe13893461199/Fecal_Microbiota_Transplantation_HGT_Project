#!/bin/bash
INPUT_XLSX="FMT_list.xlsx" #It contains information of pre-FMT recipient, post-FMT recipient, and Donor Name and Orders
TRIMMED_DIR="FMT_HGT" #The Directory where all metagenomic assembled files
BLAST_DB_DIR="blast_dbs" # The Directory where temporary files of blasting
RESULTS_DIR="blast_results" # The Directory where temporary files of analysis

# -------------------------- initial checking for necessary tools --------------------------
check_tool() {
    local tools=("python3" "makeblastdb" "blastn")
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            echo "error：necessary tools not found $tool，please install before running！"
            exit 1
        fi
    done
}

#check input
if [ ! -f "$INPUT_XLSX" ]; then  # repair
    echo "error：input file $INPUT_XLSX not exist！"
    exit 1
fi

# -------------------------- initialize environment of tasks for analysis--------------------------
check_tool
mkdir -p "$TRIMMED_DIR" "$BLAST_DB_DIR" "$RESULTS_DIR"

# -------------------------- extract sample list --------------------------
python3 - <<END  # 
import pandas as pd
import os

df = pd.read_excel(os.path.abspath("$INPUT_XLSX"), engine='openpyxl')  #reading the input FMT file 
pre_list = df['Pre-FMT'].dropna().astype(str).tolist() #pre-FMT recipient
donor_list = df['Donor'].dropna().astype(str).tolist() #donor
post_list = df['Post-FMT'].dropna().astype(str).tolist() # post-FMT recipient

with open('temp_pre.txt', 'w') as f: f.write('\n'.join(pre_list))
with open('temp_donor.txt', 'w') as f: f.write('\n'.join(donor_list))
with open('temp_post.txt', 'w') as f: f.write('\n'.join(post_list))
END

# 
pre_samples=($(cat temp_pre.txt))
donor_samples=($(cat temp_donor.txt))
post_samples=($(cat temp_post.txt))

rm -f temp_*.txt

# check sample consistency
if [ ${#pre_samples[@]} -ne ${#donor_samples[@]} ] || [ ${#pre_samples[@]} -ne ${#post_samples[@]} ]; then
    echo "Error：Sample Numbers are Inconsistent！"
    exit 1
fi

# -------------------------- main program workflow --------------------------
for i in "${!pre_samples[@]}"; do
    pre="${pre_samples[i]}"
    donor="${donor_samples[i]}"
    post="${post_samples[i]}"
    
    pre_fasta="${TRIMMED_DIR}/${pre}_filter.fasta"
    donor_fasta="${TRIMMED_DIR}/${donor}_filter.fasta"
    post_fasta="${TRIMMED_DIR}/${post}_filter.fasta"

    # 
    missing_files=()
    [ ! -f "$pre_fasta" ] && missing_files+=("Pre-FMT: $pre_fasta")
    [ ! -f "$donor_fasta" ] && missing_files+=("Donor: $donor_fasta")
    [ ! -f "$post_fasta" ] && missing_files+=("Post-FMT: $post_fasta")
    
    if [ ${#missing_files[@]} -gt 0 ]; then
        echo "warning：skip $((i+1)) th sample，file is not found："
        printf '  - %s\n' "${missing_files[@]}"
        continue
    fi

    echo "======== processing the $((i+1)) th sample ========"""
    echo "Pre-FMT: $pre | Donor: $donor | Post-FMT: $post"

    # The following codes start to do the 1st round of sequence alignment between FMT donor and post-FMT recipient (database construction)
    donor_db="${BLAST_DB_DIR}/${donor}_donor_db"
    if [ ! -d "$donor_db" ]; then
        echo "create donor database：$donor_db"
        makeblastdb -in "$donor_fasta" -dbtype nucl -out "$donor_db" || {
            echo "error：donor library construction failed！"
            continue
        }
    fi
    # The following codes start to do the 1st round of sequence alignment between pre-FMT recipient and post-FMT recipient (database construction)
    recipient_db="${BLAST_DB_DIR}/${pre}_recipient_db"
    if [ ! -d "$recipient_db" ]; then
        echo "create recipient database：$recipient_db"
        makeblastdb -in "$pre_fasta" -dbtype nucl -out "$recipient_db" || {
            echo "error：recipient library construction failed！"
            continue
        }
    fi

    # The following codes start to do the 1st round of sequence alignment between FMT donor and post-FMT recipient
    #export BLASTDB="${BLAST_DB_DIR}:${BLASTDB}"
    blastn -query "$post_fasta" \
           -db "$donor_db" \
           -outfmt 6 \
           -max_target_seqs 1 \
           -out "${RESULTS_DIR}/${post}_blast_donor.txt" || {
        echo "warning：Post-FMT1 BLAST failed！"
    }
    #The following codes start to do the 1st round of sequence alignment between pre-FMT recipient and post-FMT recipient
    #export BLASTDB="${BLAST_DB_DIR}:${BLASTDB}"
    blastn -query "$post_fasta" \
           -db "$recipient_db" \
           -outfmt 6 \
           -max_target_seqs 1 \
           -out "${RESULTS_DIR}/${post}_blast_recipient.txt" || {
        echo "warning：Post-FMT2 BLAST failed！"
    }

    # The following script "blastn_process.py" was used to generate summary of alignment, 
    # in which -t (e.g. 0.9 represents at least 90.0% of full length) controls the filtering threshold for the ratio of alignment length
    
    python blastn_process.py -i "${RESULTS_DIR}/${post}_blast_donor.txt" -o "${RESULTS_DIR}/${post}_gt_donor_results.txt" -t 0.9 
    #find those post-FMT contigs with at least 90.0% alignment coverage (identity ≥ 99.0%， e-value ≤ 10^(-10)) to donor contigs
    
    python blastn_process.py -i "${RESULTS_DIR}/${post}_blast_recipient.txt" -o "${RESULTS_DIR}/${post}_gt_recipient_results.txt" -t 0.9
    #find those post-FMT contigs with at least 90.0% alignment coverage (identity ≥ 99.0%， e-value ≤ 10^(-10)) to recipient contigs

    python classifer.py -i "$post_fasta" "${RESULTS_DIR}/${post}_gt_donor_results.txt" "${RESULTS_DIR}/${post}_gt_recipient_results.txt" -o "${RESULTS_DIR}/${donor}_contigs.fasta" "${RESULTS_DIR}/${post}_contigs.fasta" "${RESULTS_DIR}/${donor}_${post}_contigs.fasta" -s 100
    #using the two .txt file generated above to classify contigs (among contigs in post-FMT recipients) from donor, from pre-FMT recipient, both (from donor and pre-FMT recipient), 
    #and others (not specified in output code above, classifer.py will generate 4 files including the others (suspected contigs may contain HGT region)).

    #removing some redundnat files
    rm "${RESULTS_DIR}/${post}_blast_donor.txt"
    rm "${RESULTS_DIR}/${post}_blast_recipient.txt"
    rm "${RESULTS_DIR}/${post}_gt_donor_results.txt"
    rm "${RESULTS_DIR}/${post}_gt_recipient_results.txt"
    rm "${RESULTS_DIR}/${donor}_contigs.fasta"
    rm "${RESULTS_DIR}/${post}_contigs.fasta"
    rm "${RESULTS_DIR}/${donor}_${post}_contigs.fasta"
    rm "${RESULTS_DIR}/${donor}_${post}_contigs_info.txt"
    rm "${RESULTS_DIR}/${post}_contigs_info.txt"
    rm "${RESULTS_DIR}/${donor}_contigs_info.txt"
    rm "${RESULTS_DIR}/${donor}_${post}_contigs_other_donor_info.txt"
    rm "${RESULTS_DIR}/${donor}_${post}_contigs_other_recipient_info.txt"

    # The following codes start to do the 2nd round of sequence alignment between others (suspected contigs may contain HGT region) and pre-FMT recipient (database construction)
    #export BLASTDB="${BLAST_DB_DIR}:${BLASTDB}"
    blastn -query "${RESULTS_DIR}/${donor}_${post}_contigs_other.fasta" \
           -db "$recipient_db" \
           -outfmt 6 \
           -max_target_seqs 1 \
           -out "${RESULTS_DIR}/${post}_blast_recipient1.txt" || {
        echo "warning：Post-FMT2 BLAST failed！"
    }

    # The following codes start to do the 2nd round of sequence alignment between others (suspected contigs may contain HGT region) and donor (database construction)
    #export BLASTDB="${BLAST_DB_DIR}:${BLASTDB}"
    blastn -query "${RESULTS_DIR}/${donor}_${post}_contigs_other.fasta" \
           -db "$donor_db" \
           -outfmt 6 \
           -max_target_seqs 1 \
           -out "${RESULTS_DIR}/${post}_blast_donor1.txt" || {
        echo "warning：Post-FMT2 BLAST failed！"
    }

    # The following codes implement the 1st round detection of potential HGT regions on others (suspected contigs may contain HGT region)
    # using the other fasta file as root, the HGT.py search for any potential contigs shared high homologous rate (identity ≥ 99.0%， e-value ≤ 10^(-10)) to donbor contig
    # within certain length range ([500bp, 50% of full contig length]). Additionally, all suspected homologous regions must not overlap with any region within
    # the contig from pre‑FMT recipient, with flanking and non-HGT regions sufficiently aligned to contig in pre-FMT samples (≥ 99.0%, e-value ≤ 10^(-10))
    # Only contigs satisfying all criteria were retained as candidate harboring HGT events
    
    python HGT.py \
        --other_fasta "${RESULTS_DIR}/${donor}_${post}_contigs_other.fasta" \
        --recipient_blast "${RESULTS_DIR}/${post}_blast_recipient1.txt" \
        --donor_blast "${RESULTS_DIR}/${post}_blast_donor1.txt" \
        --donor_fasta "${TRIMMED_DIR}/${donor}_filter.fasta" \
        --output_prefix "$post" "$donor"
    mkdir -p HGT

    #removing some redundnat files
    mv "${post}_aligned.fasta" HGT/
    mv "${donor}_aligned.fasta" HGT/
    mv "${post}_contig.fasta" HGT/
    mv "${donor}_contig.fasta" HGT/
    mv "${post}_contig1.txt" HGT/
    mv "${donor}_contig1.txt" HGT/

    mkdir output
    mkdir HGT

    #encoding gene prediction and processing
    prodigal -i "HGT/${post}_aligned.fasta" -o "output/${post}_aligned.fasta.gbk" -a "output/${post}_aligned_proteins.faa" -d "output/${post}_aligned_nucleotides.faa" -p meta
    source activate eggnog

    #functional annotation of encoding genes predicted above
    emapper.py -i "output/${post}_aligned_proteins.faa" -o "${post}_donor_HGT" --output_dir HGT  -m diamond --data_dir data/ --cpu 8
    source activate base 
done

# quality control steps
mkdir result

# 1st round of quality control
# The qc_annotation.py script is a post‑processing quality control step in the HGT (Horizontal Gene Transfer) detection pipeline. 
# Its main function is to integrate information from multiple intermediate files and generate a structured summary
# table for each sample pair (Donor/Post‑FMT). Specifically, it has the following functions:
# 1) Parses eggnog‑mapper annotations (.emapper.annotations) to retrieve gene descriptions, COG categories, taxonomic assignments, and genomic coordinates
# 2) for each predicted gene in the recipient’s HGT‑candidate region. Reads donor BLAST results (donor_contig1.txt) to extract high‑scoring segment pairs (HSPs)
# between donor and recipient contigs.
# 3) Calculates GC content for both the HGT region (from *_aligned.fasta) and the original donor contig (from donor_contig fasta file).
# 4) Matches each gene to the overlapping HSP on the donor side, extracting alignment identity (rate) and alignment length
python qc_annotation.py 
mkdir -p final
source activate kraken2
./convert.sh
./add_species_script.sh
./root.sh
./fill_species.sh
rm final/*_HGT_full.txt
#rm -rf result/*
rm -rf output/*_aligned.fasta.gbk
source activate base
mkdir filter
python filterchecking.py
