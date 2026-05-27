#!/bin/bash
# convert.sh
# Purpose: Extract recipient (pre-FMT) and donor contig sequences for species annotation.
# usage：./convert.sh

# Path to Kraken2 database (not used directly in this script, kept for consistency)
KRAKEN2_DB="kracken2/"
TAXKIT_DATA="tax/"
TRIMMED_DIR="FMT_HGT"           # Directory containing original filtered FASTAs (*_filter.fasta)

# Input Excel file listing sample triples (Pre-FMT, Donor, Post-FMT)
INPUT_XLSX="FMT_list.xlsx"

# Helper function: Check required tools
check_tool() {
    local tools=("python3" "kraken2" "taxonkit")
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            echo "Error: $tool not found. Please install it first."
            exit 1
        fi
    done
    # seqtk is optional; if missing, fallback to awk (slower)
    if ! command -v seqtk &> /dev/null; then
        echo "Warning: seqtk not installed. Will use awk for sequence extraction (slower)."
    fi
}

check_tool

# Extract sample lists from Excel file using a Python one-liner
python3 - <<END
import pandas as pd
import os
df = pd.read_excel("$INPUT_XLSX", engine='openpyxl')
pre_list = df['Pre-FMT'].dropna().astype(str).tolist()
donor_list = df['Donor'].dropna().astype(str).tolist()
post_list = df['Post-FMT'].dropna().astype(str).tolist()
with open('temp_pre.txt', 'w') as f: f.write('\n'.join(pre_list))
with open('temp_donor.txt', 'w') as f: f.write('\n'.join(donor_list))
with open('temp_post.txt', 'w') as f: f.write('\n'.join(post_list))
END

# Read the temporary files into bash arrays
pre_samples=($(cat temp_pre.txt))
donor_samples=($(cat temp_donor.txt))
post_samples=($(cat temp_post.txt))
rm -f temp_*.txt

# Verify that all three lists have the same number of samples
if [ ${#pre_samples[@]} -ne ${#donor_samples[@] } ] || [ ${#pre_samples[@]} -ne ${#post_samples[@]} ]; then
    echo "Error: Sample count mismatch in Excel file!"
    exit 1
fi

# Main loop: process each sample triple
for i in "${!pre_samples[@]}"; do
    pre="${pre_samples[i]}"
    donor="${donor_samples[i]}"
    post="${post_samples[i]}"

    echo "========================================="
    echo "Processing sample triple: Pre=$pre, Donor=$donor, Post=$post"

    hgt_file="result/${post}_HGT_full.txt"   # Table with HGT candidates
    post_contig1="HGT/${post}_contig1.txt"   # BLAST results: post vs pre (recipient)
    if [ ! -f "$hgt_file" ]; then
        echo "Warning: $hgt_file not found. Skipping this sample."
        continue
    fi
    if [ ! -f "$post_contig1" ]; then
        echo "Warning: $post_contig1 not found. Cannot map post contigs to pre. Skipping."
        continue
    fi

    # Step 1: Extract unique donor contig base names from HGT_full.txt
    # Remove the trailing "_start-end" to get the base contig name.
    echo "Extracting donor contig base names..."
    tail -n +2 "$hgt_file" | cut -f7 | sed 's/_[0-9]*-[0-9]*$//' | sort -u > "result/${donor}_donor_names.txt"

    # Step 2: Extract unique recipient (post-FMT) contig base names
    echo "Extracting recipient (post) contig base names..."
    tail -n +2 "$hgt_file" | cut -f6 | sed 's/_[0-9]*-[0-9]*_[0-9]*$//' | sort -u > "result/${post}_post_names.txt"

    # Step 3: Map each post contig to its corresponding pre-FMT contig
    echo "Mapping post contigs to pre contigs using ${post_contig1}..."
    # Build an associative array: post_contig -> pre_contig.
    declare -A post2pre
    while IFS=$'\t' read -r post_contig pre_contig rest; do
        # Only store the first mapping encountered for each post contig
        if [[ -z "${post2pre[$post_contig]}" ]]; then
            post2pre["$post_contig"]="$pre_contig"
        fi
    done < "$post_contig1"

    # For each unique post contig name, fetch its pre contig name, then deduplicate
    > "result/${post}_pre_names.txt"
    while read post_name; do
        pre_name="${post2pre[$post_name]}"
        if [ -n "$pre_name" ]; then
            echo "$pre_name" >> "result/${post}_pre_names.txt"
        else
            echo "Warning: post contig $post_name has no mapping in $post_contig1. Skipping."
        fi
    done < "result/${post}_post_names.txt"
    sort -u "result/${post}_pre_names.txt" -o "result/${post}_pre_names.txt"

    # Step 4: Extract actual nucleotide sequences for pre-FMT and donor contigs
    pre_fasta="${TRIMMED_DIR}/${pre}_filter.fasta"
    donor_fasta="${TRIMMED_DIR}/${donor}_filter.fasta"
    pre_out="result/${post}_HGT_recipient_contig.fasta"
    donor_out="result/${donor}_HGT_donor_contig.fasta"

    # Verify that original FASTA files exist
    if [ ! -f "$pre_fasta" ]; then
        echo "Error: donor FASTA $donor_fasta not found. Skipping sample."
        continue
    fi
    if [ ! -f "$donor_fasta" ]; then
        echo "Error: pre-FMT FASTA $pre_fasta not found. Skipping sample."
        continue
    fi

    # Extract pre-FMT sequences using seqtk (preferred) or awk/sed fallback
    echo "Extracting pre-FMT sequences to $pre_out"
    if [ -s "result/${post}_pre_names.txt" ]; then
        if command -v seqtk &> /dev/null; then
            # seqtk subseq: fast extraction by exact header match
            seqtk subseq "$pre_fasta" "result/${post}_pre_names.txt" > "$pre_out"
        else
            awk 'NR==FNR{a[$1];next} /^>/ {h=substr($1,2); if(h in a) {print; flag=1; next}} flag' \
                "result/${post}_pre_names.txt" "$pre_fasta" > "$pre_out"
        fi
        # If the output is empty despite having names, try a slower sed-based method
        if [ ! -s "$pre_out" ]; then
            echo "Fallback method for pre-FMT extraction..."
            > "$pre_out"
            while read name; do
                sed -n "/^>$name$/,/^>/p" "$pre_fasta" | sed '$d' >> "$pre_out"
            done < "result/${post}_pre_names.txt"
        fi
    else
        echo "Warning: pre name list is empty. Skipping pre-FMT extraction."
        > "$pre_out"
    fi

    # Extract donor sequences using seqtk (preferred) or awk/sed fallback
    echo "Extracting donor sequences to $donor_out"
    if [ -s "result/${donor}_donor_names.txt" ]; then
        if command -v seqtk &> /dev/null; then
            seqtk subseq "$donor_fasta" "result/${donor}_donor_names.txt" > "$donor_out"
        else
            awk 'NR==FNR{a[$1];next} /^>/ {h=substr($1,2); if(h in a) {print; flag=1; next}} flag' \
                "result/${donor}_donor_names.txt" "$donor_fasta" > "$donor_out"
        fi
        if [ ! -s "$donor_out" ]; then
            echo "Fallback method for donor extraction..."
            > "$donor_out"
            while read name; do
                sed -n "/^>$name$/,/^>/p" "$donor_fasta" | sed '$d' >> "$donor_out"
            done < "result/${donor}_donor_names.txt"
        fi
    else
        echo "Warning: donor name list is empty. Skipping donor extraction."
        > "$donor_out"
    fi
done

echo "All samples processed successfully."
