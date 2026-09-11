#!/bin/bash
# add_species.sh - Add recipient and donor species information to HGT results

# Configuration
KRAKEN2_DB="kracken2/"
TAXKIT_DATA="tax/"

# Check required tools
command -v kraken2 >/dev/null 2>&1 || { echo "Error: kraken2 not found"; exit 1; }
command -v taxonkit >/dev/null 2>&1 || { echo "Error: taxonkit not found"; exit 1; }

# Read sample list from Excel file
INPUT_XLSX="FMT_list.xlsx"
if [ ! -f "$INPUT_XLSX" ]; then
    echo "Error: Input file $INPUT_XLSX does not exist!"
    exit 1
fi

# Extract sample names using a Python one-liner
python3 - <<END
import pandas as pd
df = pd.read_excel("$INPUT_XLSX", engine='openpyxl')
pre_list = df['Pre-FMT'].dropna().astype(str).tolist()
donor_list = df['Donor'].dropna().astype(str).tolist()
post_list = df['Post-FMT'].dropna().astype(str).tolist()
with open('temp_pre.txt', 'w') as f: f.write('\n'.join(pre_list))
with open('temp_donor.txt', 'w') as f: f.write('\n'.join(donor_list))
with open('temp_post.txt', 'w') as f: f.write('\n'.join(post_list))
END

pre_samples=($(cat temp_pre.txt))
donor_samples=($(cat temp_donor.txt))
post_samples=($(cat temp_post.txt))
rm -f temp_*.txt

if [ ${#pre_samples[@]} -ne ${#donor_samples[@]} ] || [ ${#pre_samples[@]} -ne ${#post_samples[@]} ]; then
    echo "Error: Sample count mismatch!"
    exit 1
fi

mkdir -p final

for i in "${!pre_samples[@]}"; do
    post="${post_samples[i]}"
    donor="${donor_samples[i]}"
    pre="${pre_samples[i]}"          # MOD 0: 需要 pre 变量
    echo "========================================="
    echo "Processing sample: $post (Donor: $donor, Pre: $pre)"

    hgt_file="result/${post}_HGT_full.txt"
    recipient_fasta="result/${post}_HGT_recipient_contig.fasta"
    donor_fasta="result/${donor}_HGT_donor_contig.fasta"

    if [ ! -f "$hgt_file" ]; then
        echo "Warning: $hgt_file not found, skipping this sample"
        continue
    fi
    if [ ! -f "$recipient_fasta" ]; then
        echo "Warning: recipient FASTA $recipient_fasta not found, skipping"
        continue
    fi
    if [ ! -f "$donor_fasta" ]; then
        echo "Warning: donor FASTA $donor_fasta not found, skipping"
        continue
    fi

    # Extract unique recipient contig base names (post-FMT)
    echo "Extracting recipient contig names..."
    tail -n +2 "$hgt_file" | cut -f6 | sed 's/_[0-9]*-[0-9]*_[0-9]*$//' | sort -u > "result/${post}_recipient_names.txt"
    # Extract unique donor contig base names
    echo "Extracting donor contig names..."
    tail -n +2 "$hgt_file" | cut -f7 | sed 's/_[0-9]*-[0-9]*$//' | sort -u > "result/${donor}_donor_names.txt"

    if [ ! -s "result/${post}_recipient_names.txt" ]; then
        echo "Warning: recipient name list empty, skipping sample"
        continue
    fi
    if [ ! -s "result/${donor}_donor_names.txt" ]; then
        echo "Warning: donor name list empty, skipping sample"
        continue
    fi

    # ============================================================
    # Recipient side: Kraken2 on pre-FMT contigs (provided by convert.sh)
    # ============================================================
    kraken2 --db "$KRAKEN2_DB" --threads 4 \
        --output "result/${post}_recipient.kraken" \
        --report "result/${post}_recipient.report" \
        "$recipient_fasta" || {
        echo "recipient kraken2 failed, skipping sample"
        continue
    }

    # MOD 1: keep contig_id (col2) + taxid (col3), filter taxid 0
    awk -F'\t' '$3 != 0 {print $2"\t"$3}' "result/${post}_recipient.kraken" \
        > "result/${post}_recipient.contig_taxid"

    # MOD 2: unique taxids -> taxonkit -> taxid -> species
    cut -f2 "result/${post}_recipient.contig_taxid" | sort -u > "result/${post}_recipient.taxids"
    echo "Annotating recipient taxIDs with species names..."
    taxonkit lineage "result/${post}_recipient.taxids" --data-dir "$TAXKIT_DATA" 2>/dev/null | \
        awk -F'\t' '{split($2,a,";"); print $1"\t"a[length(a)]}' \
        > "result/${post}_recipient.taxid2species"

    # MOD 3: join contig_taxid + taxid2species => contig -> species (pre-FMT contig)
    join -1 2 -2 1 -t $'\t' \
        <(sort -k2,2 "result/${post}_recipient.contig_taxid") \
        <(sort -k1,1 "result/${post}_recipient.taxid2species") \
        | awk -F'\t' '{print $2"\t"$3}' > "result/${post}_recipient.species.map"

    # ============================================================
    # Donor side: Kraken2 on donor contigs (already in donor FASTA)
    # ============================================================
    echo "Annotating donor taxIDs with species names..."
    kraken2 --db "$KRAKEN2_DB" --threads 4 \
        --output "result/${donor}_donor.kraken" \
        --report "result/${donor}_donor.report" \
        "$donor_fasta" || {
        echo "donor kraken2 failed, skipping sample"
        continue
    }

    awk -F'\t' '$3 != 0 {print $2"\t"$3}' "result/${donor}_donor.kraken" \
        > "result/${donor}_donor.contig_taxid"
    cut -f2 "result/${donor}_donor.contig_taxid" | sort -u > "result/${donor}_donor.taxids"
    taxonkit lineage "result/${donor}_donor.taxids" --data-dir "$TAXKIT_DATA" 2>/dev/null | \
        awk -F'\t' '{split($2,a,";"); print $1"\t"a[length(a)]}' \
        > "result/${donor}_donor.taxid2species"

    join -1 2 -2 1 -t $'\t' \
        <(sort -k2,2 "result/${donor}_donor.contig_taxid") \
        <(sort -k1,1 "result/${donor}_donor.taxid2species") \
        | awk -F'\t' '{print $2"\t"$3}' > "result/${donor}_donor.species.map"

    # ============================================================
    # MOD 4: Build species lookup arrays
    # recipient: species.map 的键是 pre-FMT contig，
    #             post -> pre 映射HGT/${post}_contig1.txt，
    #             pre -> species， post -> species
    # ============================================================
    declare -A pre2sp
    while IFS=$'\t' read -r seqid species; do
        pre2sp["$seqid"]="$species"
    done < "result/${post}_recipient.species.map"

    declare -A post2pre
    if [ -f "HGT/${post}_contig1.txt" ]; then
        while IFS=$'\t' read -r post_c pre_c rest; do
            if [ -z "${post2pre[$post_c]}" ]; then
                post2pre["$post_c"]="$pre_c"
            fi
        done < "HGT/${post}_contig1.txt"
    fi

    declare -A recipient_sp
    for post_c in "${!post2pre[@]}"; do
        pre_c="${post2pre[$post_c]}"
        recipient_sp["$post_c"]="${pre2sp[$pre_c]:-}"
    done

    declare -A donor_sp
    while IFS=$'\t' read -r seqid species; do
        donor_sp["$seqid"]="$species"
    done < "result/${donor}_donor.species.map"

    # ============================================================
    # Merge species info into HGT table
    # ============================================================
    outfile="final/${post}_HGT_full.txt"
    echo -e "GC_HGT\tGC_origin\tGene_Description\tModule_Classification\tTaxonomic_Species\tRecipient_Contig\tDonor_Contig\trate\tlength\trecipient_species\tdonor_species" > "$outfile"

    tail -n +2 "$hgt_file" | while IFS=$'\t' read -r gc_hgt gc_origin desc mod_class tax_spec rec_contig don_contig rate length _; do
        rec_base=$(echo "$rec_contig" | sed 's/_[0-9]*-[0-9]*_[0-9]*$//')
        don_base=$(echo "$don_contig" | sed 's/_[0-9]*-[0-9]*$//')
        rec_species="${recipient_sp[$rec_base]:-}"
        don_species="${donor_sp[$don_base]:-}"
        echo -e "${gc_hgt}\t${gc_origin}\t${desc}\t${mod_class}\t${tax_spec}\t${rec_contig}\t${don_contig}\t${rate}\t${length}\t${rec_species}\t${don_species}"
    done >> "$outfile"

    echo "Completed: $outfile"

    # Optional cleanup
    #rm -f "result/${post}_recipient.kraken" "result/${post}_recipient.taxids" \
    #      "result/${post}_recipient.contig_taxid" "result/${post}_recipient.taxid2species" \
    #      "result/${donor}_donor.kraken" "result/${donor}_donor.taxids" \
    #      "result/${donor}_donor.contig_taxid" "result/${donor}_donor.taxid2species"
done

echo "All samples processed successfully!"
