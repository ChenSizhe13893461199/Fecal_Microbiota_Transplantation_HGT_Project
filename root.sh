#!/bin/bash
# root.sh - Merge contig names with their assigned species names
# Output:
#   result/${post}_name.txt   -> post_contig_base <tab> species
#   result/${donor}_name.txt  -> donor_contig_base <tab> species

INPUT_XLSX="FMT_list.xlsx"
if [ ! -f "$INPUT_XLSX" ]; then
    echo "Error: Input file $INPUT_XLSX does not exist!"
    exit 1
fi

# Extract sample lists from Excel file
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

for i in "${!pre_samples[@]}"; do
    post="${post_samples[i]}"
    donor="${donor_samples[i]}"
    echo "Processing sample: $post (Donor: $donor)"

    # ============================================================
    # 1. Recipient: post_contig -> pre_contig -> species
    # ============================================================
    post_names="result/${post}_recipient_names.txt"
    post_species_map="result/${post}_recipient.species.map"
    post_contig1="HGT/${post}_contig1.txt"
    post_out="result/${post}_name.txt"

    if [ -f "$post_names" ] && [ -f "$post_species_map" ] && [ -f "$post_contig1" ]; then
        echo "  Generating recipient name->species mapping ..."

        awk -F'\t' -v OFS='\t' -v names_file="$post_names" '
        NR==FNR {
            pre2sp[$1] = $2
            next
        }
        {
            post_c = $1
            pre_c  = $2
            if (!(post_c in post2pre)) {
                post2pre[post_c] = pre_c
            }
        }
        END {
            while ((getline line < names_file) > 0) {
                split(line, a, "\t")
                pc = a[1]
                pr = (pc in post2pre) ? post2pre[pc] : ""
                sp = (pr != "" && pr in pre2sp) ? pre2sp[pr] : ""
                print pc, sp
            }
        }
        ' "$post_species_map" "$post_contig1" > "$post_out"

        echo "  Generated $post_out"
    else
        echo "  Warning: recipient files missing, creating empty $post_out"
        > "$post_out"
    fi

    # ============================================================
    # 2. Donor: donor_contig -> species
    # ============================================================
    donor_names="result/${donor}_donor_names.txt"
    donor_species_map="result/${donor}_donor.species.map"
    donor_out="result/${donor}_name.txt"

    if [ -f "$donor_names" ] && [ -f "$donor_species_map" ]; then
        echo "  Generating donor name->species mapping ..."

        awk -F'\t' -v OFS='\t' '
        NR==FNR {
            don2sp[$1] = $2
            next
        }
        {
            dc = $1
            sp = (dc in don2sp) ? don2sp[dc] : ""
            print dc, sp
        }
        ' "$donor_species_map" "$donor_names" > "$donor_out"

        echo "  Generated $donor_out"
    else
        echo "  Warning: donor files missing, creating empty $donor_out"
        > "$donor_out"
    fi
done

echo "All samples processed successfully!"
