#!/bin/bash
# root.sh - Fix species mapping (using awk to avoid bash pitfalls)

INPUT_XLSX="FMT_list.xlsx"
if [ ! -f "$INPUT_XLSX" ]; then
    echo "Error: Input file $INPUT_XLSX does not exist!"
    exit 1
fi

# Extract sample lists
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

    # ========== 1. Recipient part ==========
    post_names="result/${post}_recipient_names.txt"
    post_contig1="HGT/${post}_contig1.txt"
    pre_names="result/${post}_pre_names.txt"
    pre_species_map="result/${post}_recipient.species.map"
    post_out="result/${post}_name.txt"

    if [ -f "$post_names" ] && [ -f "$post_contig1" ] && [ -f "$pre_names" ] && [ -f "$pre_species_map" ]; then
        echo "  Generating recipient species mapping ..."
        
        # 1) Build pre_contig to species mapping file (pre2species.txt)
        #    pre_names and pre_species_map have same number of lines, in same order. Merge and take column 1 and 3 (species name)
        paste "$pre_names" "$pre_species_map" | awk -F'\t' '{print $1"\t"$3}' > "tmp_${post}_pre2species.txt"
        
        # 2) For each post contig in post_contig1.txt, take the first matching pre contig, then look up species
        #    Output post_contig -> species mapping (keep only first occurrence)
        awk -F'\t' '
        BEGIN {
            # Read pre2species mapping
            while ((getline < "tmp_'${post}'_pre2species.txt") > 0) {
                pre2sp[$1] = $2;
            }
            close("tmp_'${post}'_pre2species.txt");
        }
        {
            post = $1;
            pre = $2;
            if (!(post in seen)) {
                seen[post] = 1;
                sp = (pre in pre2sp) ? pre2sp[pre] : "";
                print post "\t" sp;
            }
        }
        ' "$post_contig1" > "tmp_${post}_post2species.txt"
        
        # 3) Left join species information with post_names.txt (empty if no match)
        awk -F'\t' '
        BEGIN {
            while ((getline < "tmp_'${post}'_post2species.txt") > 0) {
                post2sp[$1] = $2;
            }
            close("tmp_'${post}'_post2species.txt");
        }
        {
            sp = ($1 in post2sp) ? post2sp[$1] : "";
            print $1 "\t" sp;
        }
        ' "$post_names" > "$post_out"
        
        rm -f "tmp_${post}_pre2species.txt" "tmp_${post}_post2species.txt"
        echo "  Generated $post_out"
    else
        echo "  Warning: recipient files missing, skipping"
    fi

    # ========== 2. Donor part ==========
    donor_names="result/${donor}_donor_names.txt"
    donor_species_map="result/${donor}_donor.species.map"
    donor_out="result/${donor}_name.txt"

    if [ -f "$donor_names" ] && [ -f "$donor_species_map" ]; then
        # donor_species.map may contain taxID and species; need to take species corresponding to donor_names line by line
        # Assume both files have same line count; simply paste and take column 1 and 3 (species)
        paste "$donor_names" "$donor_species_map" | awk -F'\t' '{print $1"\t"$3}' > "$donor_out"
        echo "  Generated $donor_out"
    else
        echo "  Warning: donor files missing, skipping"
    fi
done

echo "All samples processed!"
