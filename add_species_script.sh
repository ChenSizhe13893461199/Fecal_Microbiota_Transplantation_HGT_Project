#!/bin/bash
# add_species.sh - Add recipient and donor species information to HGT results

# Configuration
# Path to Kraken2 database (modify according to your setup)
KRAKEN2_DB="kracken2/"
# Path to Taxonkit data directory (contains taxonomy files) (modify according to your setup)
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

# Extract sample names using a Python one‑liner
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
    echo "========================================="
    echo "处理样本: $post (Donor: $donor)"

    hgt_file="result/${post}_HGT_full.txt"
    recipient_fasta="result/${post}_HGT_recipient_contig.fasta"
    donor_fasta="result/${donor}_HGT_donor_contig.fasta"

    if [ ! -f "$hgt_file" ]; then
        echo "警告：$hgt_file 不存在，跳过该样本"
        continue
    fi
    if [ ! -f "$recipient_fasta" ]; then
        echo "警告：recipient FASTA 文件 $recipient_fasta 不存在，跳过该样本"
        continue
    fi
    if [ ! -f "$donor_fasta" ]; then
        echo "警告：donor FASTA 文件 $donor_fasta 不存在，跳过该样本"
        continue
    fi

    # 提取唯一的recipient contig基本名（去掉坐标和基因编号）
    echo "提取recipient contig名称列表..."
    tail -n +2 "$hgt_file" | cut -f6 | sed 's/_[0-9]*-[0-9]*_[0-9]*$//' | sort -u > "result/${post}_recipient_names.txt"
    # 提取唯一的donor contig基本名（去掉坐标）
    echo "提取donor contig名称列表..."
    tail -n +2 "$hgt_file" | cut -f7 | sed 's/_[0-9]*-[0-9]*$//' | sort -u > "result/${donor}_donor_names.txt"

    # 检查列表是否为空
    if [ ! -s "result/${post}_recipient_names.txt" ]; then
        echo "警告：recipient名称列表为空，跳过该样本"
        continue
    fi
    if [ ! -s "result/${donor}_donor_names.txt" ]; then
        echo "警告：donor名称列表为空，跳过该样本"
        continue
    fi

    # 对recipient FASTA运行kraken2
    echo "对recipient contig运行kraken2..."
    kraken2 --db "$KRAKEN2_DB" --threads 4 \
        --output "result/${post}_recipient.kraken" \
        --report "result/${post}_recipient.report" \
        "$recipient_fasta" || {
        echo "kraken2运行失败，跳过该样本"
        continue
    }

    # 提取taxID
    cut -f3 "result/${post}_recipient.kraken" > "result/${post}_recipient.taxids"
    # 用taxonkit获取物种名（取最后一个分号后的部分）
    echo "对recipient taxID进行物种注释..."
    taxonkit lineage "result/${post}_recipient.taxids" --data-dir "$TAXKIT_DATA" 2>/dev/null | \
        awk -F'\t' '{split($2,a,";"); print $1"\t"a[length(a)]}' > "result/${post}_recipient.species.map"

    # 对donor FASTA运行kraken2
    echo "对donor contig运行kraken2..."
    kraken2 --db "$KRAKEN2_DB" --threads 4 \
        --output "result/${donor}_donor.kraken" \
        --report "result/${donor}_donor.report" \
        "$donor_fasta" || {
        echo "kraken2运行失败，跳过该样本"
        continue
    }

    cut -f3 "result/${donor}_donor.kraken" > "result/${donor}_donor.taxids"
    echo "对donor taxID进行物种注释..."
    taxonkit lineage "result/${donor}_donor.taxids" --data-dir "$TAXKIT_DATA" 2>/dev/null | \
        awk -F'\t' '{split($2,a,";"); print $1"\t"a[length(a)]}' > "result/${donor}_donor.species.map"

    # 构建关联数组
    declare -A recipient_sp
    declare -A donor_sp
    while IFS=$'\t' read -r seqid species; do
        recipient_sp["$seqid"]="$species"
    done < "result/${post}_recipient.species.map"
    while IFS=$'\t' read -r seqid species; do
        donor_sp["$seqid"]="$species"
    done < "result/${donor}_donor.species.map"

    # 处理HGT_full.txt，追加两列
    outfile="final/${post}_HGT_full.txt"
    # 写入表头
    echo -e "GC_HGT\tGC_origin\tGene_Description\tModule_Classification\tTaxonomic_Species\tRecipient_Contig\tDonor_Contig\trate\tlength\trecipient_species\tdonor_species" > "$outfile"

    # 逐行处理
    tail -n +2 "$hgt_file" | while IFS=$'\t' read -r gc_hgt gc_origin desc mod_class tax_spec rec_contig don_contig rate length _; do
        # 提取基本名
        rec_base=$(echo "$rec_contig" | sed 's/_[0-9]*-[0-9]*_[0-9]*$//')
        don_base=$(echo "$don_contig" | sed 's/_[0-9]*-[0-9]*$//')
        # 获取物种名（如果没有则留空）
        rec_species="${recipient_sp[$rec_base]:-}"
        don_species="${donor_sp[$don_base]:-}"
        # 输出
        echo -e "${gc_hgt}\t${gc_origin}\t${desc}\t${mod_class}\t${tax_spec}\t${rec_contig}\t${don_contig}\t${rate}\t${length}\t${rec_species}\t${don_species}"
    done >> "$outfile"

    echo "完成：$outfile"

    # 清理临时文件
    #rm -f "result/${post}_recipient.kraken" "result/${post}_recipient.taxids" "result/${post}_recipient.species.map" \
          #"result/${donor}_donor.kraken" "result/${donor}_donor.taxids" "result/${donor}_donor.species.map" \
          #"result/${post}_recipient_names.txt" "result/${donor}_donor_names.txt" \
          #"result/${post}_recipient.report" "result/${donor}_donor.report"
done

echo "所有样本处理完毕！"
