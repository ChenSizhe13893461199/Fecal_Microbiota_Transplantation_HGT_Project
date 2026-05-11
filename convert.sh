#!/bin/bash
# add_species.sh
# Purpose: Extract recipient (pre-FMT) and donor contig sequences for species annotation.
# usage：./add_species.sh

# Path to Kraken2 database (not used directly in this script, kept for consistency)
KRAKEN2_DB="kracken2/"
TAXKIT_DATA="tax/"
TRIMMED_DIR="FMT_HGT"          # 原始过滤 fasta 所在目录

# Input Excel file listing sample triples (Pre-FMT, Donor, Post-FMT)
INPUT_XLSX="FMT_list.xlsx"

# Helper function: Check required tools
check_tool() {
    local tools=("python3" "kraken2" "taxonkit")
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            echo "错误：未找到 $tool，请先安装。"
            exit 1
        fi
    done
    # seqtk 可选
    if ! command -v seqtk &> /dev/null; then
        echo "警告：seqtk 未安装，将使用 awk 提取序列（速度较慢）。"
    fi
}

check_tool

# 从 Excel 提取样本列表
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

pre_samples=($(cat temp_pre.txt))
donor_samples=($(cat temp_donor.txt))
post_samples=($(cat temp_post.txt))
rm -f temp_*.txt

# 检查数量一致
if [ ${#pre_samples[@]} -ne ${#donor_samples[@] } ] || [ ${#pre_samples[@]} -ne ${#post_samples[@]} ]; then
    echo "错误：样本数量不一致！"
    exit 1
fi

# 主循环
for i in "${!pre_samples[@]}"; do
    pre="${pre_samples[i]}"
    donor="${donor_samples[i]}"
    post="${post_samples[i]}"

    echo "========================================="
    echo "处理样本对：Pre=$pre, Donor=$donor, Post=$post"

    hgt_file="result/${post}_HGT_full.txt"
    post_contig1="HGT/${post}_contig1.txt"   # post vs pre 的 BLAST 结果

    if [ ! -f "$hgt_file" ]; then
        echo "警告：$hgt_file 不存在，跳过"
        continue
    fi
    if [ ! -f "$post_contig1" ]; then
        echo "警告：$post_contig1 不存在，无法映射 pre contig，跳过"
        continue
    fi

    # ---------- 提取 donor 基本名称 ----------
    echo "提取 donor contig 基本名称..."
    tail -n +2 "$hgt_file" | cut -f7 | sed 's/_[0-9]*-[0-9]*$//' | sort -u > "result/${donor}_donor_names.txt"

    # ---------- 提取 recipient 基本名称 ----------
    echo "提取 recipient (post) contig 基本名称..."
    tail -n +2 "$hgt_file" | cut -f6 | sed 's/_[0-9]*-[0-9]*_[0-9]*$//' | sort -u > "result/${post}_post_names.txt"

    # ---------- 映射 post contig 到 pre contig ----------
    echo "从 ${post_contig1} 映射 post contig 到 pre contig..."
    # 创建一个关联数组映射 post -> pre（取第一个匹配的行）
    declare -A post2pre
    while IFS=$'\t' read -r post_contig pre_contig rest; do
        # 只取第一个匹配，如果已存在则跳过
        if [[ -z "${post2pre[$post_contig]}" ]]; then
            post2pre["$post_contig"]="$pre_contig"
        fi
    done < "$post_contig1"

    # 对于每个 post contig 名称，获取对应的 pre contig 名称，并去重
    > "result/${post}_pre_names.txt"
    while read post_name; do
        pre_name="${post2pre[$post_name]}"
        if [ -n "$pre_name" ]; then
            echo "$pre_name" >> "result/${post}_pre_names.txt"
        else
            echo "警告：post contig $post_name 在 $post_contig1 中未找到映射，跳过"
        fi
    done < "result/${post}_post_names.txt"
    sort -u "result/${post}_pre_names.txt" -o "result/${post}_pre_names.txt"

    # ---------- 提取 pre 和 donor 的序列 ----------
    pre_fasta="${TRIMMED_DIR}/${pre}_filter.fasta"
    donor_fasta="${TRIMMED_DIR}/${donor}_filter.fasta"
    pre_out="result/${post}_HGT_recipient_contig.fasta"
    donor_out="result/${donor}_HGT_donor_contig.fasta"

    # 检查原始 fasta 是否存在
    if [ ! -f "$pre_fasta" ]; then
        echo "错误：pre 原始 fasta $pre_fasta 不存在，跳过"
        continue
    fi
    if [ ! -f "$donor_fasta" ]; then
        echo "错误：donor 原始 fasta $donor_fasta 不存在，跳过"
        continue
    fi

    # 提取 pre 序列
    echo "提取 pre 序列到 $pre_out"
    if [ -s "result/${post}_pre_names.txt" ]; then
        if command -v seqtk &> /dev/null; then
            seqtk subseq "$pre_fasta" "result/${post}_pre_names.txt" > "$pre_out"
        else
            awk 'NR==FNR{a[$1];next} /^>/ {h=substr($1,2); if(h in a) {print; flag=1; next}} flag' \
                "result/${post}_pre_names.txt" "$pre_fasta" > "$pre_out"
        fi
        # 如果为空，尝试备用
        if [ ! -s "$pre_out" ]; then
            echo "备用方法提取 pre 序列..."
            > "$pre_out"
            while read name; do
                sed -n "/^>$name$/,/^>/p" "$pre_fasta" | sed '$d' >> "$pre_out"
            done < "result/${post}_pre_names.txt"
        fi
    else
        echo "警告：pre 名称列表为空，跳过 pre 序列提取"
        > "$pre_out"
    fi

    # 提取 donor 序列
    echo "提取 donor 序列到 $donor_out"
    if [ -s "result/${donor}_donor_names.txt" ]; then
        if command -v seqtk &> /dev/null; then
            seqtk subseq "$donor_fasta" "result/${donor}_donor_names.txt" > "$donor_out"
        else
            awk 'NR==FNR{a[$1];next} /^>/ {h=substr($1,2); if(h in a) {print; flag=1; next}} flag' \
                "result/${donor}_donor_names.txt" "$donor_fasta" > "$donor_out"
        fi
        if [ ! -s "$donor_out" ]; then
            echo "备用方法提取 donor 序列..."
            > "$donor_out"
            while read name; do
                sed -n "/^>$name$/,/^>/p" "$donor_fasta" | sed '$d' >> "$donor_out"
            done < "result/${donor}_donor_names.txt"
        fi
    else
        echo "警告：donor 名称列表为空，跳过 donor 序列提取"
        > "$donor_out"
    fi
done

echo "所有样本处理完成。"
