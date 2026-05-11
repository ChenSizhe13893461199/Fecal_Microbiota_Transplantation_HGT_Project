#!/bin/bash
# fill_species.sh - 填充 HGT 结果中的物种信息

INPUT_XLSX="FMT_list.xlsx"
if [ ! -f "$INPUT_XLSX" ]; then
    echo "错误：输入文件 $INPUT_XLSX 不存在！"
    exit 1
fi

# 提取样本列表
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
    echo "错误：样本数量不一致！"
    exit 1
fi

mkdir -p final

for i in "${!pre_samples[@]}"; do
    post="${post_samples[i]}"
    donor="${donor_samples[i]}"

    echo "处理样本: $post (Donor: $donor)"

    input_file="final/${post}_HGT_full.txt"
    output_file="final/${post}_HGT_statistics.txt"
    recipient_map="result/${post}_name.txt"
    donor_map="result/${donor}_name.txt"

    if [ ! -f "$input_file" ]; then
        echo "警告：输入文件 $input_file 不存在，跳过"
        continue
    fi
    if [ ! -f "$recipient_map" ]; then
        echo "警告：recipient 映射文件 $recipient_map 不存在，跳过"
        continue
    fi
    if [ ! -f "$donor_map" ]; then
        echo "警告：donor 映射文件 $donor_map 不存在，跳过"
        continue
    fi

    # 使用 awk 进行填充
    awk -F'\t' -v OFS='\t' '
    BEGIN {
        # 读取 recipient 映射
        while ((getline < "'$recipient_map'") > 0) {
            split($0, a, "\t");
            rec[a[1]] = a[2];
        }
        close("'$recipient_map'");
        # 读取 donor 映射
        while ((getline < "'$donor_map'") > 0) {
            split($0, a, "\t");
            don[a[1]] = a[2];
        }
        close("'$donor_map'");
    }
    {
        if (NR == 1) {
            # 表头，直接输出
            print;
        } else {
            # 提取 Recipient_Contig 基本名（第6列）和 Donor_Contig 基本名（第7列）
            rec_full = $6;
            don_full = $7;
            # 去掉坐标和基因编号部分
            gsub(/_[0-9]*-[0-9]*_[0-9]*$/, "", rec_full);
            gsub(/_[0-9]*-[0-9]*$/, "", don_full);
            # 获取物种名（如果不存在则留空）
            rec_species = (rec_full in rec) ? rec[rec_full] : "";
            don_species = (don_full in don) ? don[don_full] : "";
            # 更新第10列和第11列
            $10 = rec_species;
            $11 = don_species;
            print;
        }
    }' "$input_file" > "$output_file"

    echo "已生成 $output_file"
done

echo "所有样本处理完毕！"