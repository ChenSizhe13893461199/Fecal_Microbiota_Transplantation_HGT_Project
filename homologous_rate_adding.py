# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 11:33:40 2026

@author: mrrec
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 14:09:52 2026
@author: mrrec
"""

import os
import re
import pandas as pd
from collections import defaultdict

def extract_base_and_length(contig_str, is_recipient=True):
    """
    从 contig 字符串中提取基础标识和（如果是受体）区域长度。
    输入格式示例: "NODE_369_length_73589_cov_48.599043_41175-62863_1"
    返回: (base, length)
        base: 去掉末尾 "_数字" 的部分，如 "NODE_369_length_73589_cov_48.599043_41175-62863"
        length: 如果是受体，根据坐标计算长度；供体则返回 0
    """
    parts = contig_str.rsplit('_', 1)
    if len(parts) == 2 and parts[1].isdigit():
        base = parts[0]
    else:
        base = contig_str

    length = 0
    if is_recipient:
        coord_part = base.rsplit('_', 1)[-1]
        match = re.match(r'(\d+)-(\d+)$', coord_part)
        if match:
            start = int(match.group(1))
            end = int(match.group(2))
            length = abs(end - start) + 1
    return base, length

def extract_pure_recipient_base(contig_str):
    """
    从受体 contig 字符串中提取纯 contig 标识（不含坐标和末尾数字）。
    输入示例: "NODE_1550_length_9207_cov_4.840691_1-3934_1"
    返回: "NODE_1550_length_9207_cov_4.840691"
    """
    parts = contig_str.rsplit('_', 1)
    if len(parts) == 2 and parts[1].isdigit():
        base_with_coord = parts[0]
    else:
        base_with_coord = contig_str
    pure_base = base_with_coord.rsplit('_', 1)[0]
    return pure_base

def extract_length_from_base(base_str):
    """从纯 contig 标识中提取全长数字，如 'NODE_1550_length_9207_cov_4.840691' → 9207"""
    match = re.search(r'length_(\d+)', base_str)
    return int(match.group(1)) if match else None

def parse_coordinate_pair(coord_str, sep='-'):
    """
    解析坐标字符串，返回 (start, end) 且保证 start ≤ end。
    支持分隔符 '-' 或 '_'。
    """
    parts = coord_str.split(sep)
    if len(parts) != 2:
        return None, None
    try:
        a, b = int(parts[0]), int(parts[1])
        return (a, b) if a <= b else (b, a)
    except:
        return None, None

def calculate_judgment(rec_base_with_coord, pure_rec_base, pre_recipient_str):
    """
    根据规则计算判定值：
    - 若受体区域不在边缘 → 1
    - 若受体在边缘：
        - 无 Pre_Recipient 或解析失败 → 1
        - Pre_Recipient 不在边缘 → 1
        - Pre_Recipient 在边缘：
            - 同侧 → 0
            - 不同侧 → 1
    """
    # 1. 解析受体全长及坐标
    len_rec = extract_length_from_base(pure_rec_base)
    if len_rec is None:
        return 1  # 无法获取全长，保守返回1

    coord_part_rec = rec_base_with_coord.rsplit('_', 1)[-1]  # 如 "1-3934"
    rec_start, rec_end = parse_coordinate_pair(coord_part_rec, sep='-')
    if rec_start is None or rec_end is None:
        return 1

    rec_left_edge = (rec_start == 1)
    rec_right_edge = (rec_end == len_rec)

    # 受体不在边缘 → 1
    if not rec_left_edge and not rec_right_edge:
        return 1

    # 受体在边缘，检查 Pre_Recipient
    if not pre_recipient_str:
        return 1  # 无信息，保守

    # 2. 解析 Pre_Recipient
    # 格式： "NODE_893_length_7882_cov_32.011499_1563_6837"
    pre_parts = pre_recipient_str.rsplit('_', 2)
    if len(pre_parts) != 3:
        return 1
    pre_base = pre_parts[0]
    pre_coord_str = f"{pre_parts[1]}_{pre_parts[2]}"  # 合并为 "1563_6837"
    pre_start, pre_end = parse_coordinate_pair(pre_coord_str, sep='_')
    if pre_start is None or pre_end is None:
        return 1

    len_pre = extract_length_from_base(pre_base)
    if len_pre is None:
        return 1

    pre_left_edge = (pre_start == 1)
    pre_right_edge = (pre_end == len_pre)

    # Pre_Recipient 不在边缘 → 1
    if not pre_left_edge and not pre_right_edge:
        return 1

    # 两者都在边缘，判断是否同侧
    if rec_left_edge and pre_left_edge:
        return 0
    if rec_right_edge and pre_right_edge:
        return 0
    return 1

def get_pre_recipient(blast_file_path, query_pure_base):
    """
    在 blast_recipient 文件中查找 query_pure_base 对应的参考序列信息。
    返回格式: "subject_s_start_s_end"，如果找不到则返回空字符串。
    """
    if not os.path.isfile(blast_file_path):
        return ""
    try:
        with open(blast_file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 12:
                    continue
                query = parts[0].strip()
                if query == query_pure_base:
                    subject = parts[1]
                    s_start = parts[8]
                    s_end = parts[9]
                    return f"{subject}_{s_start}_{s_end}"
    except Exception as e:
        print(f"读取 blast 文件 {blast_file_path} 时出错: {e}")
    return ""

def process_fmt_hgt1_folder(fmt_path, hgt1_path):
    """
    处理一个 FMT 的 HGT1 文件夹，返回该 FMT 下所有事件的统计信息。
    事件定义为 (recipient_base, donor_base) 对。
    返回一个字典，键为 (rec_base, don_base)，值为：
        {'length': length, 'count': gene_count,
         'recipient_species': species, 'donor_species': species,
         'pure_rec_base': pure_rec_base, 'source_file': source_file,
         'pre_recipient': pre_recipient, 'rate': rate}
    """
    event_stats = defaultdict(lambda: {
        'length': 0,
        'count': 0,
        'recipient_species': None,
        'donor_species': None,
        'pure_rec_base': None,
        'source_file': None,
        'pre_recipient': "",
        'rate': None
    })

    txt_files = [f for f in os.listdir(hgt1_path) if f.endswith('.txt')]
    for txt_file in txt_files:
        filepath = os.path.join(hgt1_path, txt_file)
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                start_idx = 0
                # 跳过可能的表头（如果第一列不是数字）
                if lines and not lines[0].split('\t')[0].replace('.', '').isdigit():
                    start_idx = 1
                for line in lines[start_idx:]:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) < 13:  # 至少13列
                        continue
                    recipient_contig = parts[7]   # 第8列
                    donor_contig = parts[8]       # 第9列
                    recipient_species = parts[11] # 第12列
                    donor_species = parts[12]     # 第13列
                    rate = parts[9]               # 第10列，rate

                    rec_base, length = extract_base_and_length(recipient_contig, is_recipient=True)
                    don_base, _ = extract_base_and_length(donor_contig, is_recipient=False)
                    pure_rec_base = extract_pure_recipient_base(recipient_contig)

                    key = (rec_base, don_base)
                    stats = event_stats[key]
                    stats['count'] += 1

                    if stats['length'] == 0:
                        stats['length'] = length
                    elif stats['length'] != length:
                        print(f"警告: 事件 {key} 的长度不一致，之前记录为 {stats['length']}，当前为 {length}，将保留第一次的值。")

                    if stats['recipient_species'] is None:
                        stats['recipient_species'] = recipient_species
                    elif stats['recipient_species'] != recipient_species:
                        print(f"警告: 事件 {key} 的受体物种不一致，之前记录为 {stats['recipient_species']}，当前为 {recipient_species}，将保留第一次的值。")

                    if stats['donor_species'] is None:
                        stats['donor_species'] = donor_species
                    elif stats['donor_species'] != donor_species:
                        print(f"警告: 事件 {key} 的供体物种不一致，之前记录为 {stats['donor_species']}，当前为 {donor_species}，将保留第一次的值。")

                    if stats['source_file'] is None:
                        stats['source_file'] = txt_file
                        stats['pure_rec_base'] = pure_rec_base
                        stats['rate'] = rate   # 记录第一次出现的rate
                    else:
                        # 检查rate是否一致，不一致时警告并保留第一次的值
                        if stats['rate'] != rate:
                            print(f"警告: 事件 {key} 的rate不一致，之前记录为 {stats['rate']}，当前为 {rate}，将保留第一次的值。")
        except Exception as e:
            print(f"处理文件 {filepath} 时出错: {e}")

    # 获取 Pre-Recipient
    blast_dir = os.path.join(fmt_path, "blast_results")
    if not os.path.isdir(blast_dir):
        print(f"警告: {fmt_path} 中不存在 blast_results 文件夹，Pre-Recipient 将为空")
        blast_dir = None

    for key, stats in event_stats.items():
        pure_rec_base = stats['pure_rec_base']
        source_file = stats['source_file']
        if source_file is None or pure_rec_base is None:
            continue
        blast_filename = source_file.replace("_HGT_statistics", "_blast_recipient")
        if blast_dir:
            blast_file_path = os.path.join(blast_dir, blast_filename)
            pre_recipient = get_pre_recipient(blast_file_path, pure_rec_base)
            stats['pre_recipient'] = pre_recipient
        else:
            stats['pre_recipient'] = ""

    return event_stats

def main():
    root_dir = r"D:\最后可以提交的算法文件"
    output_excel = os.path.join(root_dir, "HGT_event_details111111.xlsx")

    fmt_results = {}
    all_records = []

    for fmt_name in os.listdir(root_dir):
        fmt_path = os.path.join(root_dir, fmt_name)
        if not os.path.isdir(fmt_path):
            continue

        hgt1_path = os.path.join(fmt_path, "HGT1")
        if not os.path.isdir(hgt1_path):
            continue

        print(f"正在处理 {fmt_name} ...")
        event_stats = process_fmt_hgt1_folder(fmt_path, hgt1_path)

        if not event_stats:
            print(f"{fmt_name} 中没有有效事件数据")
            continue

        records = []
        for (rec_base, don_base), stats in event_stats.items():
            judgment = calculate_judgment(rec_base, stats['pure_rec_base'], stats['pre_recipient'])
            if judgment == 1:  # 只保留 Judgment=1 的事件
                # 将 rate 转换为数值，保留4位小数（若为字符串则直接保留）
                try:
                    rate_val = float(stats['rate']) if stats['rate'] is not None else None
                except:
                    rate_val = stats['rate']
                record = {
                    "Recipient_Base": rec_base,
                    "Pre_Recipient": stats['pre_recipient'],
                    "Donor_Base": don_base,
                    "Judgment": judgment,
                    "Region_Length": stats['length'],
                    "Homologous_Rate": rate_val,          # 新增列
                    "Gene_Count": stats['count'],
                    "Recipient_Species": stats['recipient_species'],
                    "Donor_Species": stats['donor_species'],
                    "File": stats['source_file']
                }
                records.append(record)
                all_records.append({
                    "FMT": fmt_name,
                    "Recipient_Base": rec_base,
                    "Pre_Recipient": stats['pre_recipient'],
                    "Donor_Base": don_base,
                    "Judgment": judgment,
                    "Region_Length": stats['length'],
                    "Homologous_Rate": rate_val,
                    "Gene_Count": stats['count'],
                    "Recipient_Species": stats['recipient_species'],
                    "Donor_Species": stats['donor_species'],
                    "File": stats['source_file']
                })

        records.sort(key=lambda x: (x["Recipient_Base"], x["Donor_Base"]))
        fmt_results[fmt_name] = records

    # 写入 Excel
    with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
        # 调整列顺序，将 Homologous_Rate 放在 Region_Length 之后
        column_order = ["Recipient_Base", "Pre_Recipient", "Donor_Base", "Judgment",
                        "Region_Length", "Homologous_Rate", "Gene_Count",
                        "Recipient_Species", "Donor_Species", "File"]

        for fmt_name, records in fmt_results.items():
            if not records:
                continue
            df = pd.DataFrame(records)
            total_events = len(df)
            total_genes = df['Gene_Count'].sum()
            total_row = pd.DataFrame([{
                "Recipient_Base": "总计",
                "Pre_Recipient": "",
                "Donor_Base": "",
                "Judgment": "",
                "Region_Length": "",
                "Homologous_Rate": "",
                "Gene_Count": f"{total_genes} (共{total_events}个事件)",
                "Recipient_Species": "",
                "Donor_Species": "",
                "File": ""
            }])
            df = pd.concat([df, total_row], ignore_index=True)
            df = df[column_order]
            sheet_name = fmt_name[:31]
            df.to_excel(writer, sheet_name=sheet_name, index=False)
            print(f"已写入 {fmt_name} 共 {total_events} 个事件到工作表 {sheet_name}")

        if all_records:
            overall_df = pd.DataFrame(all_records)
            overall_df.sort_values(by=["FMT", "Recipient_Base", "Donor_Base"], inplace=True)
            overall_columns = ["FMT"] + column_order
            overall_df = overall_df[overall_columns]
            overall_df.to_excel(writer, sheet_name="Overall", index=False)
            print(f"已写入总体统计共 {len(overall_df)} 个事件到工作表 Overall")

    print(f"\n所有结果已保存至: {output_excel}")

if __name__ == "__main__":
    main()