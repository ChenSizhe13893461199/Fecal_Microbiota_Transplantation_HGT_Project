#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
add_gc_nonhgt.py

为已有的 *_HGT_statistics.txt 文件添加两列：
    - GC_nonHGT      ：受体 contig 中非 HGT 区域的 GC 含量（百分比）
    - GC_nonHGTdonor ：供体 contig 中非 HGT 区域的 GC 含量（百分比）
计算依据：
    受体：从 HGT/{post}_contig.fasta 获得 contig 全长序列，
          从 HGT/{post}_aligned.fasta 获得受体上的 HGT 区域坐标。
    供体：从 HGT/{donor}_contig.fasta 获得 contig 全长序列，
          从 HGT/{donor}_contig1.txt 获得供体上的 HGT 区域坐标（BLAST HSP）。
输出文件保存在 HGT1/ 目录下，文件名为 {post}_HGT_statistics1.txt。

用法：
    python add_gc_nonhgt.py
"""

import os
import pandas as pd
from Bio import SeqIO

def parse_aligned_fasta(aligned_file):
    """
    解析 {post}_aligned.fasta，提取每个受体 contig 上的 HGT 区域（合并重叠区间）。
    返回字典 {contig_id: [(start, end), ...]}，坐标已排序且不重叠。
    """
    contig_intervals = {}
    for rec in SeqIO.parse(aligned_file, 'fasta'):
        parts = rec.id.split('|')
        if len(parts) != 3:
            continue
        qid = parts[0]
        recipient_part = parts[2]
        if not recipient_part.startswith('recipient:'):
            continue
        coord_str = recipient_part.replace('recipient:', '')
        if '-' not in coord_str:
            continue
        start, end = map(int, coord_str.split('-'))
        if start > end:
            start, end = end, start
        contig_intervals.setdefault(qid, []).append((start, end))

    # 合并重叠区间
    for contig, intervals in contig_intervals.items():
        intervals.sort()
        merged = []
        for s, e in intervals:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        contig_intervals[contig] = [(s, e) for s, e in merged]
    return contig_intervals


def parse_donor_contig1(contig1_file):
    """
    解析 donor_contig1.txt（BLAST outfmt 6 格式），
    提取供体 contig 上的 HSP 坐标（sstart/send），合并重叠区间。
    返回字典 {donor_contig: [(start, end), ...]}。
    """
    donor_intervals = {}
    if not os.path.exists(contig1_file) or os.path.getsize(contig1_file) == 0:
        return donor_intervals

    with open(contig1_file, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            don_contig = parts[1]          # subject 是供体 contig
            try:
                sstart = int(parts[8])
                send = int(parts[9])
            except ValueError:
                continue
            # 确保 start <= end
            if sstart > send:
                sstart, send = send, sstart
            donor_intervals.setdefault(don_contig, []).append((sstart, send))

    # 合并每个供体 contig 上的重叠区间
    for contig, intervals in donor_intervals.items():
        intervals.sort()
        merged = []
        for s, e in intervals:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        donor_intervals[contig] = [(s, e) for s, e in merged]
    return donor_intervals


def calculate_nonhgt_gc(seq, intervals):
    """
    计算序列 seq（Bio.Seq 对象）中除 intervals 覆盖区域外的 GC 含量。
    intervals 是已合并的列表，每个元素为 (start, end)，1-based 闭区间。
    返回 GC 百分比（浮点数）。
    """
    seq_str = str(seq).upper()
    total_len = len(seq_str)
    if not intervals:
        gc_count = seq_str.count('G') + seq_str.count('C')
        nonhgt_len = total_len
    else:
        gc_count = 0
        pos = 0  # 0-based 当前位置
        for s, e in intervals:
            if s > pos + 1:
                segment = seq_str[pos:s-1]
                gc_count += segment.count('G') + segment.count('C')
            pos = e  # 更新到区间结束（1-based）
        if pos < total_len:
            segment = seq_str[pos:]
            gc_count += segment.count('G') + segment.count('C')
        nonhgt_len = total_len - sum(e - s + 1 for s, e in intervals)

    if nonhgt_len == 0:
        return 0.0
    return (gc_count / nonhgt_len) * 100


def process_sample(post, donor, input_dir='.', output_dir='HGT1'):
    """处理单个样本 post，donor 为对应的供体样本名"""
    stat_file = os.path.join(input_dir, f"{post}_HGT_statistics.txt")
    if not os.path.exists(stat_file):
        print(f"跳过 {post}：{stat_file} 不存在")
        return

    # 受体相关文件
    rec_contig_fasta = f"HGT/{post}_contig.fasta"
    rec_aligned_fasta = f"HGT/{post}_aligned.fasta"
    if not os.path.exists(rec_contig_fasta) or not os.path.exists(rec_aligned_fasta):
        print(f"跳过 {post}：缺少 {rec_contig_fasta} 或 {rec_aligned_fasta}")
        return

    # 供体相关文件
    don_contig_fasta = f"HGT/{donor}_contig.fasta"
    don_contig1_file = f"HGT/{donor}_contig1.txt"
    if not os.path.exists(don_contig_fasta) or not os.path.exists(don_contig1_file):
        print(f"跳过 {post}：缺少 {don_contig_fasta} 或 {don_contig1_file}")
        return

    # ---------- 1. 计算受体非 HGT GC ----------
    rec_intervals = parse_aligned_fasta(rec_aligned_fasta)
    rec_seqs = {rec.id: rec.seq for rec in SeqIO.parse(rec_contig_fasta, 'fasta')}
    rec_nongc = {}
    for contig, seq in rec_seqs.items():
        intervals = rec_intervals.get(contig, [])
        gc = calculate_nonhgt_gc(seq, intervals)
        rec_nongc[contig] = gc

    # ---------- 2. 计算供体非 HGT GC ----------
    don_intervals = parse_donor_contig1(don_contig1_file)
    don_seqs = {rec.id: rec.seq for rec in SeqIO.parse(don_contig_fasta, 'fasta')}
    don_nongc = {}
    for contig, seq in don_seqs.items():
        intervals = don_intervals.get(contig, [])
        gc = calculate_nonhgt_gc(seq, intervals)
        don_nongc[contig] = gc

    # ---------- 3. 读取原统计文件，添加列 ----------
    df = pd.read_csv(stat_file, sep='\t')

    # 从 Recipient_Contig 提取受体 contig 名称（格式：contig_start-end_genenum）
    def extract_rec_contig(rc):
        parts = rc.rsplit('_', 2)
        return parts[0] if len(parts) == 3 else rc
    df['rec_contig_id'] = df['Recipient_Contig'].apply(extract_rec_contig)
    df['GC_nonHGT'] = df['rec_contig_id'].map(rec_nongc).fillna(0.0)

    # 从 Donor_Contig 提取供体 contig 名称（格式：don_contig_start-end）
    def extract_don_contig(dc):
        # 假设格式为 "don_contig_start-end"
        if '_' in dc:
            # 去掉最后一个下划线及之后的部分
            return dc.rsplit('_', 1)[0]
        return dc
    df['don_contig_id'] = df['Donor_Contig'].apply(extract_don_contig)
    df['GC_nonHGTdonor'] = df['don_contig_id'].map(don_nongc).fillna(0.0)

    # 整理列顺序：把两个 GC 列放在最前面
    cols = ['GC_nonHGT', 'GC_nonHGTdonor'] + [c for c in df.columns if c not in ['GC_nonHGT', 'GC_nonHGTdonor', 'rec_contig_id', 'don_contig_id']]
    df = df[cols]

    # ---------- 4. 输出 ----------
    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, f"{post}_HGT_statistics1.txt")
    df.to_csv(out_file, sep='\t', index=False)
    print(f"已生成 {out_file}")


def main():
    excel_file = "FMT_list.xlsx"
    if not os.path.exists(excel_file):
        print(f"错误：找不到 {excel_file}")
        return

    df_excel = pd.read_excel(excel_file, engine='openpyxl')
    required_cols = ['Pre-FMT', 'Donor', 'Post-FMT']
    for col in required_cols:
        if col not in df_excel.columns:
            print(f"错误：Excel 文件中缺少列 '{col}'")
            return

    for idx, row in df_excel.iterrows():
        pre = str(row['Pre-FMT'])
        donor = str(row['Donor'])
        post = str(row['Post-FMT'])
        process_sample(post, donor)

if __name__ == "__main__":
    main()