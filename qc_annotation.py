#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qc_annotation.py

Integrates HGT analysis results by matching each predicted gene with its corresponding BLAST HSP
(from donor_contig1.txt). Generates a complete table containing GC content, gene annotations,
and alignment coordinates.

Usage (run directly, no arguments):
    python qc_annotation.py
"""

import pandas as pd
import os
import sys

# ---------- 解析函数 ----------
def parse_emapper_annotations(filepath):
    """
    解析 eggnog-mapper 的 .annotations 文件。
    返回五个字典：
        - gene_info: {recipient_gene_id: (description, cog_category, max_annot_lvl)}
        - gene_to_region: {recipient_gene_id: region_id}
        - gene_coords: {recipient_gene_id: (rec_start, rec_end)}   # 排序后的坐标
        - gene_contigs: {recipient_gene_id: (rec_contig, don_contig)}
        - region_to_genes: {region_id: list_of_recipient_gene_ids}
    """
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        return {}, {}, {}, {}, {}

    gene_info = {}
    gene_to_region = {}
    gene_coords = {}
    gene_contigs = {}
    region_to_genes = {}

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            query = parts[0]          # e.g. "A|B|recipient:start-end_geneNum"
            max_annot_lvl = parts[5]   # Taxonomic_Species
            cog_category = parts[6]     # Module_Classification
            description = parts[7]      # Gene_Description

            # 解析 query 格式：recipient_contig|donor_contig|recipient:start-end_geneNum
            seg = query.split('|')
            if len(seg) != 3:
                continue
            rec_contig = seg[0]          # e.g. NODE_2754_length_5606_cov_2.954242
            don_contig = seg[1]           # e.g. NODE_5840_length_5382_cov_12.879670
            rest = seg[2]                  # e.g. recipient:3105-5606_1
            if not rest.startswith('recipient:'):
                continue
            coord_gene = rest.replace('recipient:', '')  # e.g. 3105-5606_1
            if '_' not in coord_gene:
                continue
            coord, gene_num = coord_gene.rsplit('_', 1)  # coord = 3105-5606, gene_num = 1
            # 解析坐标
            if '-' not in coord:
                continue
            start_str, end_str = coord.split('-')
            try:
                rec_start = int(start_str)
                rec_end = int(end_str)
            except ValueError:
                continue
            # 确保 rec_start <= rec_end 用于区间判断
            if rec_start > rec_end:
                rec_start, rec_end = rec_end, rec_start

            # 构建 recipient_gene_id
            recipient_gene_id = f"{rec_contig}_{coord}_{gene_num}"
            # 构建 region_id (用于匹配 aligned.fasta)
            region_id = f"{rec_contig}|{don_contig}|recipient:{coord}"

            gene_info[recipient_gene_id] = (description, cog_category, max_annot_lvl)
            gene_to_region[recipient_gene_id] = region_id
            gene_coords[recipient_gene_id] = (rec_start, rec_end)
            gene_contigs[recipient_gene_id] = (rec_contig, don_contig)
            region_to_genes.setdefault(region_id, []).append(recipient_gene_id)

    return gene_info, gene_to_region, gene_coords, gene_contigs, region_to_genes


def parse_contig1_blast(filepath):
    """
    解析 donor_contig1.txt（BLAST outfmt 6 格式），
    按 (recipient_contig, donor_contig) 分组存储每个 HSP 的详细信息。
    返回字典：
        {(rec_contig, don_contig): list_of_dict}
        每个 dict 包含：
            'qstart': int, 'qend': int,         # 原始顺序
            'sstart': int, 'send': int,          # 原始顺序
            'qmin': int, 'qmax': int,            # 排序后的 query 区间（用于重叠判断）
            'rate': str, 'length': str
    """
    mapping = {}
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        return mapping

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            rec_contig = parts[0]
            don_contig = parts[1]
            rate = parts[2]
            length = parts[3]
            try:
                qstart = int(parts[6])
                qend = int(parts[7])
                sstart = int(parts[8])
                send = int(parts[9])
            except ValueError:
                continue
            # 排序 query 区间用于重叠判断
            qmin = min(qstart, qend)
            qmax = max(qstart, qend)

            key = (rec_contig, don_contig)
            mapping.setdefault(key, []).append({
                'qstart': qstart,
                'qend': qend,
                'sstart': sstart,
                'send': send,
                'qmin': qmin,
                'qmax': qmax,
                'rate': rate,
                'length': length
            })
    return mapping


def parse_fasta(filepath):
    """解析 fasta 文件，返回 {header: sequence}"""
    seqs = {}
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        return seqs
    with open(filepath, 'r', encoding='utf-8') as f:
        cur_id = None
        cur_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if cur_id:
                    seqs[cur_id] = ''.join(cur_seq)
                cur_id = line[1:]   # 去掉 '>'
                cur_seq = []
            else:
                if line:
                    cur_seq.append(line)
        if cur_id:
            seqs[cur_id] = ''.join(cur_seq)
    return seqs


def gc_content(seq):
    """计算 DNA 序列的 GC 含量（百分比）"""
    seq = seq.upper()
    g = seq.count('G')
    c = seq.count('C')
    total = len(seq)
    return (g + c) / total * 100 if total > 0 else 0.0


def intervals_overlap(start1, end1, start2, end2):
    """判断两个闭区间是否重叠（假设 start <= end）"""
    return max(start1, start2) <= min(end1, end2)


# ---------- 主程序 ----------
def main():
    excel_file = "FMT_list.xlsx"
    if not os.path.exists(excel_file):
        print(f"错误：找不到 {excel_file}")
        sys.exit(1)

    df = pd.read_excel(excel_file, engine='openpyxl')
    required_cols = ['Pre-FMT', 'Donor', 'Post-FMT']
    for col in required_cols:
        if col not in df.columns:
            print(f"错误：Excel 文件中缺少列 '{col}'")
            sys.exit(1)

    os.makedirs("result", exist_ok=True)

    for idx, row in df.iterrows():
        pre = str(row['Pre-FMT'])
        donor = str(row['Donor'])
        post = str(row['Post-FMT'])
        print(f"处理样本 {post} (Donor: {donor})...")

        # 定义文件路径
        ann_file = f"HGT/{post}_donor_HGT.emapper.annotations"
        aligned_fasta = f"HGT/{post}_aligned.fasta"
        # 使用 donor 的 contig1 文件（包含 recipient vs donor 的 HSP）
        contig1_file = f"HGT/{donor}_contig1.txt"
        donor_contig_fasta = f"HGT/{donor}_contig.fasta"

        # 检查必需文件
        if not os.path.exists(ann_file):
            print(f"  警告：注释文件 {ann_file} 不存在，跳过该样本")
            continue
        if not os.path.exists(aligned_fasta):
            print(f"  警告：HGT 区域序列文件 {aligned_fasta} 不存在，跳过该样本")
            continue
        if not os.path.exists(contig1_file):
            print(f"  警告：donor contig1 文件 {contig1_file} 不存在，跳过该样本")
            continue
        if not os.path.exists(donor_contig_fasta):
            print(f"  警告：donor contig fasta 文件 {donor_contig_fasta} 不存在，跳过该样本")
            continue

        # 1. 解析注释
        gene_info, gene_to_region, gene_coords, gene_contigs, region_to_genes = parse_emapper_annotations(ann_file)
        if not gene_info:
            print(f"  警告：注释文件解析失败或无有效数据，跳过该样本")
            continue

        # 2. 解析 donor contig1 BLAST 文件
        blast_hsps = parse_contig1_blast(contig1_file)
        if not blast_hsps:
            print(f"  警告：donor contig1 文件解析失败或无有效数据，跳过该样本")
            continue

        # 3. 解析 aligned.fasta
        region_seqs = parse_fasta(aligned_fasta)
        region_gc = {rid: gc_content(seq) for rid, seq in region_seqs.items()}

        # 4. 解析 donor_contig.fasta
        donor_seqs = parse_fasta(donor_contig_fasta)
        donor_gc = {contig: gc_content(seq) for contig, seq in donor_seqs.items()}

        # 5. 为每个基因组装输出行
        output_rows = []
        for gene_id, (desc, cog_cat, max_lvl) in gene_info.items():
            rec_contig, don_contig = gene_contigs.get(gene_id, ("", ""))
            if not rec_contig or not don_contig:
                continue

            rec_start, rec_end = gene_coords.get(gene_id, (0, 0))
            if rec_start == 0 and rec_end == 0:
                continue

            region_id = gene_to_region.get(gene_id, "")
            gc_hgt_val = region_gc.get(region_id, 0.0)
            gc_origin_val = donor_gc.get(don_contig, 0.0)

            # 在 donor contig1 的 HSP 中查找与基因重叠的 HSP
            key = (rec_contig, don_contig)
            hsps = blast_hsps.get(key, [])
            matched_hsp = None
            for hsp in hsps:
                if intervals_overlap(rec_start, rec_end, hsp['qmin'], hsp['qmax']):
                    matched_hsp = hsp
                    break

            if matched_hsp is None:
                print(f"    警告：基因 {gene_id} 在 {contig1_file} 中找不到与之重叠的 HSP，跳过")
                continue

            # 构建 donor 侧坐标字符串（保留原始顺序）
            donor_coord_str = f"{matched_hsp['sstart']}-{matched_hsp['send']}"
            donor_contig_full = f"{don_contig}_{donor_coord_str}"

            row = {
                'GC_HGT': f"{gc_hgt_val:.4f}",
                'GC_origin': f"{gc_origin_val:.4f}",
                'Gene_Description': desc,
                'Module_Classification': cog_cat,
                'Taxonomic_Species': max_lvl,
                'Recipient_Contig': gene_id,
                'Donor_Contig': donor_contig_full,
                'rate': matched_hsp['rate'],
                'length': matched_hsp['length'],
                'recipient_species': "",
                'donor_species': ""
            }
            output_rows.append(row)

        if output_rows:
            out_file = f"result/{post}_HGT_full.txt"
            header = ['GC_HGT', 'GC_origin', 'Gene_Description', 'Module_Classification',
                      'Taxonomic_Species', 'Recipient_Contig', 'Donor_Contig', 'rate',
                      'length', 'recipient_species', 'donor_species']
            with open(out_file, 'w', encoding='utf-8') as f:
                f.write('\t'.join(header) + '\n')
                for r in output_rows:
                    line = '\t'.join(str(r[col]) for col in header)
                    f.write(line + '\n')
            print(f"  已生成结果文件：{out_file}，共 {len(output_rows)} 条记录")
        else:
            print(f"  警告：没有有效的输出行，跳过样本 {post}")

if __name__ == "__main__":
    main()
