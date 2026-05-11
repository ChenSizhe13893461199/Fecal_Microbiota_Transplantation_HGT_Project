#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
add_gc_nonhgt.py

Adds two columns to existing *_HGT_statistics.txt files:
    - GC_nonHGT      : GC content (%) of the non-HGT regions in the recipient contig
    - GC_nonHGTdonor : GC content (%) of the non-HGT regions in the donor contig

Calculation basis:
    Recipient: full contig sequence from HGT/{post}_contig.fasta,
               HGT region coordinates from HGT/{post}_aligned.fasta.
    Donor: full contig sequence from HGT/{donor}_contig.fasta,
           HGT region coordinates from HGT/{donor}_contig1.txt (BLAST HSPs).

Output files are saved in the HGT1/ directory as {post}_HGT_statistics1.txt.

Usage:
    python add_gc_nonhgt.py
"""

import os
import pandas as pd
from Bio import SeqIO

def parse_aligned_fasta(aligned_file):
    """
    Parse {post}_aligned.fasta to extract HGT regions on each recipient contig,
    merging overlapping intervals.

    Args:
        aligned_file (str): Path to the aligned FASTA file.

    Returns:
        dict: {contig_id: [(start, end), ...]} where coordinates are sorted and non-overlapping.
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

    # Merge overlapping or touching intervals for each contig
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
    Parse donor_contig1.txt (BLAST outfmt 6) to extract HSP coordinates on donor contigs,
    merging overlapping intervals.

    Args:
        contig1_file (str): Path to donor_contig1.txt.

    Returns:
        dict: {donor_contig: [(start, end), ...]} with merged, non-overlapping intervals.
    """
    donor_intervals = {}
    if not os.path.exists(contig1_file) or os.path.getsize(contig1_file) == 0:
        return donor_intervals

    with open(contig1_file, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            don_contig = parts[1]          # Subject (donor) contig ID
            try:
                sstart = int(parts[8])     # Subject start (1-based)
                send = int(parts[9])       # Subject end
            except ValueError:
                continue
            # Ensure start <= end for interval representation
            if sstart > send:
                sstart, send = send, sstart
            donor_intervals.setdefault(don_contig, []).append((sstart, send))

    # Merge overlapping intervals for each donor contig
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
    Calculate GC content (%) of a sequence excluding specified intervals (HGT regions).

    Args:
        seq (Bio.Seq.Seq): The full sequence (recipient or donor contig).
        intervals (list of tuple): List of (start, end) intervals (1‑based inclusive)
                                   corresponding to HGT regions to exclude.

    Returns:
        float: GC percentage of the non-HGT portions. Returns 0.0 if no non‑HGT region exists.
    """
    seq_str = str(seq).upper()
    total_len = len(seq_str)
    if not intervals:
        # No HGT regions → entire sequence is non-HGT
        gc_count = seq_str.count('G') + seq_str.count('C')
        nonhgt_len = total_len
    else:
        gc_count = 0
        pos = 0  # 0‑based position tracking where the previous interval ended
        for s, e in intervals:
            # Region before this interval (if any)
            if s > pos + 1:
                segment = seq_str[pos:s-1]
                gc_count += segment.count('G') + segment.count('C')
            pos = e   # move to the end of this interval (1‑based)
        # Region after the last interval
        if pos < total_len:
            segment = seq_str[pos:]
            gc_count += segment.count('G') + segment.count('C')
        nonhgt_len = total_len - sum(e - s + 1 for s, e in intervals)

    if nonhgt_len == 0:
        return 0.0
    return (gc_count / nonhgt_len) * 100


def process_sample(post, donor, input_dir='.', output_dir='HGT1'):
    """
    Process a single sample pair (post‑FMT recipient and donor).

    Args:
        post (str): Post‑FMT sample name.
        donor (str): Donor sample name.
        input_dir (str): Directory containing the input *_HGT_statistics.txt file.
        output_dir (str): Directory where the enriched output file will be written.
    """
    stat_file = os.path.join(input_dir, f"{post}_HGT_statistics.txt")
    if not os.path.exists(stat_file):
        print(f"Skipping {post}: {stat_file} does not exist")
        return

    # ----- Recipient (post‑FMT) files -----
    rec_contig_fasta = f"HGT/{post}_contig.fasta"
    rec_aligned_fasta = f"HGT/{post}_aligned.fasta"
    if not os.path.exists(rec_contig_fasta) or not os.path.exists(rec_aligned_fasta):
        print(f"Skipping {post}: missing {rec_contig_fasta} or {rec_aligned_fasta}")
        return

    # ----- Donor files -----
    don_contig_fasta = f"HGT/{donor}_contig.fasta"
    don_contig1_file = f"HGT/{donor}_contig1.txt"
    if not os.path.exists(don_contig_fasta) or not os.path.exists(don_contig1_file):
        print(f"Skipping {post}: missing {don_contig_fasta} or {don_contig1_file}")
        return

    # ---------- 1. Compute recipient non‑HGT GC ----------
    rec_intervals = parse_aligned_fasta(rec_aligned_fasta)
    rec_seqs = {rec.id: rec.seq for rec in SeqIO.parse(rec_contig_fasta, 'fasta')}
    rec_nongc = {}
    for contig, seq in rec_seqs.items():
        intervals = rec_intervals.get(contig, [])
        gc = calculate_nonhgt_gc(seq, intervals)
        rec_nongc[contig] = gc

    # ---------- 2. Compute donor non‑HGT GC ----------
    don_intervals = parse_donor_contig1(don_contig1_file)
    don_seqs = {rec.id: rec.seq for rec in SeqIO.parse(don_contig_fasta, 'fasta')}
    don_nongc = {}
    for contig, seq in don_seqs.items():
        intervals = don_intervals.get(contig, [])
        gc = calculate_nonhgt_gc(seq, intervals)
        don_nongc[contig] = gc

    # ---------- 3. Read the original statistics file and add columns ----------
    df = pd.read_csv(stat_file, sep='\t')

    # Extract recipient contig base name (strip "_start-end_genenum")
    def extract_rec_contig(rc):
        parts = rc.rsplit('_', 2)
        return parts[0] if len(parts) == 3 else rc
    df['rec_contig_id'] = df['Recipient_Contig'].apply(extract_rec_contig)
    df['GC_nonHGT'] = df['rec_contig_id'].map(rec_nongc).fillna(0.0)

    # Extract donor contig base name (strip "_start-end")
    def extract_don_contig(dc):
        # 假设格式为 "don_contig_start-end"
        if '_' in dc:
            # Remove last underscore and everything after it
            return dc.rsplit('_', 1)[0]
        return dc
    df['don_contig_id'] = df['Donor_Contig'].apply(extract_don_contig)
    df['GC_nonHGTdonor'] = df['don_contig_id'].map(don_nongc).fillna(0.0)

    # Reorder columns: put the two new GC columns at the front, drop temporary keys
    cols = ['GC_nonHGT', 'GC_nonHGTdonor'] + [c for c in df.columns if c not in ['GC_nonHGT', 'GC_nonHGTdonor', 'rec_contig_id', 'don_contig_id']]
    df = df[cols]

    # ---------- 4. Write output ----------
    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, f"{post}_HGT_statistics1.txt")
    df.to_csv(out_file, sep='\t', index=False)
    print(f"已生成 {out_file}")


def main():
    excel_file = "FMT_list.xlsx"
    if not os.path.exists(excel_file):
        print(f"Error: {excel_file} not found.")
        return

    df_excel = pd.read_excel(excel_file, engine='openpyxl')
    required_cols = ['Pre-FMT', 'Donor', 'Post-FMT']
    for col in required_cols:
        if col not in df_excel.columns:
            print(f"Error: Excel file missing column '{col}'")
            return

    for idx, row in df_excel.iterrows():
        pre = str(row['Pre-FMT'])
        donor = str(row['Donor'])
        post = str(row['Post-FMT'])
        process_sample(post, donor)

if __name__ == "__main__":
    main()
