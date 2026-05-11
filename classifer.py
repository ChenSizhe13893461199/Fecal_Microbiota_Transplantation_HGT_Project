#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function: 
- Automatically finds original BLAST files (by replacing '_gt_' with '_blast_')
- For "other" contigs, generates two separate info files using original BLAST data
- Preserves original outputs for donor/recipient/both categories
"""

import argparse
import os
from Bio import SeqIO

def read_blast_records(blast_file):
    """Read BLAST output file, return dict: query -> list of full lines"""
    records = {}
    if not os.path.exists(blast_file):
        return records
    with open(blast_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            query = line.split('\t')[0]
            records.setdefault(query, []).append(line)
    return records

def write_fasta_and_info(fasta_path, info_path, seq_records, blast_dict):
    """Write FASTA file and corresponding info file (blast lines for each query)."""
    with open(fasta_path, 'w') as f:
        SeqIO.write(seq_records, f, 'fasta')

    info_lines = []
    for rec in seq_records:
        qid = rec.id
        if qid in blast_dict:
            info_lines.extend(blast_dict[qid])

    with open(info_path, 'w') as f:
        for line in info_lines:
            f.write(line + '\n')

def main(contigs_fasta, donor_blast_filtered, recipient_blast_filtered,
         donor_fasta, recipient_fasta, both_fasta, n):
    # n is ignored (kept for compatibility)

    # Locate original BLAST files by removing '_gt_' from filtered paths
    def get_original_path(filtered_path):
        return filtered_path.replace('_gt_', '_blast_')

    donor_blast_original = get_original_path(donor_blast_filtered)
    recipient_blast_original = get_original_path(recipient_blast_filtered)

    # Read filtered BLAST records (for classification and main info files)
    donor_filtered = read_blast_records(donor_blast_filtered)
    recipient_filtered = read_blast_records(recipient_blast_filtered)

    # Read original BLAST records (for other contigs)
    donor_original = read_blast_records(donor_blast_original)
    recipient_original = read_blast_records(recipient_blast_original)

    donor_queries = set(donor_filtered.keys())
    recipient_queries = set(recipient_filtered.keys())

    # Read all contigs
    all_contigs = {rec.id: rec for rec in SeqIO.parse(contigs_fasta, 'fasta')}
    all_queries = set(all_contigs.keys())

    # Classification
    both_q = donor_queries & recipient_queries
    donor_only_q = donor_queries - both_q
    recipient_only_q = recipient_queries - both_q
    other_q = all_queries - donor_queries - recipient_queries

    # Build record lists
    donor_records = [all_contigs[q] for q in donor_only_q if q in all_contigs]
    recipient_records = [all_contigs[q] for q in recipient_only_q if q in all_contigs]
    both_records = [all_contigs[q] for q in both_q if q in all_contigs]
    other_records = [all_contigs[q] for q in other_q if q in all_contigs]

    # --- Output files ---

    # 1. Main outputs (using filtered BLAST)
    donor_info = os.path.splitext(donor_fasta)[0] + '_info.txt'
    recipient_info = os.path.splitext(recipient_fasta)[0] + '_info.txt'
    both_info = os.path.splitext(both_fasta)[0] + '_info.txt'

    write_fasta_and_info(donor_fasta, donor_info, donor_records, donor_filtered)
    write_fasta_and_info(recipient_fasta, recipient_info, recipient_records, recipient_filtered)

    # Both: merge lines from filtered donor and filtered recipient
    combined_filtered = {}
    for q in both_q:
        combined_filtered[q] = donor_filtered.get(q, []) + recipient_filtered.get(q, [])
    write_fasta_and_info(both_fasta, both_info, both_records, combined_filtered)

    # 2. Other outputs (using original BLAST)
    if other_records:
        # Other FASTA (already generated above, but we need the path)
        base_both = os.path.splitext(both_fasta)[0]
        other_fasta = base_both + '_other.fasta'
        # Write other FASTA (if not already written; but we can write again safely)
        with open(other_fasta, 'w') as f:
            SeqIO.write(other_records, f, 'fasta')

        # Other donor info (from original donor BLAST)
        other_donor_info = base_both + '_other_donor_info.txt'
        # Collect lines for other queries from original donor BLAST
        other_donor_lines = []
        for q in other_q:
            if q in donor_original:
                other_donor_lines.extend(donor_original[q])
        with open(other_donor_info, 'w') as f:
            for line in other_donor_lines:
                f.write(line + '\n')

        # Other recipient info (from original recipient BLAST)
        other_recipient_info = base_both + '_other_recipient_info.txt'
        other_recipient_lines = []
        for q in other_q:
            if q in recipient_original:
                other_recipient_lines.extend(recipient_original[q])
        with open(other_recipient_info, 'w') as f:
            for line in other_recipient_lines:
                f.write(line + '\n')
    else:
        # No other contigs, create empty files or skip
        base_both = os.path.splitext(both_fasta)[0]
        other_fasta = base_both + '_other.fasta'
        # Ensure empty files exist (optional)
        open(other_fasta, 'w').close()
        open(base_both + '_other_donor_info.txt', 'w').close()
        open(base_both + '_other_recipient_info.txt', 'w').close()

    print("Classification finished.")
    print(f"  Donor     : {donor_fasta}  +  {donor_info}")
    print(f"  Recipient : {recipient_fasta}  +  {recipient_info}")
    print(f"  Both      : {both_fasta}  +  {both_info}")
    print(f"  Other     : {other_fasta}")
    print(f"    donor info   : {base_both + '_other_donor_info.txt'}")
    print(f"    recipient info: {base_both + '_other_recipient_info.txt'}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Classify contigs based on filtered BLAST results, with original BLAST for other contigs")
    parser.add_argument('-i', '--input', nargs=3, required=True,
                        help='input: contigs.fasta, donor_filtered.txt, recipient_filtered.txt')
    parser.add_argument('-o', '--output', nargs=3, required=True,
                        help='output: donor_contigs.fasta, recipient_contigs.fasta, both_contigs.fasta')
    parser.add_argument('-s', '--similarity', nargs=1, required=True,
                        help='similarity cut-off (ignored, kept for compatibility)')
    args = parser.parse_args()

    main(args.input[0], args.input[1], args.input[2],
         args.output[0], args.output[1], args.output[2],
         float(args.similarity[0]))
