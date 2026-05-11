#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modified version: keeps all alignments meeting:
    alignment_length / min(query_len, subject_len) >= threshold
    pident >= 99.0
    evalue <= 1e-10
"""

import pandas as pd
import argparse
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process Blast Results')
    parser.add_argument('-i', '--input', required=True, help='Input Blast File')
    parser.add_argument('-o', '--output', required=True, help='Output New File')
    parser.add_argument('-t', '--threshold', type=float, default=0.9,
                        help='Threshold for alignment length proportion (default: 0.9)')
    return parser.parse_args()

def extract_length(contig_id):
    """Extract length from contig ID like 'NODE_1_length_683391_cov_19.286475'"""
    parts = contig_id.split('_')
    try:
        idx = parts.index('length')
        if idx + 1 < len(parts):
            return float(parts[idx+1])
    except ValueError:
        pass
    # If format unexpected, return 0 (will cause filtering to skip)
    return 0.0

def main():
    args = parse_arguments()
    threshold = args.threshold

    blast_results = pd.read_csv(args.input, sep="\t", header=None)

    filtered_lines = []
    for index, row in blast_results.iterrows():
        query = str(row[0])
        subject = str(row[1])
        pident = float(row[2])
        length = float(row[3])
        evalue = float(row[10])

        q_len = extract_length(query)
        s_len = extract_length(subject)
        if q_len == 0.0 or s_len == 0.0:
            # Skip if length extraction failed
            continue

        min_len = min(q_len, s_len)
        alignment_ratio = length / min_len

        if alignment_ratio >= threshold and pident >= 99.0 and evalue <= 1e-10:
            filtered_lines.append(row)

    with open(args.output, "w") as f:
        for row in filtered_lines:
            f.write("\t".join(map(str, row)) + "\n")

    print(f"Filtered {len(filtered_lines)} alignments with ratio>={threshold}, pident>=99.0, evalue<=1e-10 -> {args.output}")

if __name__ == "__main__":
    main()