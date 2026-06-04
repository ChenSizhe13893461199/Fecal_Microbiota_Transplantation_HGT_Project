#!/usr/bin/env python3
"""
HGT.py – Extract contigs that align to both recipient and donor with specific coverage and identity thresholds.

This script identifies candidate horizontal gene transfer (HGT) contigs from an "other" category of contigs
(those not uniquely assigned to donor or recipient) by requiring:
    1. High-identity alignment to the recipient (pre-FMT) covering at least a user-defined fraction of the contig.
    2. A high-quality alignment to the donor that does NOT overlap the recipient-aligned regions.
It outputs full contig sequences, aligned regions, and BLAST records for downstream annotation.
"""
import argparse
import os
from collections import defaultdict
from Bio import SeqIO

def parse_blast(blast_file, min_pident=None):
    hits = defaultdict(list)
    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            qseqid = parts[0]               # query sequence ID
            sseqid = parts[1]               # subject (database) sequence ID
            pident = float(parts[2])        # percent identity
            length = int(parts[3])          # alignment length
            qstart = int(parts[6])          # query start (1-based)
            qend = int(parts[7])            # query end
            sstart = int(parts[8])          # subject start
            send = int(parts[9])            # subject end
            evalue = float(parts[10])       # e-value
            if min_pident is not None and pident < min_pident:
                continue
            hits[qseqid].append((sseqid, pident, length, qstart, qend, sstart, send, evalue, line))
    return hits

def merge_intervals(intervals):
    """
    Merge list of (start, end) intervals (1-based inclusive).
    Returns list of merged intervals.
    """
    if not intervals:
        return []
    intervals.sort()
    merged = [list(intervals[0])]              # start with the first interval
    for start, end in intervals[1:]:
                                               # If current interval starts <= previous end+1, they overlap or touch → merge
        if start <= merged[-1][1] + 1:  # allow touching
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]

def intervals_overlap(interval, other_intervals):
    """Check if interval overlaps with any in other_intervals."""
    s, e = interval
    for os, oe in other_intervals:
        if max(s, os) <= min(e, oe):
            return True
    return False

def main():
    parser = argparse.ArgumentParser(description='Extract HGT candidates from secondary BLAST results.')
    parser.add_argument('--other_fasta', required=True, help='FASTA file of other contigs (e.g., donor_post_contigs_other.fasta)')
    parser.add_argument('--recipient_blast', required=True, help='BLAST output against recipient (e.g., XXX_blast_recipient1.txt)')
    parser.add_argument('--donor_blast', required=True, help='BLAST output against donor (e.g., XXX_blast_donor1.txt)')
    parser.add_argument('--donor_fasta', required=True, help='Donor original FASTA file (e.g., XXX_filter.fasta)')
    parser.add_argument('--output_prefix', nargs=2, required=True, help='Post sample name and donor sample name, e.g., XXX XXX')
    parser.add_argument('--min_recipient_cov', type=float, default=0.5, 
                        help='Minimum coverage of recipient-homologous regions on the contig (fraction, default=0.5)')
    args = parser.parse_args()

    post, donor = args.output_prefix
    min_cov = args.min_recipient_cov

    # Output file names
    # Define output file names (all will be written in the current working directory)
    post_contig_fasta = f"{post}_contig.fasta"        # Full sequences of selected recipient contigs
    post_aligned_fasta = f"{post}_aligned.fasta"      # Recipient-side HGT regions only
    post_info_txt = f"{post}_contig1.txt"             # BLAST lines (recipient hits) for selected contigs
    donor_contig_fasta = f"{donor}_contig.fasta"      # Full donor contigs that participated in HGT
    donor_aligned_fasta = f"{donor}_aligned.fasta"    # Donor-side HGT regions only
    donor_info_txt = f"{donor}_contig1.txt"           # BLAST lines (donor hits) for selected contigs

    # Read other contigs
    contigs = {rec.id: rec for rec in SeqIO.parse(args.other_fasta, 'fasta')}
    # Read donor sequences
    donor_seqs = {rec.id: rec for rec in SeqIO.parse(args.donor_fasta, 'fasta')}

    # Parse BLAST hits
    recip_hits = parse_blast(args.recipient_blast, min_pident=99.0)
    donor_hits = parse_blast(args.donor_blast, min_pident=99.0)

    # Prepare output collections
    selected_queries = set()
    post_aligned_records = []        # for post_aligned.fasta (query regions)
    donor_aligned_records = []       # for donor_aligned.fasta (subject regions)
    donor_contig_ids = set()          # unique donor subjects to output complete sequences
    post_info_lines = []              # for post_info.txt
    donor_info_lines = []             # for donor_info.txt
    
    # Main loop: evaluate each "other" contig
    for qid, contig in contigs.items():
        # ----- Step 1: Must have at least one hit to recipient -----
        if qid not in recip_hits:
            continue

        # ----- Step 2: Collect recipient intervals with high homologous identity (≥98.0%) -----
        recip_intervals = []
        for hit in recip_hits[qid]:
            _, pident, _, qs, qe, _, _, evalue, line = hit
            # pident already filtered, but double-check
            if pident >= 98.0 and evalue <= 1e-10:                 # only consider nearly identical matches
                recip_intervals.append((qs, qe))
                post_info_lines.append(line)   # store line for later output (only if contig is selected)
        if not recip_intervals:
            continue
        
        # Step 3: Merge overlapping/touching intervals and compute coverage
        merged_recip = merge_intervals(recip_intervals)
        total_covered = sum(e - s + 1 for s, e in merged_recip)
        contig_len = len(contig.seq)
        coverage = total_covered / contig_len

        # Step 4: Check recipient coverage threshold (user-defined)
        if coverage < min_cov:
            continue
        # It ensures that the recipient-coverage region occupies at least X % of the contig, as subsequent codes searching for potential HGT regions with length (> 500, <（100-X）% of the contig)
        # balancing sensitivity (capturing genuine HGT insertions) and specificity (excluding donor-only contigs)

        
        # ----- Step 5: Must have at least one hit to donor -----
        # ----- Donor check (pident >= 99%, length > 500, evalue <= 1e-10, non-overlap) -----
        if qid not in donor_hits:
            continue
        donor_good = False

        # ----- Step 6: Evaluate donor hits for satisfying HGT criteria -----
        for hit in donor_hits[qid]:
            sseqid, pident, length, qs, qe, ss, se, evalue, line = hit
            if pident >= 99.0 and length > 500 and evalue <= 1e-10:
                # Strict filters for donor side: ≥99% identity, alignment length >500 bp, e-value ≤1e-10
                if not intervals_overlap((qs, qe), merged_recip):
                    # Crucial: the donor-aligned region must NOT overlap any recipient-aligned region
                    donor_good = True
                    # Collect donor info line
                    donor_info_lines.append(line)
                    # Prepare aligned region from donor subject
                    if sseqid in donor_seqs:
                        # Extract subject region (ss to se)
                        sstart, send = sorted((ss, se))
                        subseq = donor_seqs[sseqid].seq[sstart-1:send]
                        desc = f"{qid}|{sseqid}|donor:{sstart}-{send}"
                        donor_aligned_records.append(SeqIO.SeqRecord(subseq, id=desc, description=""))
                        donor_contig_ids.add(sseqid)
                    # Also prepare aligned region from query (for post_aligned)
                    qstart, qend = sorted((qs, qe))
                    qsubseq = contig.seq[qstart-1:qend]
                    desc2 = f"{qid}|{sseqid}|recipient:{qstart}-{qend}"
                    post_aligned_records.append(SeqIO.SeqRecord(qsubseq, id=desc2, description=""))
        # Step 7: If all conditions satisfied, add the contig to the final set
        if donor_good:
            selected_queries.add(qid)

    # ----- Write output files -----
    # 1. Write full sequences of selected recipient contigs
    with open(post_contig_fasta, 'w') as f:
        for qid in selected_queries:
            SeqIO.write(contigs[qid], f, 'fasta')

    # 2. post_aligned.fasta (query aligned regions)
    with open(post_aligned_fasta, 'w') as f:
        SeqIO.write(post_aligned_records, f, 'fasta')

    # 3. post_info.txt (recipient blast lines)
    with open(post_info_txt, 'w') as f:
        f.writelines(post_info_lines)

    # 4. donor_contig.fasta (full donor subject sequences)
    with open(donor_contig_fasta, 'w') as f:
        for sid in donor_contig_ids:
            if sid in donor_seqs:
                SeqIO.write(donor_seqs[sid], f, 'fasta')

    # 5. donor_aligned.fasta (subject aligned regions)
    with open(donor_aligned_fasta, 'w') as f:
        SeqIO.write(donor_aligned_records, f, 'fasta')

    # 6. donor_info.txt (donor blast lines)
    with open(donor_info_txt, 'w') as f:
        f.writelines(donor_info_lines)

    print(f"Done. Selected {len(selected_queries)} contigs (max_HGT_cov={min_cov}).")
    print(f"Output files: {post_contig_fasta}, {post_aligned_fasta}, {post_info_txt}, {donor_contig_fasta}, {donor_aligned_fasta}, {donor_info_txt}")

if __name__ == "__main__":
    main()
