# -*- coding: utf-8 -*-

import os
import re
import pandas as pd
from collections import defaultdict

def extract_base_and_length(contig_str, is_recipient=True):
    """
    Extract base identifier and, for recipient contigs, the region length.

    Input format example (recipient): "e.g. NODE_369_length_73589_cov_48.599043_41175-62863_1"
    Returns:
        base (str): The contig ID without the trailing gene number (e.g., "..._41175-62863")
        length (int): For recipient, the length of the coordinate interval; for donor, 0.
    """
    # Split off the trailing "_<gene_num>" for recipient or "_<coord>" for donor
    parts = contig_str.rsplit('_', 1)
    if len(parts) == 2 and parts[1].isdigit():
        base = parts[0]
    else:
        base = contig_str

    length = 0
    if is_recipient:
        # The last part after the last underscore should be the coordinate range, e.g., "41175-62863"
        coord_part = base.rsplit('_', 1)[-1]
        match = re.match(r'(\d+)-(\d+)$', coord_part)
        if match:
            start = int(match.group(1))
            end = int(match.group(2))
            length = abs(end - start) + 1
    return base, length

def extract_pure_recipient_base(contig_str):
    """
    Extract the pure contig base name (without coordinates and final gene number).

    Example: "NODE_1550_length_9207_cov_4.840691_1-3934_1" -> "NODE_1550_length_9207_cov_4.840691"
    """
    parts = contig_str.rsplit('_', 1)
    if len(parts) == 2 and parts[1].isdigit():
        base_with_coord = parts[0]
    else:
        base_with_coord = contig_str
    pure_base = base_with_coord.rsplit('_', 1)[0]
    return pure_base

def extract_length_from_base(base_str):
     """
    Extract the full contig length from a base identifier string.

    Example: "NODE_1550_length_9207_cov_4.840691" -> 9207
    Returns int or None if not found.
    """
    match = re.search(r'length_(\d+)', base_str)
    return int(match.group(1)) if match else None

def parse_coordinate_pair(coord_str, sep='-'):
    """
    Parse a coordinate string into a tuple (start, end) with start <= end.
    Supports '-' or '_' as separator.

    Returns (None, None) if parsing fails.
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
    Calculate a judgment value (0 or 1) based on edge positions.

    Rules:
        - If recipient region is not at an edge of its contig -> 1 (keep)
        - If recipient is at an edge:
            - Pre_Recipient not at an edge -> 1
            - Both at edges:
                - Same side (both left or both right) -> 0 (likely false positive, discard)
                - Different sides -> 1

    Judgment = 1 indicates a potential genuine HGT event.
    Judgment = 0 indicates a likely false positive (conserved edge region present in pre‑FMT).
    """
    # Parse recipient full length and coordinates
    len_rec = extract_length_from_base(pure_rec_base)
    if len_rec is None:
        return 1  # Cannot determine length, conservative keep

    coord_part_rec = rec_base_with_coord.rsplit('_', 1)[-1]  # 如 "1-3934"
    rec_start, rec_end = parse_coordinate_pair(coord_part_rec, sep='-')
    if rec_start is None or rec_end is None:
        return 1

    rec_left_edge = (rec_start == 1)
    rec_right_edge = (rec_end == len_rec)

    # Recipient not at edge -> keep
    if not rec_left_edge and not rec_right_edge:
        return 1

    # Recipient at edge, need pre‑recipient info
    if not pre_recipient_str:
        return 1  # No info, conservative keep

    # Parse Pre_Recipient string, format: "subject_s_start_s_end" e.g., "NODE_893_length_7882_cov_32.011499_1563_6837"
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

    # Pre‑recipient not at edge -> keep
    if not pre_left_edge and not pre_right_edge:
        return 1

    # Both at edges: check if same side
    if rec_left_edge and pre_left_edge:
        return 0         # Both left edge → likely false positive
    if rec_right_edge and pre_right_edge:
        return 0         # Both right edge → likely false positive
    return 1             # Different edges → keep

def get_pre_recipient(blast_file_path, query_pure_base):
    """
    Search a recipient BLAST file (post vs pre) for the given query contig base name.
    Return a string formatted as "subject_s_start_s_end", or empty string if not found.
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
        print(f"Error reading BLAST file {blast_file_path}: {e}")
    return ""

def process_fmt_hgt1_folder(fmt_path, hgt1_path):
    """
    Process a single FMT sample folder containing HGT1/ subdirectory.
    Aggregates individual gene records into HGT events defined by (recipient_base, donor_base).

    Returns:
        dict: event_stats, key = (rec_base, don_base), value = dict with fields:
            'length', 'count', 'recipient_species', 'donor_species',
            'pure_rec_base', 'source_file', 'pre_recipient', 'rate'
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
    # Find all .txt files in HGT1 directory (each corresponds to a processed HGT statistics file)
    txt_files = [f for f in os.listdir(hgt1_path) if f.endswith('.txt')]
    for txt_file in txt_files:
        filepath = os.path.join(hgt1_path, txt_file)
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                start_idx = 0
                # Skip header line if present (first column not numeric)
                if lines and not lines[0].split('\t')[0].replace('.', '').isdigit():
                    start_idx = 1
                for line in lines[start_idx:]:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) < 13:  
                        continue
                    recipient_contig = parts[7]  
                    donor_contig = parts[8]       
                    recipient_species = parts[11] 
                    donor_species = parts[12]     
                    rate = parts[9]              

                    rec_base, length = extract_base_and_length(recipient_contig, is_recipient=True)
                    don_base, _ = extract_base_and_length(donor_contig, is_recipient=False)
                    pure_rec_base = extract_pure_recipient_base(recipient_contig)

                    key = (rec_base, don_base)
                    stats = event_stats[key]
                    stats['count'] += 1

                    if stats['length'] == 0:
                        stats['length'] = length
                    elif stats['length'] != length:
                        print(f"Warning: event {key} length mismatch")

                    if stats['recipient_species'] is None:
                        stats['recipient_species'] = recipient_species
                    elif stats['recipient_species'] != recipient_species:
                        print(f"Warning: event {key} recipient species mismatch")

                    if stats['donor_species'] is None:
                        stats['donor_species'] = donor_species
                    elif stats['donor_species'] != donor_species:
                        print(f"Warning: event {key} donor species mismatch")

                    if stats['source_file'] is None:
                        stats['source_file'] = txt_file
                        stats['pure_rec_base'] = pure_rec_base
                        stats['rate'] = rate   
                    else:
                        if stats['rate'] != rate:
                            print(f"Warning: event {key} rate mismatch: {stats['rate']} vs {rate}")
        except Exception as e:
            print(f"Error processing file {filepath}: {e}")

    # Retrieve Pre‑Recipient information from original BLAST results (post vs pre)
    blast_dir = os.path.join(fmt_path, "blast_results")
    if not os.path.isdir(blast_dir):
        print(f"Warning: {fmt_path} has no blast_results folder; Pre‑Recipient will be empty")
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
    # Root directory containing individual FMT sample folders (each with HGT1/ and blast_results/)
    root_dir = r"your directory" #please replace it with your directory containing all cohorts directories, in which the directory of HGT1 in each cohort contains the relevant results .txt files
    output_excel = os.path.join(root_dir, "HGT_event_details.xlsx")

    fmt_results = {}
    all_records = []

    for fmt_name in os.listdir(root_dir):
        fmt_path = os.path.join(root_dir, fmt_name)
        if not os.path.isdir(fmt_path):
            continue

        hgt1_path = os.path.join(fmt_path, "HGT1")
        if not os.path.isdir(hgt1_path):
            continue

        print(f"Processing {fmt_name} ...")
        event_stats = process_fmt_hgt1_folder(fmt_path, hgt1_path)

        if not event_stats:
            print(f"No valid event data in {fmt_name}")
            continue

        records = []
        for (rec_base, don_base), stats in event_stats.items():
            judgment = calculate_judgment(rec_base, stats['pure_rec_base'], stats['pre_recipient'])
            if judgment == 1:  # Only keep events with Judgment == 1 (likely true HGT)
                # Convert rate to float if possible
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
                    "Homologous_Rate": rate_val,          
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

    # Write results to an Excel workbook with one sheet per FMT sample and an Overall sheet
    with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
        # Define column order for the output
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
            print(f"Written {fmt_name}: {total_events} events to sheet {sheet_name}")

        if all_records:
            overall_df = pd.DataFrame(all_records)
            overall_df.sort_values(by=["FMT", "Recipient_Base", "Donor_Base"], inplace=True)
            overall_columns = ["FMT"] + column_order
            overall_df = overall_df[overall_columns]
            overall_df.to_excel(writer, sheet_name="Overall", index=False)
            print(f"Written overall statistics: {len(overall_df)} events to sheet Overall")

    print(f"\nAll results saved to: {output_excel}")

if __name__ == "__main__":
    main()
