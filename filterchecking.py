# -*- coding: utf-8 -*-
"""
filterchecking.py - Final Quality Control Filter for HGT Events

Processes all *_HGT_statistics.txt files in the 'final' directory.
Filters out records where recipient and donor species are identical
or where either species is human (Homo sapiens). Saves the filtered
results to the 'filter' directory while preserving subfolder structure.

Usage: python filterchecking.py
"""

import os
import pandas as pd
from pathlib import Path

def extract_species_name(full_name):
    """
    从完整的物种描述中提取物种名称（属+种）
    例如: "Bifidobacterium pseudocatenulatum DSM 20438 = JCM 1200 = LMG 10505" -> "Bifidobacterium pseudocatenulatum"
    """
    if pd.isna(full_name) or full_name == "-" or full_name == "":
        return ""
    
    parts = str(full_name).strip().split()
    if len(parts) >= 2:
        # Check if the second part indicates an unidentified species ("sp.", "sp", "spp.")
        # Return "Genus sp." (and optionally a third part if present)
        if parts[1] in ["sp.", "sp", "spp.", "spp"]:
            return parts[0] + " " + parts[1] + ((" " + parts[2]) if len(parts) > 2 else "")
        else:
        # Return "Genus species"
            return parts[0] + " " + parts[1]
    else:
        return str(full_name).strip()

def is_human_species(species_name):
    """
    Check whether a species name indicates human (Homo sapiens).

    Args:
        species_name (str or NaN): The species name to check.

    Returns:
        bool: True if the species matches human-related patterns, False otherwise.
    """
    if pd.isna(species_name) or species_name == "-" or species_name == "":
        return False
    
    species_lower = str(species_name).lower().strip()
    
    human_patterns = [
        "homo sapiens",
        "homo_sapiens",
        "h. sapiens",
        "h.sapiens",
        "human",
        "homo"
    ]
    
    for pattern in human_patterns:
        if pattern in species_lower:
            return True
    
    if species_lower.startswith("homo "):
        return True
    
    return False

def filter_hgt_files(input_base_dir, output_base_dir):
    """
    Recursively process all *_HGT_statistics.txt files, apply filters,
    and save the filtered data to the output directory.

    Filters applied:
        1. Recipient and donor species must not be identical (exact and simplified).
        2. Neither species may be human.
        3. Both species must be non-empty and not "-".

    Args:
        input_base_dir (str): Root directory containing the input files (typically 'final').
        output_base_dir (str): Directory where filtered files will be written (typically 'filter').
    """
    input_base_path = Path(input_base_dir)
    output_base_path = Path(output_base_dir)
    output_base_path.mkdir(parents=True, exist_ok=True)
    
    total_files = 0
    processed_files = 0
    filtered_records = 0
    kept_records = 0
    human_records = 0
    
    # Recursively traverse the input directory
    for root, dirs, files in os.walk(input_base_path):
        root_path = Path(root)
        
        for file_name in files:
            # Only process files ending with _HGT_statistics.txt
            if file_name.endswith('_HGT_statistics.txt'):
                input_file_path = root_path / file_name
                total_files += 1
                
                try:
                    # Read the tab-separated file into a DataFrame
                    df = pd.read_csv(input_file_path, sep='\t')
                    
                    if df.empty:
                        print(f"File is empty: {input_file_path}")
                        continue
                    
                    # Verify that required columns exist
                    required_columns = ['recipient_species', 'donor_species']
                    if not all(col in df.columns for col in required_columns):
                        print(f"文件缺少必要列: {input_file_path}")
                        continue
                    
                    # Work on a copy to avoid modifying the original
                    df_filtered = df.copy()
                    
                    df_filtered['recipient_species_simple'] = df_filtered['recipient_species'].apply(extract_species_name)
                    df_filtered['donor_species_simple'] = df_filtered['donor_species'].apply(extract_species_name)
                    
                    # Build filter conditions
                    condition_same_exact = df_filtered['recipient_species'] == df_filtered['donor_species']
                    condition_same_species = df_filtered['recipient_species_simple'] == df_filtered['donor_species_simple']
                    
                    condition_human_recipient = df_filtered['recipient_species'].apply(is_human_species)
                    condition_human_donor = df_filtered['donor_species'].apply(is_human_species)
                    condition_human = condition_human_recipient | condition_human_donor
                    
                    # Exclude empty or invalid species values
                    condition_valid_recipient = ~df_filtered['recipient_species_simple'].isin(["", "-", None])
                    condition_valid_donor = ~df_filtered['donor_species_simple'].isin(["", "-", None])
                    
                    # Apply all filters: keep records where:
                    #   - species are not identical (exact nor simplified)
                    #   - neither is human
                    #   - both have valid species names
                    df_filtered = df_filtered[
                        (~condition_same_exact) & 
                        (~condition_same_species) & 
                        (~condition_human) &
                        condition_valid_recipient & 
                        condition_valid_donor
                    ].copy()
                    
                    # Remove temporary columns used for filtering
                    if 'recipient_species_simple' in df_filtered.columns:
                        df_filtered = df_filtered.drop(['recipient_species_simple', 'donor_species_simple'], axis=1)
                    
                    # Gather statistics
                    original_count = len(df)
                    filtered_count = len(df_filtered)
                    human_count = condition_human.sum()
                    
                    filtered_records += (original_count - filtered_count - human_count)
                    human_records += human_count
                    kept_records += filtered_count
                    
                    # Write output if any records remain
                    if not df_filtered.empty:
                        # Preserve relative path structure from input directory
                        relative_path = input_file_path.relative_to(input_base_path)
                        output_file_path = output_base_path / relative_path
                        
                        # Create parent directories if needed
                        output_file_path.parent.mkdir(parents=True, exist_ok=True)
                        
                        # Save filtered data
                        df_filtered.to_csv(output_file_path, sep='\t', index=False)
                        processed_files += 1
                        print(f"Processed: {input_file_path} -> {output_file_path}")
                        print(f"Original records: {original_count}, Filtered: {filtered_count}, Human records removed: {human_count}")
                    else:
                        print(f"File became empty after filtering, skipping: {input_file_path})
                        if human_count > 0:
                            print(f"Included {human_count} human records")
                        
                except Exception as e:
                    print(f"Error processing file {input_file_path}: {e}"")
                    import traceback
                    traceback.print_exc()
    
    # Print final summary
    print("\n" + "="*50)
    print("Processing complete!")
    print(f"Total files found: {total_files}")
    print(f"Successfully processed (non-empty output): {processed_files}")
    print(f"Total original records: {filtered_records + kept_records + human_records}")
    print(f"Records filtered (identical species): {filtered_records}")
    print(f"Records filtered (human): {human_records}")
    print(f"Records kept: {kept_records}")
    print("="*50)

def main():
    # Define input and output base directories
    input_base_dir = r"final"       # Directory containing raw HGT statistics files
    output_base_dir = r"filter"     # Directory for filtered results
    
    # Run the filtering function
    filter_hgt_files(input_base_dir, output_base_dir)
    
    print(f"\nAll files processed!")
    print(f"Filtered results saved to: {output_base_dir}")

if __name__ == "__main__":
    main()
