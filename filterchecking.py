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
import re
from pathlib import Path
import shutil

def extract_species_name(full_name):
    """
    Extract species name (genus + species) from the full species description.
    Example: "Bifidobacterium pseudocatenulatum DSM 20438 = JCM 1200 = LMG 10505" -> "Bifidobacterium pseudocatenulatum"
    """
    if pd.isna(full_name) or full_name == "-" or full_name == "":
        return ""

    # Split the string and take the first two words as species name (genus + species)
    parts = str(full_name).strip().split()
    if len(parts) >= 2:
        # Check if the second part is "sp." or similar, indicating unidentified species
        if parts[1] in ["sp.", "sp", "spp.", "spp"]:
            return parts[0] + " " + parts[1] + ((" " + parts[2]) if len(parts) > 2 else "")
        else:
            return parts[0] + " " + parts[1]
    else:
        return str(full_name).strip()

def is_human_species(species_name):
    """
    Check if the species name contains human (Homo sapiens).
    Supports multiple possible representations.
    """
    if pd.isna(species_name) or species_name == "-" or species_name == "":
        return False

    species_lower = str(species_name).lower().strip()

    # Check various possible representations of human species
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

    # Check if the name starts with "Homo "
    if species_lower.startswith("homo "):
        return True

    return False

def filter_hgt_files(input_base_dir, output_base_dir):
    """
    Traverse the input directory, process all .txt files (HGT results),
    filter and save to the output directory.

    Args:
        input_base_dir: Input root directory, e.g., "HGT1"
        output_base_dir: Output root directory, e.g., "HGT1_filtered"
    """
    input_base_path = Path(input_base_dir)
    output_base_path = Path(output_base_dir)

    # Ensure the output directory exists
    output_base_path.mkdir(parents=True, exist_ok=True)

    # Statistics
    total_files = 0
    processed_files = 0
    filtered_records = 0
    kept_records = 0
    human_records = 0

    # Recursively find all .txt files in the input directory (including subdirectories)
    for txt_file in input_base_path.rglob("*.txt"):
        total_files += 1
        input_file_path = txt_file

        try:
            # Read the file (assuming tab-separated)
            df = pd.read_csv(input_file_path, sep='\t')

            if df.empty:
                print(f"File is empty: {input_file_path}")
                continue

            # Check if required columns exist
            required_columns = ['recipient_species', 'donor_species']
            if not all(col in df.columns for col in required_columns):
                print(f"File missing required columns: {input_file_path}")
                continue

            # Create a copy for filtering
            df_filtered = df.copy()

            # Extract species name (to species level)
            df_filtered['recipient_species_simple'] = df_filtered['recipient_species'].apply(extract_species_name)
            df_filtered['donor_species_simple'] = df_filtered['donor_species'].apply(extract_species_name)

            # Create filtering conditions
            condition_same_exact = df_filtered['recipient_species'] == df_filtered['donor_species']
            condition_same_species = df_filtered['recipient_species_simple'] == df_filtered['donor_species_simple']

            # Check for human species
            condition_human_recipient = df_filtered['recipient_species'].apply(is_human_species)
            condition_human_donor = df_filtered['donor_species'].apply(is_human_species)
            condition_human = condition_human_recipient | condition_human_donor

            # Exclude empty or "-" values
            condition_valid_recipient = ~df_filtered['recipient_species_simple'].isin(["", "-", None])
            condition_valid_donor = ~df_filtered['donor_species_simple'].isin(["", "-", None])

            # Apply filters: keep records with different species and not human
            df_filtered = df_filtered[
                (~condition_same_exact) &
                (~condition_same_species) &
                (~condition_human) &
                condition_valid_recipient &
                condition_valid_donor
            ].copy()

            # Remove auxiliary columns
            if 'recipient_species_simple' in df_filtered.columns:
                df_filtered = df_filtered.drop(['recipient_species_simple', 'donor_species_simple'], axis=1)

            # Statistics
            original_count = len(df)
            filtered_count = len(df_filtered)
            human_count = condition_human.sum()

            filtered_records += (original_count - filtered_count - human_count)
            human_records += human_count
            kept_records += filtered_count

            # Save filtered data if not empty
            if not df_filtered.empty:
                # Build output path while preserving the directory structure
                relative_path = input_file_path.relative_to(input_base_path)
                output_file_path = output_base_path / relative_path

                # Ensure the output directory exists
                output_file_path.parent.mkdir(parents=True, exist_ok=True)

                # Save filtered data
                df_filtered.to_csv(output_file_path, sep='\t', index=False)
                processed_files += 1
                print(f"Processed: {input_file_path} -> {output_file_path}")
                print(f"  Original records: {original_count}, after filter: {filtered_count}, human records: {human_count}")
            else:
                print(f"File empty after filtering, skipped: {input_file_path}")
                if human_count > 0:
                    print(f"  Including human records: {human_count}")

        except Exception as e:
            print(f"Error processing file {input_file_path}: {e}")
            import traceback
            traceback.print_exc()

    # Print summary
    print("\n" + "=" * 50)
    print("Processing completed!")
    print(f"Total files: {total_files}")
    print(f"Successfully processed files: {processed_files}")
    print(f"Total original records: {filtered_records + kept_records + human_records}")
    print(f"Filtered out records (same species): {filtered_records}")
    print(f"Filtered out records (human): {human_records}")
    print(f"Kept records: {kept_records}")
    print("=" * 50)

def main():
    # Set paths
    input_base_dir = "HGT1"
    output_base_dir = "HGT1_filtered"

    # Execute filtering
    filter_hgt_files(input_base_dir, output_base_dir)

    print(f"\nAll files processed!")
    print(f"Filtered results saved in: {output_base_dir}")

if __name__ == "__main__":
    main()
