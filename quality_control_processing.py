# -*- coding: utf-8 -*-
"""
quality_control_processing.py

Aggregates HGT events across all FMT samples, applies an edge-based
judgment (0 or 1) using the corresponding pre-FMT recipient BLAST mapping,
then filters records by species level (recipient_species must differ from
donor_species; human species excluded).

Output:
    <root_dir>/HGT_event_details.xlsx  (single sheet: "Overall")

Supported directory layouts
---------------------------
Flat layout:
    <root_dir>/
        HGT1_filtered/                       (or HGT1/)
            <sample>_HGT_statistics1.txt
        blast_results/                       (or result/, optional)
            <sample>_blast_recipient1.txt

Nested layout:
    <root_dir>/
        <cohort>/
            HGT1_filtered/                   (or HGT1/)
                <sample>_HGT_statistics1.txt
            blast_results/
                <sample>_blast_recipient1.txt
"""

import os
import re
import traceback
import pandas as pd
from collections import defaultdict


# ---------------------------------------------------------------------------
# Sheet-name helper
# ---------------------------------------------------------------------------
def safe_sheet_name(name, max_len=31):
    """Convert an arbitrary string into a valid Excel sheet name."""
    name = str(name)
    name = re.sub(r'[\[\]:*?/\\]', '_', name)
    name = name.strip().strip("'")
    name = name[:max_len]
    return name or "Sheet"


# ---------------------------------------------------------------------------
# Species-level helpers
# ---------------------------------------------------------------------------
def extract_species_name(full_name):
    """
    Extract the species-level name (genus + species) from a full description.

    Examples:
        "Faecalibacterium prausnitzii SL3/3"                 -> "Faecalibacterium prausnitzii"
        "Bifidobacterium pseudocatenulatum DSM 20438 = JCM"  -> "Bifidobacterium pseudocatenulatum"
        "Ruminococcus sp. ABC"                               -> "Ruminococcus sp."
        "-" or "" or NaN                                     -> ""
    """
    if full_name is None:
        return ""
    if isinstance(full_name, float) and pd.isna(full_name):
        return ""

    s = str(full_name).strip()
    if s in ("", "-"):
        return ""

    parts = s.split()
    if len(parts) >= 2:
        # Keep "Genus sp." together with a trailing identifier if present
        if parts[1].lower() in ("sp.", "sp", "spp.", "spp"):
            base = parts[0] + " " + parts[1]
            if len(parts) > 2:
                base += " " + parts[2]
            return base
        return parts[0] + " " + parts[1]
    return s


def is_human_species(species_name):
    """Return True if the species string looks like human."""
    if species_name is None:
        return False
    if isinstance(species_name, float) and pd.isna(species_name):
        return False

    s = str(species_name).lower().strip()
    if s in ("", "-"):
        return False

    patterns = ("homo sapiens", "homo_sapiens", "h. sapiens", "h.sapiens",
                "human", "homo")
    for p in patterns:
        if p in s:
            return True
    return s.startswith("homo ")


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------
def extract_base_and_length(contig_str, is_recipient=True):
    """
    Extract the contig base identifier (without the trailing gene number)
    and, for recipient contigs, the length of the coordinate interval.

    Input example (recipient):
        NODE_369_length_73589_cov_48.599043_41175-62863_1
    """
    parts = contig_str.rsplit('_', 1)
    if len(parts) == 2 and parts[1].isdigit():
        base = parts[0]
    else:
        base = contig_str

    length = 0
    if is_recipient:
        coord_part = base.rsplit('_', 1)[-1]
        m = re.match(r'(\d+)-(\d+)$', coord_part)
        if m:
            s, e = int(m.group(1)), int(m.group(2))
            length = abs(e - s) + 1
    return base, length


def extract_pure_recipient_base(contig_str):
    """
    Strip both the trailing gene number and the coordinate suffix.

    Example:
        NODE_1550_length_9207_cov_4.840691_1-3934_1
        -> NODE_1550_length_9207_cov_4.840691
    """
    parts = contig_str.rsplit('_', 1)
    base_with_coord = parts[0] if (len(parts) == 2 and parts[1].isdigit()) else contig_str
    return base_with_coord.rsplit('_', 1)[0]


def extract_length_from_base(base_str):
    m = re.search(r'length_(\d+)', base_str)
    return int(m.group(1)) if m else None


def parse_coordinate_pair(coord_str, sep='-'):
    parts = coord_str.split(sep)
    if len(parts) != 2:
        return None, None
    try:
        return int(parts[0]), int(parts[1])
    except ValueError:
        return None, None


# ---------------------------------------------------------------------------
# Judgment
# ---------------------------------------------------------------------------
def calculate_judgment(rec_base_with_coord, pure_rec_base, pre_recipient_str):
    """
    Return 1 for a potential genuine HGT event, 0 for a likely false positive.
    """
    len_rec = extract_length_from_base(pure_rec_base)
    rec_start, rec_end = parse_coordinate_pair(
        rec_base_with_coord.rsplit('_', 1)[-1], sep='-'
    )
    if None in (rec_start, rec_end, len_rec):
        return 1

    if rec_start < rec_end:
        rec_left = (rec_start == 1)
        rec_right = (rec_end == len_rec)
    else:
        rec_left = (rec_start == len_rec)
        rec_right = (rec_end == 1)

    if not rec_left and not rec_right:
        return 1

    if not pre_recipient_str:
        return 1

    pre_parts = pre_recipient_str.rsplit('_', 2)
    if len(pre_parts) != 3:
        return 1
    pre_base = pre_parts[0]
    pre_start, pre_end = parse_coordinate_pair(
        f"{pre_parts[1]}_{pre_parts[2]}", sep='_'
    )
    if pre_start is None or pre_end is None:
        return 1

    len_pre = extract_length_from_base(pre_base)
    if len_pre is None:
        return 1

    pre_left = (pre_start == 1)
    pre_right = (pre_end == len_pre)

    if not pre_left and not pre_right:
        return 1

    if pre_start < pre_end:
        if rec_left and pre_left:
            return 0
        if rec_right and pre_right:
            return 0
        return 1
    else:
        if rec_left and pre_right:
            return 0
        if rec_right and pre_left:
            return 0
        return 1


def get_pre_recipient(blast_file_path, query_pure_base):
    """
    Scan a post-vs-pre BLAST file for a given query contig base name.
    Returns 'subject_s_start_s_end' or '' if not found.
    """
    if not blast_file_path or not os.path.isfile(blast_file_path):
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
                if parts[0].strip() == query_pure_base:
                    return f"{parts[1]}_{parts[8]}_{parts[9]}"
    except Exception as e:
        print(f"  Error reading BLAST file {blast_file_path}: {e}")
    return ""


# ---------------------------------------------------------------------------
# Per-file processing
# ---------------------------------------------------------------------------
def parse_hgt_stat_file(stat_filepath, stat_basename, blast_dir):
    """
    Read one *_HGT_statistics*.txt file, aggregate gene rows into events
    keyed by (recipient_base, donor_base), and attach pre-recipient info.
    """
    event_stats = defaultdict(lambda: {
        'length': 0,
        'count': 0,
        'recipient_species': None,
        'donor_species': None,
        'pure_rec_base': None,
        'pre_recipient': "",
        'rate': None,
    })

    try:
        with open(stat_filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"  Error reading {stat_filepath}: {e}")
        return event_stats

    start_idx = 0
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
        rate = parts[9]
        recipient_species = parts[11]
        donor_species = parts[12]

        rec_base, length = extract_base_and_length(recipient_contig, True)
        don_base, _ = extract_base_and_length(donor_contig, False)
        pure_rec_base = extract_pure_recipient_base(recipient_contig)

        stats = event_stats[(rec_base, don_base)]
        stats['count'] += 1
        if stats['length'] == 0:
            stats['length'] = length
        if stats['recipient_species'] is None:
            stats['recipient_species'] = recipient_species
        if stats['donor_species'] is None:
            stats['donor_species'] = donor_species
        if stats['pure_rec_base'] is None:
            stats['pure_rec_base'] = pure_rec_base
            stats['rate'] = rate

    # Attach Pre_Recipient from the corresponding BLAST file, if available
    blast_filename = stat_basename.replace("_HGT_statistics", "_blast_recipient")
    blast_file_path = os.path.join(blast_dir, blast_filename) if blast_dir else None

    for _, stats in event_stats.items():
        if stats['pure_rec_base'] and blast_file_path:
            stats['pre_recipient'] = get_pre_recipient(
                blast_file_path, stats['pure_rec_base']
            )

    return event_stats


# ---------------------------------------------------------------------------
# Directory discovery
# ---------------------------------------------------------------------------
def find_hgt_stat_dirs(root_dir):
    """Return all directories that directly contain *_HGT_statistics*.txt files."""
    hits = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        dirnames[:] = [d for d in dirnames if not d.startswith('.')]
        if any(f.endswith('.txt') and 'HGT_statistics' in f for f in filenames):
            hits.append(dirpath)
    return hits


def find_blast_dir_for(hgt_stat_dir, root_dir):
    """Try to locate the blast_results/ (or result/) directory."""
    candidates = [
        os.path.join(hgt_stat_dir, "blast_results"),
        os.path.join(hgt_stat_dir, "result"),
        os.path.join(os.path.dirname(hgt_stat_dir), "blast_results"),
        os.path.join(os.path.dirname(hgt_stat_dir), "result"),
        os.path.join(root_dir, "blast_results"),
        os.path.join(root_dir, "result"),
    ]
    for c in candidates:
        if os.path.isdir(c):
            return c
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    root_dir = r"."

    hgt_stat_dirs = find_hgt_stat_dirs(root_dir)
    if not hgt_stat_dirs:
        print(f"Error: no directory containing *_HGT_statistics*.txt found "
              f"under {os.path.abspath(root_dir)}")
        return

    print("Found HGT statistics directories:")
    for d in hgt_stat_dirs:
        print(f"  {d}")

    output_excel = os.path.join(root_dir, "HGT_event_details.xlsx")

    all_records = []
    total_before_species_filter = 0
    total_after_judgment = 0

    for stat_dir in hgt_stat_dirs:
        blast_dir = find_blast_dir_for(stat_dir, root_dir)
        if blast_dir is None:
            print(f"Warning: no blast_results/ found for {stat_dir}; "
                  f"Pre_Recipient will be empty.")
        else:
            print(f"Using BLAST directory: {blast_dir}")

        stat_files = sorted(
            f for f in os.listdir(stat_dir)
            if f.endswith('.txt') and 'HGT_statistics' in f
        )

        for stat_file in stat_files:
            sample = re.sub(r'_HGT_statistics.*\.txt$', '', stat_file)
            stat_path = os.path.join(stat_dir, stat_file)

            print(f"Processing sample {sample} ({stat_file}) ...")
            event_stats = parse_hgt_stat_file(stat_path, stat_file, blast_dir)
            if not event_stats:
                print(f"  No valid event data in {stat_file}")
                continue

            total_before_species_filter += len(event_stats)

            for (rec_base, don_base), stats in event_stats.items():
                # ---- 1) Edge-based judgment ----
                judgment = calculate_judgment(
                    rec_base, stats['pure_rec_base'], stats['pre_recipient']
                )
                if judgment != 1:
                    continue
                total_after_judgment += 1

                # ---- 2) Species-level filter ----
                rec_species_full = stats['recipient_species']
                don_species_full = stats['donor_species']

                rec_species_simple = extract_species_name(rec_species_full)
                don_species_simple = extract_species_name(don_species_full)

                # Drop if recipient == donor at species level
                if rec_species_simple and don_species_simple \
                        and rec_species_simple == don_species_simple:
                    continue

                # Drop if either side is human
                if is_human_species(rec_species_full) or is_human_species(don_species_full):
                    continue

                # Drop if species information is missing / invalid
                if rec_species_simple in ("", "-") or don_species_simple in ("", "-"):
                    continue

                try:
                    rate_val = float(stats['rate']) if stats['rate'] is not None else None
                except ValueError:
                    rate_val = stats['rate']

                all_records.append({
                    "FMT": sample,
                    "Recipient_Base": rec_base,
                    "Pre_Recipient": stats['pre_recipient'],
                    "Donor_Base": don_base,
                    "Judgment": judgment,
                    "Region_Length": stats['length'],
                    "Homologous_Rate": rate_val,
                    "Gene_Count": stats['count'],
                    "Recipient_Species": rec_species_full,
                    "Donor_Species": don_species_full,
                    "Recipient_Species_Simple": rec_species_simple,
                    "Donor_Species_Simple": don_species_simple,
                    "File": stat_file,
                })

    # ---------------- Summary ----------------
    print("\n" + "=" * 60)
    print(f"Total events parsed            : {total_before_species_filter}")
    print(f"Passed Judgment=1              : {total_after_judgment}")
    print(f"Passed species-level filter    : {len(all_records)}")
    print("=" * 60)

    if not all_records:
        print("No events passed all filters; nothing to write.")
        return

    # ---------------- Write Excel (Overall only) ----------------
    column_order = [
        "FMT", "Recipient_Base", "Pre_Recipient", "Donor_Base", "Judgment",
        "Region_Length", "Homologous_Rate", "Gene_Count",
        "Recipient_Species", "Donor_Species",
        "Recipient_Species_Simple", "Donor_Species_Simple",
        "File",
    ]

    overall_df = pd.DataFrame(all_records)
    overall_df.sort_values(
        by=["FMT", "Recipient_Base", "Donor_Base"], inplace=True
    )
    overall_df = overall_df[column_order].reset_index(drop=True)

    try:
        with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
            overall_df.to_excel(writer, sheet_name="Overall", index=False)
            try:
                writer.book.active = 0
                for ws in writer.book.worksheets:
                    ws.sheet_state = 'visible'
            except Exception:
                pass
        print(f"\nAll results saved to: {output_excel} "
              f"(sheet: Overall, {len(overall_df)} rows)")

    except Exception as e:
        print(f"[ERROR] Failed to write Excel ({e}). Falling back to TSV ...")
        traceback.print_exc()
        tsv_path = os.path.join(root_dir, "HGT_event_details_Overall.tsv")
        overall_df.to_csv(tsv_path, index=False, sep='\t')
        print(f"  Wrote {tsv_path}")


if __name__ == "__main__":
    main()
