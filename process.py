#!/usr/bin/env python3
"""
Process a tab-separated file: for each unique value in the 2nd column,
keep only the row with the maximum value in the 4th column.
Usage: python process.py --input <input_file> --output <output_file>
"""

import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Keep rows with max 4th column for each unique 2nd column."
    )
    parser.add_argument("--input", required=True, help="Input file path")
    parser.add_argument("--output", required=True, help="Output file path")
    args = parser.parse_args()

    input_path = args.input
    output_path = args.output

    # Dictionary: key = value of 2nd column, value = (line, max_4th_column)
    best = {}
    # Preserve order of first appearance of each 2nd column value
    order = []

    try:
        with open(input_path, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                if not line.strip():  # skip empty lines
                    continue
                fields = line.split('\t')
                if len(fields) < 4:
                    # Not enough columns, skip or raise warning
                    print(f"Warning: skipping malformed line: {line}", file=sys.stderr)
                    continue
                col2 = fields[1]
                try:
                    col4 = int(fields[3])
                except ValueError:
                    print(f"Warning: non-integer in 4th column: {line}", file=sys.stderr)
                    continue

                if col2 not in best:
                    best[col2] = (line, col4)
                    order.append(col2)
                else:
                    current_best_line, current_best_val = best[col2]
                    if col4 > current_best_val:
                        best[col2] = (line, col4)
    except FileNotFoundError:
        print(f"Error: input file '{input_path}' not found.", file=sys.stderr)
        sys.exit(1)

    # Write results in original order of first appearance
    with open(output_path, 'w') as f:
        for col2 in order:
            line, _ = best[col2]
            f.write(line + '\n')

if __name__ == "__main__":
    main()