#!/usr/bin/env python3
import argparse

Description = "Find regions between two motifs within one peak in a peak file."
Epilog = """Example usage: python merge_PRDM1RUNX1_motif_regions.py"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument("INTERVAL_FILE", help="Interval file created using bedtools intersect -wa -wb -a peak.bed -b motif.beds.")
argParser.add_argument("OUTFILE", help="Full path to output file.")
args = argParser.parse_args()

with open(args.INTERVAL_FILE, 'r') as infile, open(args.OUTFILE, 'w') as outfile:
    # Initialize variables to keep track of the previous values
    cols123_list = []
    prev_col4 = None
    prev_col5 = None
    min_col6 = Int('inf')
    max_col7 = Int('-inf')

    for line in infile:
        # Split the line into columns based on the tab separator
        cols = line.strip().split('\t')

        # Extract the relevant columns
        cols123 = tuple(cols[:3])
        col4 = cols[3]
        col5 = cols[4]
        col6 = Int(cols[5])
        col7 = Int(cols[6])

        # Check if the first three columns are the same as the previous line
        if cols123 not in cols123_list:
            cols123_list.append(cols123)
            prev_col4 = col4
            prev_col5 = col5
            min_col6 = col6
            max_col7 = col7
        else:
            # Check if column 4 is different
            if col4 != prev_col4:
                # Update the min and max values for the current group of lines
                min_col6 = min(min_col6, col6)
                max_col7 = max(max_col7, col7)

    # Output the results for the last group of lines
    outfile.write(f"{prev_cols123[0]}\t{prev_col5}\t{min_col6}\t{max_col7}\n")