#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

def find_lines_and_write_to_file(input_file, output_file, target_value):
    """Find line numbers containing the target value and write to the output file in sorted order."""
    
    # Open the input file and read all lines
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Find line numbers that contain the target value
    line_numbers = sorted(set(i for i, line in enumerate(lines, start=1) if line.strip() == target_value))

    if not line_numbers:
        return  # If no matching lines, skip writing

    # Ensure previous entries in the output file are also considered
    existing_numbers = set()
    if os.path.exists(output_file):
        with open(output_file, 'r') as existing_file:
            existing_numbers = {int(line.strip()) for line in existing_file if line.strip().isdigit()}
    
    # Merge and sort unique values
    all_numbers = sorted(existing_numbers | set(line_numbers))

    # Write the sorted numbers to the output file
    with open(output_file, 'w') as outfile:  # Use 'w' to overwrite with sorted values
        outfile.writelines(f"{num}\n" for num in all_numbers)

    print(f"Found '{target_value}' on lines: {line_numbers} -> Written in sorted order to {output_file}")

# Main function to handle inputs
def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <boundary_file>")
        sys.exit(1)

    # Get the input boundary nodes file from the command line argument
    input_file = sys.argv[1]

    # Define the directory path (assumes files are in the same directory as the input file)
    run_dir = os.path.dirname(input_file)

    # Ask the user which borehole they set up the interval for
    print('\n')        
    borehole_choice = input("Which borehole did you set up the interval for? (TC/TL/TN/TS/TU) [Default: TC]: \n").strip().upper()
    borehole_mapping = {
        "TC": "8",
        "TL": "9",
        "TN": "12",
        "TS": "10",
        "TU": "11"
    }
    interval_value = borehole_mapping.get(borehole_choice, '8')
    # List of target values and corresponding output file names
    targets = [
        (interval_value, 'well_non_interval.vset'),
        ('8', 'TC.vset'),
        ('9', 'TL.vset'),
        ('10', 'TS.vset'),
        ('11', 'TU.vset'),
        ('12', 'TN.vset'),
        ('1', 'top.vset'),
        ('2', 'bottom.vset'),
        ('3', 'north.vset'),
        ('4', 'south.vset'),
        ('5', 'east.vset'),
        ('6', 'west.vset'),
        ('7', 'drift.vset'),
        ('13', 'well_interval.vset'),
        ('14', 'TL_cell.vset'),
        ('15', 'TS_cell.vset'),
        ('16', 'TU_cell.vset'),
        ('17', 'TN_cell.vset'),
        ('18', 'TC_cell.vset'),
        ('101', 'top.vset'),
        ('101', 'east.vset'),
        ('102', 'top.vset'),
        ('102', 'west.vset'),
        ('103', 'top.vset'),
        ('103', 'north.vset'),
        ('104', 'top.vset'),
        ('104', 'south.vset'),
        ('105', 'bottom.vset'),
        ('105', 'east.vset'),
        ('106', 'bottom.vset'),
        ('106', 'west.vset'),
        ('107', 'bottom.vset'),
        ('107', 'north.vset'),
        ('108', 'bottom.vset'),
        ('108', 'south.vset'),
        ('109', 'north.vset'),
        ('109', 'east.vset'),
        ('110', 'north.vset'),
        ('110', 'west.vset'),
        ('111', 'south.vset'),
        ('111', 'east.vset'),
        ('112', 'south.vset'),
        ('112', 'west.vset'),
        ('113', 'top.vset'),
        ('113', 'north.vset'),
        ('113', 'east.vset'),
        ('114', 'top.vset'),
        ('114', 'north.vset'),
        ('114', 'west.vset'),
        ('115', 'top.vset'),
        ('115', 'south.vset'),
        ('115', 'east.vset'),
        ('116', 'top.vset'),
        ('116', 'south.vset'),
        ('116', 'west.vset'),
        ('117', 'bottom.vset'),
        ('117', 'north.vset'),
        ('117', 'east.vset'),
        ('118', 'bottom.vset'),
        ('118', 'north.vset'),
        ('118', 'west.vset'),
        ('119', 'bottom.vset'),
        ('119', 'south.vset'),
        ('119', 'east.vset'),
        ('120', 'bottom.vset'),
        ('120', 'south.vset'),
        ('120', 'west.vset'),
    ]

    # Process each target
    for target_value, output_file in targets:
        output_path = os.path.join(run_dir, output_file)
        find_lines_and_write_to_file(input_file, output_path, target_value)

if __name__ == "__main__":
    main()
