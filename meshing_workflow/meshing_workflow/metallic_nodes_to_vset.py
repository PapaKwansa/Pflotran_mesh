#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 09:00:00 2024
@author: yao
"""

import os
import sys

def copy_file_exactly(input_file, output_file):
    # Open the input file and read all lines
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Write the lines to the output file without modification
    with open(output_file, 'w') as outfile:
        outfile.writelines(lines)

    print(f"Copied contents from '{input_file}' to '{output_file}' exactly as they are.")

# Main function to handle inputs
def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input_file>")
        sys.exit(1)

    # Get the input file from the command line argument
    input_file = sys.argv[1]

    # Define the output file path
    run_dir = os.path.dirname(input_file)
    output_file = os.path.join(run_dir, 'metallic.vset')

    # Copy contents exactly from the input file to metallic.vset
    copy_file_exactly(input_file, output_file)

if __name__ == "__main__":
    main()


