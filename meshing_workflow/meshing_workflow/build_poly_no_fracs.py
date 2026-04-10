#!/usr/bin/env python3

import numpy as np
import numpy.linalg as lg
import subprocess
import sys
import os
import shutil

###############################################################################
# 1) Function to copy .poly file from needed_files/ to the current directory
###############################################################################

def copy_poly_file(filename):
    """
    Copy the .poly file from the 'needed_files' subfolder to the current directory.
    """
    source_path = os.path.join("needed_files", filename)  # Source file in needed_files/
    destination_path = filename  # Destination file in current directory

    if not os.path.exists(source_path):
        print(f"Error: File '{source_path}' not found.")
        sys.exit(1)

    try:
        shutil.copy(source_path, destination_path)
        print(f"Copied '{source_path}' to '{destination_path}'.")
    except Exception as e:
        print(f"Error copying file: {e}")
        sys.exit(1)

###############################################################################
# 2) The main routine
###############################################################################

def main():
    """
    1) Use provided filename argument or prompt user for input.
    2) Copy the .poly file from needed_files/ to current directory.
    3) Run TetGen.
    """
    # Check if filename is passed as a command-line argument
    if len(sys.argv) > 1:
        filename = sys.argv[1].strip()
        tetgen_exe = sys.argv[2].strip()
        print(f"\nUsing provided .poly filename: {filename}")
    else:
        # Prompt user for input if no argument is provided
        filename = input("\nEnter .poly filename ('surf2.poly' in our case): ").strip()
    
    # Ensure the filename has the correct extension
    if not filename.endswith('.poly'):
        filename += '.poly'

    print(f"--> Processing .poly file: {filename}")

    # Step 1: Copy the .poly file from needed_files/ to current directory
    copy_poly_file(filename)

    # Step 2: Run TetGen
    print("\n--> Running TetGen...")
    cmd = f"{tetgen_exe} -pnq1.3a1e12aAA {filename}"
    subprocess.run(cmd, shell=True)
    print("--> Done TetGen!")

if __name__ == "__main__":
    main()
