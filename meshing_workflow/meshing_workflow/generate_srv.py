import os
import shutil
import sys

# Directory containing the needed files
needed_files_dir = "needed_files"

# Check for mesh_name argument
if len(sys.argv) < 2:
    script_name = os.path.basename(sys.argv[0])
    print(f"Usage: python {script_name} mesh_name")
    sys.exit(1)

# Read mesh_name from the command line
mesh_name = sys.argv[1].strip()

# Define source and destination file paths
src_srv_file = os.path.join(needed_files_dir, "example.srv")
dst_srv_file = f"{mesh_name}.srv"

# Copy and rename the file
shutil.copy(src_srv_file, dst_srv_file)

print(f"Copied '{src_srv_file}' to current directory as '{dst_srv_file}'.")
