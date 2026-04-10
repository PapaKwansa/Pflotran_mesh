import subprocess
import os
import sys
import glob
import multiprocessing as mp

# Define paths
voronoi_exe = "/home/yao/Research/voronoi/src/voronoi" # PLEASE UPDATE IT WITH YOUR LOCAL PATH
tetgen_exe = "/home/yao/Desktop/E4D/E4D-master/third_party/tetgen1.5.0/tetgen" # PLEASE UPDATE IT WITH YOUR LOCAL PATH

vset_files = [
    "top.vset", "bottom.vset", "north.vset", "south.vset", 
    "east.vset", "west.vset", "drift.vset", "well_interval.vset",
    "TC.vset", "TL.vset", "TS.vset", "TU.vset", "TN.vset", 
    "well_non_interval.vset", "TL_cell.vset", "TS_cell.vset", 
    "TU_cell.vset", "TN_cell.vset", "TC_cell.vset"
]

print('\n' + '*' * 80 + '\n')
print(f"--> Thank you for using the meshing tool from CUSSP!")
print('\n' + '*' * 80 + '\n')
print(f"--> Starting meshing workflow for the 4100-level testbed at SURF...")
print('\n' + '*' * 80 + '\n')

# Delete specific files in the current directory (not in subfolders)
print(f"--> Deleting previously generated files...")
file_extensions = [
    "*.vset", "*.ex", "*.edge", "*.ele", "*.face", "*.neigh", "*.node", "*.h5", "*.srv",
    "*.inp", "*.mapping", "*.poly", "*.uge", "*.ugi", "*.txt", "*.pvtp", "*.vtp", "*.vtu"
]

for ext in file_extensions:
    for file in glob.glob(ext):  # Only matches files in the current folder
        try:
            os.remove(file)
            print(f"Deleted: {file}")
        except Exception as e:
            print(f"Error deleting {file}: {e}")
print('\n')
print(f"--> Files deleted...")
print('\n' + '*' * 80 + '\n')

# Ask user to confirm voronoi.exe path
path_confirm = input("Have you installed TetGen and VORONOI? Have you specified your local paths for voronoi.exe and tetgen.exe in workflow.py? (Y/N): ").strip().lower()
if path_confirm != 'y':
    print('\n' + '*' * 80 + '\n')
    print("Please install TetGen and VORONOI, and ensure the paths for voronoi.exe in workflow.py and tetgen.exe in build_poly_~.py are correctly set before proceeding.")
    print('\n' + '*' * 80 + '\n')
    sys.exit()
print('\n' + '*' * 80 + '\n')
# Ask user if there are fractures in the system
has_fractures = input("Do you want to have fractures in the system? (Y/N): ").strip().lower()
print('\n' + '*' * 80 + '\n')

if has_fractures in ['Y', 'y']:
    mesh_name = "surf2"
    poly_script = "build_poly_fracs.py"
    tag_script = "frac_get_material_boundary_tags.py"
else:
    mesh_name = "surf2" #in case the poly file is different
    poly_script = "build_poly_no_fracs.py"
    tag_script = "no_frac_get_material_boundary_tags.py"

# File definitions based on mesh_name
node_file = f"{mesh_name}.1.node"
ele_file = f"{mesh_name}.1.ele"
avs_file = f"{mesh_name}.inp"
h5_file = f"{mesh_name}_material_ids.h5"
uge_file = f"{mesh_name}.uge"
output_ugi = f"{mesh_name}.ugi"
boundary_file = f"{mesh_name}_boundaries_updated.txt"
material_file = f"{mesh_name}_materials.txt"

# Run the appropriate poly script to generate TetGen files
print(f"--> Running {poly_script} to process the .poly file...")
print('\n')
command = f"python3 {poly_script} {mesh_name} {tetgen_exe}"
subprocess.run(command, shell=True)
print('\n' + '*' * 80 + '\n')

# Generate boundary tag and material ids files
print(f"--> Running {tag_script} to generate boundary tag and material ids files...")
command = f"python3 {tag_script} {mesh_name}"
subprocess.run(command, shell=True)
print('\n' + '*' * 80 + '\n')

# Delete a row in .node
print(f"--> Deleting a non-data row from the .node file to prevent potential bugs in subsequent steps...")
command = f"python3 delete_node_row.py {mesh_name}"
subprocess.run(command, shell=True)
print('\n' + '*' * 80 + '\n')

# Generate AVS mesh from TetGen files
print(f"--> Generating AVS mesh from TetGen files...")
command = f"python3 tetgen_to_avs.py {node_file} {ele_file} {avs_file}"
subprocess.run(command, shell=True)
print('\n' + '*' * 80 + '\n')

# Generate PFLOTRAN material h5 file
print(f"--> Generating PFLOTRAN material h5 file...")
command = f"python3 material_h5_from_txt.py {node_file} {h5_file} {material_file}"
subprocess.run(command, shell=True)
print('*' * 80 + '\n')

# Correct material files
print(f"--> Correcting borehole material tags...")
command = f"python3 correct_materials.py {mesh_name}"
subprocess.run(command, shell=True)
print('*' * 80 + '\n')

# Correct boundary files
print(f"--> Correcting borehole boundary tags...")
command = f"python3 correct_boundary_tags.py {mesh_name}"
subprocess.run(command, shell=True)
print('*' * 80 + '\n')

# Generate PFLOTRAN uge mesh file
print(f"--> Generating PFLOTRAN uge file...")
command = f"{voronoi_exe} -avs {avs_file} -type pflotran -o {uge_file}"
subprocess.run(command, shell=True)
print('*' * 80 + '\n')

# Extract borehole nodes as vertex sets
print(f"--> Extracting vertex sets for borehole nodes...")
command = f"python3 boundary_nodes_to_vset.py {boundary_file}"
subprocess.run(command, shell=True)
print('*' * 80 + '\n')

# Convert vertex sets to PFLOTRAN ex files
print(f"--> Converting vertex sets to PFLOTRAN ex files...")
for vset_file in vset_files:
    command = f"python3 convert_vset_to_ex.py {uge_file} {vset_file}"
    subprocess.run(command, shell=True)
print('*' * 80 + '\n')

# Generate metallic files
print(f"--> Generating metallic .vset and .ex files...")
command = f"python3 metallic.py {mesh_name}"
subprocess.run(command, shell=True)
print('*' * 80 + '\n')

# Generate .ugi files
print(f"--> Generating .ugi files...")
command = f"python3 generate_ugi.py {ele_file} {node_file} {output_ugi}"
subprocess.run(command, shell=True)
print('*' * 80 + '\n')

# Delete a row in .ugi
print(f"--> Deleting a non-data row from the .ugi file to prevent potential bugs in subsequent steps...")
command = f"python3 delete_ugi_row.py {mesh_name}"
subprocess.run(command, shell=True)
print('*' * 80 + '\n')

import subprocess

# Generate .mapping and _all.vset files
print(f"--> Generating .mapping and _all.vset files...")
command = f"python3 mapping.py {mesh_name}"
subprocess.run(command, shell=True)
print('*' * 80 + '\n')

# Generate .srv file
print(f"--> Generating .srv file...")
command = f"python3 generate_srv.py {mesh_name}"
subprocess.run(command, shell=True)
print('*' * 80 + '\n')

print(f"--> Meshing workflow is done!")
print('*' * 80 + '\n')

# Voronoi visualization
user_input = input("Would you like to generate visualization files for the Voronoi mesh? (This process may take several hours...) (Y/N): ").strip().lower()
if user_input in ['yes', 'y']:
    print(f"--> Generating visualization files...")
    command = f"python3 new_tet_voro_vtu_manual.py"
    result = subprocess.run(command, shell=True)

    # Check for errors
    if result.returncode != 0:
        print(f"Error executing {command}:")
        print(result.stderr)
        sys.exit(1)

    print('*' * 80 + '\n')
else:
    print("Visualization generation skipped. Exiting program.")
    sys.exit()