import numpy as np
import os
import sys

# Directory where needed files are stored
needed_files_dir = "needed_files"

# Build full path to "original_metallic.ex" within needed_files/
metallic_ex_file = os.path.join(needed_files_dir, "original_metallic.ex")

# Ensure a mesh_name argument is provided
if len(sys.argv) < 2:
    script_name = os.path.basename(sys.argv[0])
    sys.exit(1)

# Read mesh_name from the command line
mesh_name = sys.argv[1].strip()

# Construct the .1.node file path in the current directory
node_file = f"{mesh_name}.1.node"

# Example debug prints (optional)
print(f"Reading metallic ex file from: {metallic_ex_file}")
print(f"Reading node file from current directory: {node_file}")
updated_metallic_ex_file = "metallic.ex"
updated_vset_file = "metallic.vset"

def find_nearest_nodes(metallic_coords, metallic_data, node_coords, node_indices):
    """
    Find the nearest node in .node file for each metallic node and update only the index.
    """
    matched_data = []  # Store updated metallic data
    matched_indices = []  # Store new indices

    for i, metallic_coord in enumerate(metallic_coords):
        distances = np.linalg.norm(node_coords - metallic_coord, axis=1)  # Compute Euclidean distances
        min_index = np.argmin(distances)  # Find the index of the nearest node
        closest_node_index = node_indices[min_index]  # Get the corresponding node index

        # Update only the index in metallic.ex, keeping coordinates unchanged
        metallic_data[i][0] = closest_node_index  
        matched_data.append(metallic_data[i])
        matched_indices.append(closest_node_index)

    return matched_indices, matched_data

def sort_by_new_index(matched_indices, updated_metallic_data):
    """
    Sort metallic data based on the updated indices.
    """
    sorted_data = sorted(zip(matched_indices, updated_metallic_data), key=lambda x: int(x[0]))
    sorted_indices = [item[0] for item in sorted_data]
    sorted_metallic_data = [" ".join(item[1]) + "\n" for item in sorted_data]
    return sorted_indices, sorted_metallic_data

# Read metallic.ex file and extract metallic node coordinates
with open(metallic_ex_file, "r") as f:
    metallic_lines = f.readlines()

# Ignore the first line, extract coordinates from columns 2-4
metallic_coords = []
metallic_data = []
for line in metallic_lines[1:]:
    parts = line.split()
    if len(parts) >= 5:
        coord = np.array(list(map(float, parts[1:4])))  # Extract XYZ coordinates
        metallic_coords.append(coord)
        metallic_data.append(parts)  # Store full line data for updating

# Read TC_hollow.1.node file and extract all node coordinates and indices
node_coords = []
node_indices = []

with open(node_file, "r") as f:
    node_lines = f.readlines()

for line in node_lines:
    parts = line.split()
    if len(parts) >= 4:
        node_index = parts[0]
        coord = np.array(list(map(float, parts[1:4])))
        node_coords.append(coord)
        node_indices.append(node_index)

# Convert to NumPy arrays for efficient calculations
node_coords = np.array(node_coords)
metallic_coords = np.array(metallic_coords)

# Step 1: Find the nearest nodes and update only the indices
matched_indices, updated_metallic_data = find_nearest_nodes(metallic_coords, metallic_data, node_coords, node_indices)

# Step 2: Sort by new indices
sorted_indices, sorted_metallic_data = sort_by_new_index(matched_indices, updated_metallic_data)

# Step 3: Save metallic_updated.ex with updated indices but original coordinates
with open(updated_metallic_ex_file, "w") as f:
    f.write("CONNECTIONS\t3250\n")  # Rewrite header
    f.writelines(sorted_metallic_data)

# Step 4: Save metallic_updated.vset (only sorted indices)
with open(updated_vset_file, "w") as f:
    f.writelines("\n".join(sorted_indices) + "\n")

print(f"Updated metallic.ex saved as {updated_metallic_ex_file}")
print(f"Sorted matched node indices saved as {updated_vset_file}")
