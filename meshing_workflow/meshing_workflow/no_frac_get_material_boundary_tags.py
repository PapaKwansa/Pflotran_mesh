#!/usr/bin/env python3

"""
Replicate Tim's logic and codes, and generalize it for the new workflow,
but include the well injection interval assignment (mat ID=11),
then the fractures start from 12 onwards.
"""

import numpy as np
import subprocess
import os
import sys

print(f"\nGenerating boundary tags text file and material ID text file for the workflow...")

# Ask the user for the prefix name
if len(sys.argv) > 1:
    pre = sys.argv[1]
    print(f"\nUsing provided mesh prefix: {pre}")
else:
    pre = input("\nEnter the prefix for the mesh files ('surf2' in our case): ").strip() 

# Read node positions and flags
pos = np.genfromtxt(pre + '.1.node', usecols=(1,2,3,4), skip_header=1)
flags = np.genfromtxt(pre + '.1.node', usecols=5, skip_header=1)

# Initialize boundary tags and material IDs with default values
btags = np.full(len(pos), -999, dtype=int)
matids = np.ones(len(pos), dtype=int)  # Default to 1 (host rock)

tol = 0.5  # Adjust tolerance for boundary detection

# --------------------------------------------------------------------
# .face => boundary markers
# --------------------------------------------------------------------
face_file = pre + '.1.face'
try:
    faces = np.genfromtxt(face_file, dtype=int, skip_header=1, usecols=(1,2,3,4))
    for f in faces:
        nodeA, nodeB, nodeC, bmark = f
        btags[nodeA - 1] = bmark
        btags[nodeB - 1] = bmark
        btags[nodeC - 1] = bmark
    print(f"Read boundary markers from '{face_file}' and assigned them to nodes.")
except IOError:
    print(f"Warning: Could not find '{face_file}'. Skipping face-based boundary assignment.")

# --------------------------------------------------------------------
# Identify boundary nodes using coordinate conditions with tolerance
# --------------------------------------------------------------------

# Compute min/max values for boundary detection
min_x, max_x = np.min(pos[:, 0]), np.max(pos[:, 0])
min_y, max_y = np.min(pos[:, 1]), np.max(pos[:, 1])
min_z, max_z = np.min(pos[:, 2]), np.max(pos[:, 2])

# Logical masks for boundary detection
is_top    = np.isclose(pos[:, 2], max_z, atol=tol)
is_bottom = np.isclose(pos[:, 2], min_z, atol=tol)
is_north  = np.isclose(pos[:, 1], max_y, atol=tol)
is_south  = np.isclose(pos[:, 1], min_y, atol=tol)
is_east   = np.isclose(pos[:, 0], max_x, atol=tol)
is_west   = np.isclose(pos[:, 0], min_x, atol=tol)

# Define single-face, two-face, and three-face combinations
two_face_map = { 
    (1, 5): 101, (1, 6): 102, (1, 3): 103, (1, 4): 104, 
    (2, 5): 105, (2, 6): 106, (2, 3): 107, (2, 4): 108, 
    (3, 5): 109, (3, 6): 110, (4, 5): 111, (4, 6): 112 
}

three_face_map = { 
    (1, 3, 5): 113, (1, 3, 6): 114, (1, 4, 5): 115, (1, 4, 6): 116, 
    (2, 3, 5): 117, (2, 3, 6): 118, (2, 4, 5): 119, (2, 4, 6): 120 
}

# Assign boundary tags
for i in range(len(pos)):
    boundaries = []
    if is_top[i]:    boundaries.append(1)
    if is_bottom[i]: boundaries.append(2)
    if is_north[i]:  boundaries.append(3)
    if is_south[i]:  boundaries.append(4)
    if is_east[i]:   boundaries.append(5)
    if is_west[i]:   boundaries.append(6)

    if len(boundaries) == 1:
        btags[i] = boundaries[0]
    elif len(boundaries) == 2:
        pair = tuple(sorted(boundaries))
        btags[i] = two_face_map.get(pair, 999)  
    elif len(boundaries) == 3:
        triple = tuple(sorted(boundaries))
        btags[i] = three_face_map.get(triple, 999)  
    else:
        btags[i] = -999

# Mark drift nodes (11, -129, -130) as 7
btags[np.where((flags == 11) | (flags == -129) | (flags == -130))[0]] = 7

# Save to _boundaries.txt
np.savetxt(pre + "_boundaries.txt", btags, fmt='%d')
print("\nGenerated boundary tags file.")
print('\n')

# --------------------------------------------------------------------
# Read the .ele file and assign material IDs
# --------------------------------------------------------------------
regions = np.genfromtxt(pre + '.1.ele', usecols=(1,2,3,4,5), dtype=int, skip_header=1)

# Predefined material mapping
material_mapping = {
    1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7,
    8: 8, 9: 8, 10: 8, 11: 9, 12: 9, 13: 10
}

for region in regions:
    node_ids = np.array(region[:4]) - 1  
    material_tag = region[4]
    new_material = material_mapping.get(material_tag, 1)  
    unclassified_nodes = node_ids[matids[node_ids] == 1]
    matids[unclassified_nodes] = new_material

# --------------------------------------------------------------------
# Identify wellbore injection interval and assign mat ID=11
# --------------------------------------------------------------------
# Function to calculate distance from points to a plane
def dist_point_to_plane(center, normal, points):
    return np.abs(np.dot(points - center, normal)) / np.linalg.norm(normal)

# Ask user for borehole selection
borehole_options = {"TL": "tl", "TS": "ts", "TU": "tu","TN": "tn", "TC": "tc"}
borehole_choice = input("Which borehole should the injection/pressurization interval be set on? (TC/TL/TN/TS/TU) [Default: TC]: ").strip().upper()
print('\n')
borehole_key = borehole_options.get(borehole_choice, "tc")

# Load borehole data
borehole_file = f'needed_files/{borehole_key}_pos_vs_depth.txt'
boreholes = {borehole_key: np.genfromtxt(borehole_file)}

# Ask user for fracture plane parameters
use_default_plane = input("To locate the borehole interval, imagine a fracture plane intersecting it. Which fracture plane do you prefer? Use the default plane? (Y/N) [Default: Y]: ").strip().lower()
print('\n')
if use_default_plane == "n":
    cen1 = np.array([float(x) for x in input("Enter fracture plane center coordinates (x, y, z) separated by spaces: ").split()])
    print('\n')
    norm1 = np.array([float(x) for x in input("Enter fracture plane normal vector (x, y, z) separated by spaces: ").split()])
    print('\n')
    norm1 /= np.linalg.norm(norm1)
else:
    cen1 = np.array([1245.0, -884.0, 334.58])
    norm1 = np.array([1.0, 0.37, 0.18])
    norm1 /= np.linalg.norm(norm1)

# Ask user for interval length
use_default_length = input("Would you like to use the default interval length (2 m)? (Y/N) [Default: Y]: ").strip().lower()
print('\n')
if use_default_length == "n":
    interval_length = float(input("Enter the desired interval length in meters: "))
    print('\n')
else:
    interval_length = 20.0
interval_half = interval_length / 2.0

# Define material ID mapping
borehole_materials = {"tl": 6,"ts": 7, "tu": 8, "tn": 9, "tc": 10}
material_id = borehole_materials.get(borehole_key, 10)

# Process borehole points
borehole_data = boreholes.get(borehole_key, np.array([]))
if borehole_data.size > 0:
    borehole_nodes = np.where(matids == material_id)[0]
    borehole_dists_to_plane = dist_point_to_plane(cen1, norm1, pos[borehole_nodes, :3])
    borehole_update_inds = borehole_nodes[np.where(borehole_dists_to_plane < interval_half)[0]]
    matids[borehole_update_inds] = 11

# Save Updated _materials.txt File
np.savetxt(pre + "_materials.txt", matids, fmt='%d')

# Assign permeable lens in sandstone
lens_radius = 1000.0  # 1km radius around wellbore
for i in range(len(pos)):
    if matids[i] == 2:  # bartlesville_sand
        dist_xy = np.sqrt((pos[i,0] - 5000.0)**2 + (pos[i,1] - 5000.0)**2)
        if dist_xy < lens_radius:
            matids[i] = 6
print("\nGenerated materials ID file.")

# --------------------------------------------------------------------
# Debugging output: count occurrences of each boundary tag and material ID
# --------------------------------------------------------------------
print("\nBoundary tag distribution:")
unique_btags, counts = np.unique(btags, return_counts=True)
for tag, count in zip(unique_btags, counts):
    print(f"Boundary type {tag}: {count} cells")

print("\nMaterial ID distribution:")
unique_matids, counts = np.unique(matids, return_counts=True)
for mid, count in zip(unique_matids, counts):
    print(f"Material ID {mid}: {count} cells")

print("\n------------------------------------------------------")
