#!/usr/bin/env python3
"""
Replicate Tim's logic and codes, and generalize it for the new workflow,
but include the well injection interval assignment (mat ID = 11),
then the fractures start from 12 onwards.

Modified to LOAD the fracture definitions from 'fractures.json'
instead of prompting for them again, and to read **all** borehole
trajectory files so that near‑borehole assignment works for every well.
"""

import numpy as np
import subprocess
import os
import sys
import json
import glob        # <-- NEW: to collect all borehole files

print("\nGenerating boundary tags text file and material ID text file for the workflow...")

# --------------------------------------------------------------------
# 0) Mesh prefix
# --------------------------------------------------------------------
if len(sys.argv) > 1:
    pre = sys.argv[1]
    print(f"\nUsing provided mesh prefix: {pre}")
else:
    pre = input("\nEnter the prefix for the mesh files ('surf2' in our case): ").strip()

# --------------------------------------------------------------------
# 1) Read node positions and flags
# --------------------------------------------------------------------
pos   = np.genfromtxt(pre + '.1.node', usecols=(1, 2, 3, 4), skip_header=1)
flags = np.genfromtxt(pre + '.1.node', usecols=5,           skip_header=1)

# Initialize boundary tags and material IDs with default values
btags  = np.full(len(pos), -999, dtype=int)
matids = np.ones(len(pos), dtype=int)     # Default to 1 (host rock)

tol = 0.5  # Boundary detection tolerance

# --------------------------------------------------------------------
# 2) .face => boundary markers
# --------------------------------------------------------------------
face_file = pre + '.1.face'
try:
    faces = np.genfromtxt(face_file, dtype=int, skip_header=1, usecols=(1, 2, 3, 4))
    for f in faces:
        nodeA, nodeB, nodeC, bmark = f
        btags[nodeA - 1] = bmark
        btags[nodeB - 1] = bmark
        btags[nodeC - 1] = bmark
    print(f"Read boundary markers from '{face_file}' and assigned them to nodes.")
except IOError:
    print(f"Warning: Could not find '{face_file}'. Skipping face‑based boundary assignment.")

# --------------------------------------------------------------------
# 3) Identify boundary nodes using coordinate conditions with tolerance
# --------------------------------------------------------------------
min_x, max_x = np.min(pos[:, 0]), np.max(pos[:, 0])
min_y, max_y = np.min(pos[:, 1]), np.max(pos[:, 1])
min_z, max_z = np.min(pos[:, 2]), np.max(pos[:, 2])

is_top    = np.isclose(pos[:, 2], max_z, atol=tol)
is_bottom = np.isclose(pos[:, 2], min_z, atol=tol)
is_north  = np.isclose(pos[:, 1], max_y, atol=tol)
is_south  = np.isclose(pos[:, 1], min_y, atol=tol)
is_east   = np.isclose(pos[:, 0], max_x, atol=tol)
is_west   = np.isclose(pos[:, 0], min_x, atol=tol)

two_face_map = {
    (1, 5): 101, (1, 6): 102, (1, 3): 103, (1, 4): 104,
    (2, 5): 105, (2, 6): 106, (2, 3): 107, (2, 4): 108,
    (3, 5): 109, (3, 6): 110, (4, 5): 111, (4, 6): 112,
}

three_face_map = {
    (1, 3, 5): 113, (1, 3, 6): 114, (1, 4, 5): 115, (1, 4, 6): 116,
    (2, 3, 5): 117, (2, 3, 6): 118, (2, 4, 5): 119, (2, 4, 6): 120,
}

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
        btags[i] = two_face_map.get(pair, -999)
    elif len(boundaries) == 3:
        triple = tuple(sorted(boundaries))
        btags[i] = three_face_map.get(triple, -999)
    else:
        btags[i] = -999

# Mark drift nodes (11, -129, -130) as 7
btags[np.where((flags == 11) | (flags == -129) | (flags == -130))[0]] = 7

# Save to _boundaries.txt
np.savetxt(pre + "_boundaries.txt", btags, fmt='%d')
print("\nGenerated boundary tags file.\n")

# --------------------------------------------------------------------
# 4) Read the .ele file and assign material IDs
# --------------------------------------------------------------------
regions = np.genfromtxt(pre + '.1.ele', usecols=(1, 2, 3, 4, 5), dtype=int, skip_header=1)

material_mapping = {
    1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7,
    8: 8, 9: 8, 10: 8, 11: 9, 12: 9, 13: 10,
}

for region in regions:
    node_ids     = np.array(region[:4]) - 1
    material_tag = region[4]
    new_material = material_mapping.get(material_tag, 1)
    unclassified_nodes = node_ids[matids[node_ids] == 1]
    matids[unclassified_nodes] = new_material

# --------------------------------------------------------------------
# 5) Collect **all** borehole trajectories
# --------------------------------------------------------------------
boreholes = {}
for path in glob.glob("needed_files/*_pos_vs_depth.txt"):
    key = os.path.basename(path).split("_")[0].lower()   # e.g. "tc"
    try:
        boreholes[key] = np.genfromtxt(path)
    except Exception as e:
        print(f"  (skip {path}: {e})")

# --------------------------------------------------------------------
# 6) Ask user which borehole is the injection/pressurization well
# --------------------------------------------------------------------
borehole_options = {"TL": "tl", "TS": "ts", "TU": "tu", "TN": "tn", "TC": "tc"}
borehole_choice  = input("Which borehole should the injection/pressurization interval be set on? (TC/TL/TN/TS/TU) [Default: TC]: ").strip().upper()
print()
borehole_key = borehole_options.get(borehole_choice, "tc")  # default "tc"

# --------------------------------------------------------------------
# 7) Set injection interval plane and length
# --------------------------------------------------------------------
def dist_point_to_plane(center, normal, points):
    return np.abs(np.dot(points - center, normal)) / np.linalg.norm(normal)

use_default_plane = input("Use the default fracture plane to locate the borehole interval? (Y/N) [Default: Y]: ").strip().lower()
print()
if use_default_plane == "n":
    cen1  = np.array([float(x) for x in input("Enter fracture plane center (x y z): ").split()])
    print()
    norm1 = np.array([float(x) for x in input("Enter fracture plane normal (x y z): ").split()])
    print()
    norm1 /= np.linalg.norm(norm1)
else:
    cen1  = np.array([1246.0, -881.0, 334.0])
    norm1 = np.array([0.60, 0.40, -0.70])
    norm1 /= np.linalg.norm(norm1)

use_default_length = input("Use the default interval length (2 m)? (Y/N) [Default: Y]: ").strip().lower()
print()
interval_length = 2.0 if use_default_length != "n" else float(input("Enter interval length in meters: "))
print()
interval_half = interval_length / 2.0

# --------------------------------------------------------------------
# 8) Assign mat ID = 11 to injection interval nodes
# --------------------------------------------------------------------
borehole_materials = {"tl": 6, "ts": 7, "tu": 8, "tn": 9, "tc": 10}
material_id        = borehole_materials.get(borehole_key, 10)

borehole_data = boreholes.get(borehole_key, np.array([]))
if borehole_data.size > 0:
    borehole_nodes          = np.where(matids == material_id)[0]
    borehole_dists_to_plane = dist_point_to_plane(cen1, norm1, pos[borehole_nodes, :3])
    borehole_update_inds    = borehole_nodes[np.where(borehole_dists_to_plane < interval_half)[0]]
    matids[borehole_update_inds] = 11

# Save _materials.txt after injection‑well tagging
np.savetxt(pre + "_materials.txt", matids, fmt='%d')
print("Generated materials ID file (after injection‑well assignment).")

########################################################################
# 9) Assign fracture material IDs, starting from 12
########################################################################
with open("fractures.json", "r") as f:
    fractures = json.load(f)

fracture_flags        = np.unique(flags[flags >= 1001])      # 1001, 1002, …
fracture_material_map = {flag: 12 + i for i, flag in enumerate(fracture_flags)}

for flag, mat_id in fracture_material_map.items():
    matids[np.where(flags == flag)[0]] = mat_id

########################################################################
# 10) Near‑Borehole assignment to fracture IDs
########################################################################
dr = 0.2  # plane‑distance half width

for frac in fractures:
    frac_flag    = frac["flag"]
    frac_center  = np.array(frac["center"])
    frac_normal  = np.array(frac["normal"]) / np.linalg.norm(frac["normal"])
    frac_radius  = frac["radius"]

    frac_material = fracture_material_map.get(frac_flag, None)
    if frac_material is None:
        continue

    for name, wb in boreholes.items():
        if wb.size == 0:
            continue

        dp     = dist_point_to_plane(frac_center, frac_normal, wb[:, 1:4])
        dp_min = dp.min()
        if dp_min > 1.0:
            continue

        i1 = np.where(dp == dp_min)[0]
        if len(i1) == 0:
            continue

        pt = wb[i1, 1:4][0]

        dist1 = np.sqrt(
            (pos[:, 0] - pt[0]) ** 2
          + (pos[:, 1] - pt[1]) ** 2
          + (pos[:, 2] - pt[2]) ** 2
        )
        inds = np.where(dist1 < 1.0)[0]

        dist2 = dist_point_to_plane(frac_center, frac_normal, pos[inds, 0:3])
        inds  = inds[np.where(dist2 < (dr / 2.0))[0]]

        rvals = np.sqrt(
            (pos[inds, 0] - frac_center[0]) ** 2
          + (pos[inds, 1] - frac_center[1]) ** 2
          + (pos[inds, 2] - frac_center[2]) ** 2
        )
        inds = inds[np.where(rvals < frac_radius)[0]]

        allowed_ids = [1]  
        inds  = inds[np.where(np.isin(matids[inds], allowed_ids))[0]]
        matids[inds] = frac_material

########################################################################
# 11) Resolve nodes claimed by multiple fractures
########################################################################
if len(fractures) > 1:
    plane_dist_list  = []
    radial_dist_list = []

    for frac in fractures:
        fc = np.array(frac["center"])
        fn = np.array(frac["normal"]) / np.linalg.norm(frac["normal"])

        radial_dist = np.sqrt(
            (pos[:, 0] - fc[0]) ** 2
          + (pos[:, 1] - fc[1]) ** 2
          + (pos[:, 2] - fc[2]) ** 2
        )
        radial_dist_list.append(radial_dist)

        plane_dist = dist_point_to_plane(fc, fn, pos[:, 0:3])
        plane_dist_list.append(plane_dist)

    plane_dist_array  = np.array(plane_dist_list)
    radial_dist_array = np.array(radial_dist_list)
    fracture_radii    = [f["radius"] for f in fractures]
    fracture_mat      = [fracture_material_map[f["flag"]] for f in fractures]

    for nn in range(len(pos)):
        in_fracs = []
        for i, frac in enumerate(fractures):
            if (plane_dist_array[i, nn] < (dr / 2.0)) and (radial_dist_array[i, nn] < fracture_radii[i]):
                in_fracs.append(i)

        if len(in_fracs) >= 2:
            best_idx = min(in_fracs, key=lambda j: plane_dist_array[j, nn])
            if matids[nn] == 1 or matids[nn] in fracture_mat:
                matids[nn] = fracture_mat[best_idx]

########################################################################
# 12) Convert drift nodes (btags == 7) to the closest fracture material
########################################################################
drift_inds = np.where(btags == 7)[0]
if len(drift_inds) > 0 and len(fractures) > 0:

    plane_dist_list  = []
    radial_dist_list = []
    fracture_radii   = []
    for frac in fractures:
        fc = np.array(frac["center"])
        fn = np.array(frac["normal"]) / np.linalg.norm(frac["normal"])

        plane_dist = dist_point_to_plane(fc, fn, pos[:, :3])
        plane_dist_list.append(plane_dist)

        radial_dist = np.sqrt(
            (pos[:, 0] - fc[0]) ** 2 +
            (pos[:, 1] - fc[1]) ** 2 +
            (pos[:, 2] - fc[2]) ** 2
        )
        radial_dist_list.append(radial_dist)

        fracture_radii.append(frac["radius"])

    plane_dist_array  = np.array(plane_dist_list)
    radial_dist_array = np.array(radial_dist_list)
    fracture_mat      = [fracture_material_map[f["flag"]] for f in fractures]

    for node_idx in drift_inds:

        close_fracs = [
            i for i in range(len(fractures))
            if (plane_dist_array[i, node_idx]  < dr) and
               (radial_dist_array[i, node_idx] < fracture_radii[i])
        ]
        if not close_fracs:
            continue

        best_idx = min(close_fracs, key=lambda j: plane_dist_array[j, node_idx])

        if matids[node_idx] == 1 or matids[node_idx] in fracture_mat:
            matids[node_idx] = fracture_mat[best_idx]

# --------------------------------------------------------------------
# 13) Save final _materials.txt
# --------------------------------------------------------------------
np.savetxt(pre + "_materials.txt", matids, fmt='%d')
print("\nGenerated materials ID file (final).")

# --------------------------------------------------------------------
# 14) Debugging statistics
# --------------------------------------------------------------------
print("\nBoundary tag distribution:")
for tag, count in zip(*np.unique(btags, return_counts=True)):
    print(f"  Boundary type {tag:>4}: {count}")

print("\nMaterial ID distribution:")
for mid, count in zip(*np.unique(matids, return_counts=True)):
    print(f"  Material ID {mid:>4}: {count}")

print("\n------------------------------------------------------")

# --------------------------------------------------------------------
# 15) Generate .trn file for px.py visualization
# --------------------------------------------------------------------
print("Generating .trn file for px.py")
with open(pre + ".trn", "w") as fp:
    fp.write("0.000000000000E+00  0.000000000000E+00  0.000000000000E+00\n")

print("\nGenerating material tag .xmf visualization file")
flags_summary_file = pre + "_material_flags_summary.txt"
with open(flags_summary_file, "w") as f:
    f.write(f"{len(matids)} 1\n")
    np.savetxt(f, matids, fmt="%d")
print(f"Copied {pre+'_materials.txt'} to {flags_summary_file} with an added header.")

print("\nRunning px.py on flags summary file...")
try:
    subprocess.run(
        f"python3 px.py -f {pre} {flags_summary_file} meshtags 0.0 tags",
        shell=True,
        check=True,
    )
    print("px.py executed successfully.")
except subprocess.CalledProcessError as e:
    print(f"Error occurred while running px.py: {e}")

print("\n------------------------------------------------------")
