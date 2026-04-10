#!/usr/bin/env python3

""" 
Learned from and replicate Tim's build_mesh_wtih_2fracs.py code.
"""

import numpy as np
import numpy.linalg as lg
import subprocess
import shutil
import sys
import os
import json  # <-- ADDED for saving fracture definitions

###############################################################################
# 1) Get mesh name from user or command line
###############################################################################

if len(sys.argv) > 1:
    mesh_name = sys.argv[1].strip()
    tetgen_exe = sys.argv[2].strip()
    print(f"\nUsing provided mesh name: {mesh_name}")
else:
    mesh_name = input("\nEnter mesh name (without .poly extension): ").strip()

poly_filename = f"{mesh_name}.poly"

shutil.copy(f"needed_files/{poly_filename}", poly_filename)
print("\n------------------------------------------------------")

print(f"\nBegin the process of building a new {poly_filename} file...")

###############################################################################
# 2) Function to copy .poly file from needed_files/ to the current directory
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
# 1) Helper functions for plane distance & checking point inside triangle
###############################################################################

def dist_point_to_plane(plane_pt, plane_normal, frac_points):
    """
    Learned from Tim's code
    Computes signed distance from each fracture point in 'frac_points' (N,4)
    to the plane defined by 'plane_pt' (3,) and 'plane_normal' (3,).
    Returns a shape (N,) array of distances.
    """
    plane_pt = np.array(plane_pt, dtype=np.float64)
    plane_normal = np.array(plane_normal, dtype=np.float64)

    pts_xyz = frac_points[:, 0:3].astype(np.float64)  # (N,3)
    vecs = pts_xyz - plane_pt
    return np.einsum('ij,j->i', vecs, plane_normal)

def point_in_triangle(pt, tri_pts):
    """
    Checks if 3D point 'pt' is inside or on edges of triangle 'tri_pts' (3x3).
    Assumes 'pt' is already projected onto the plane of 'tri_pts'.
    Returns True/False.
    """
    v0 = tri_pts[2] - tri_pts[0]
    v1 = tri_pts[1] - tri_pts[0]
    v2 = pt - tri_pts[0]

    dot00 = np.dot(v0, v0)
    dot01 = np.dot(v0, v1)
    dot02 = np.dot(v0, v2)
    dot11 = np.dot(v1, v1)
    dot12 = np.dot(v1, v2)

    denom = (dot00 * dot11 - dot01 * dot01)
    if abs(denom) < 1e-14:
        return False  # degenerate or near-degenerate
    inv_denom = 1.0 / denom
    u = (dot11 * dot02 - dot01 * dot12) * inv_denom
    v = (dot00 * dot12 - dot01 * dot02) * inv_denom
    return (u >= 0) and (v >= 0) and (u + v <= 1.0)

###############################################################################
# 2) Build fracture points
###############################################################################

def build_fracture_points(center, radius, normal, spacing, flag,
                          max_horizontal_distance=1e6,
                          upper_bound=1e6, lower_bound=-1e6):
    """
    Replication of Tim's 'build_frac_pts' logic:
      1) r = radius + 5, rkeep = radius
      2) Build rotation from [0,0,1] to 'normal'.
      3) Generate grid [-r, r] x [-r, r], step=spacing.
      4) Rotate + translate
      5) 3D spherical cull => dist <= rkeep
      6) Horizontal distance & Z-bounds cull
    """
    normal = np.array(normal, dtype=float)
    normal /= lg.norm(normal)
    m_norm = np.array([0.0, 0.0, 1.0], dtype=float)

    r = radius + 5.0
    rkeep = radius

    costheta = np.dot(m_norm, normal) / (lg.norm(m_norm)*lg.norm(normal))
    axis = np.cross(m_norm, normal)
    axis /= np.linalg.norm(axis)

    c = costheta
    s = np.sqrt(1.0 - c*c)
    C = 1.0 - c
    x, y, z = axis
    rmat = np.zeros((3,3), dtype=float)

    rmat[0,0] = x*x*C + c
    rmat[0,1] = x*y*C - z*s
    rmat[0,2] = x*z*C + y*s
    rmat[1,0] = y*x*C + z*s
    rmat[1,1] = y*y*C + c
    rmat[1,2] = y*z*C - x*s
    rmat[2,0] = z*x*C - y*s
    rmat[2,1] = z*y*C + x*s
    rmat[2,2] = z*z*C + c

    pts_list = []
    for yy in np.arange(-r, r, spacing):
        for xx in np.arange(-r, r, spacing):
            pts_list.append([xx, yy, 0.0, flag])
    pts = np.array(pts_list, dtype=float)

    pts_rot = pts.copy()
    for i in range(len(pts_rot)):
        old_xyz = pts[i,0:3]
        new_xyz = np.dot(rmat, old_xyz)
        pts_rot[i,0:3] = new_xyz + center

    dist_3d = np.sqrt(
        (pts_rot[:,0] - center[0])**2 +
        (pts_rot[:,1] - center[1])**2 +
        (pts_rot[:,2] - center[2])**2
    )
    keep_idx = np.where(dist_3d <= rkeep)[0]
    pts_rot = pts_rot[keep_idx]

    dist_horiz = np.sqrt(
        (pts_rot[:,0] - center[0])**2 +
        (pts_rot[:,1] - center[1])**2
    )
    keep_idx = np.where(
        (dist_horiz <= max_horizontal_distance) &
        (pts_rot[:,2] > lower_bound) &
        (pts_rot[:,2] < upper_bound)
    )[0]
    pts_final = pts_rot[keep_idx]

    print(f"Generated {len(pts_final)} valid fracture points for flag {flag}.")
    return pts_final

###############################################################################
# 3) add_fracture_points_to_poly that REMOVES points near facets
###############################################################################

def add_fracture_points_to_poly(poly_file, fractures):
    """
    Replicate Tim's codes
    1) Reads entire .poly file
    2) Parse node lines => build node_id->(x,y,z) for plane checks
    3) Store facet-and-beyond lines verbatim (including any user PLCs, zones, etc.)
    4) Parse triangular facets in memory => remove near-facet points (0.5*dr)
    5) Deduplicate & append final fracture points
    6) Write everything back => the facet section is appended verbatim
    """

    with open(poly_file, 'r') as f:
        lines = f.readlines()

    # 3.1) parse total_nodes from first line
    header = lines[0].strip().replace(',', '').split()
    total_nodes = int(header[0])

    # find "Total number of facets"
    facet_line_index = None
    for idx, line in enumerate(lines):
        if "Total number of facets" in line:
            facet_line_index = idx
            break
    if facet_line_index is None:
        raise ValueError("Could not find 'Total number of facets' in the .poly file.")

    # node lines => [1 : facet_line_index]
    node_lines = lines[1:facet_line_index]
    # everything from "Total number of facets" onward => keep verbatim
    facet_and_beyond = lines[facet_line_index:]

    # build node_id->(x,y,z)
    node_coords = {}
    for line in node_lines:
        toks = line.strip().split()
        if len(toks) < 4:
            continue
        try:
            node_id = int(toks[0])
            xx = np.float64(toks[1])
            yy = np.float64(toks[2])
            zz = np.float64(toks[3])
            node_coords[node_id] = (xx, yy, zz)
        except:
            pass

    # parse how many facets
    # e.g. "10339           1   #Total number of facets..."
    facet_header = facet_and_beyond[0].strip().split()
    total_facets = int(facet_header[0])

    facet_def_lines = facet_and_beyond[1 : 1 + total_facets]
    rest_lines = facet_and_beyond[1 + total_facets : ]

    # parse triangular facets from facet_def_lines
    tri_facets = []
    idx_line = 0
    while idx_line < len(facet_def_lines):
        raw_line = facet_def_lines[idx_line].strip()
        idx_line += 1
        # skip comment lines or empty lines
        if not raw_line or raw_line.startswith('#'):
            continue
        parts = raw_line.split()
        if len(parts) < 2:
            continue
        # subpoly_count, hole_count
        try:
            subpoly_count = int(parts[0])
            hole_count = int(parts[1])
        except ValueError:
            continue

        # next subpoly_count lines => corners
        for sp in range(subpoly_count):
            if idx_line >= len(facet_def_lines):
                break
            poly_line = facet_def_lines[idx_line].strip().split()
            idx_line += 1
            if not poly_line or len(poly_line) < 4:
                continue
            corner_count = int(poly_line[0])
            if corner_count == 3:
                nA = int(poly_line[1])
                nB = int(poly_line[2])
                nC = int(poly_line[3])
                tri_facets.append((nA, nB, nC))

        # skip hole_count lines if they exist
        for _ in range(hole_count):
            idx_line += 1

    # tri_facets => list of (nodeA, nodeB, nodeC)
    # remove near-facet points
    for (nA, nB, nC) in tri_facets:
        if (nA not in node_coords) or (nB not in node_coords) or (nC not in node_coords):
            continue
        p0 = np.array(node_coords[nA], dtype=np.float64)
        p1 = np.array(node_coords[nB], dtype=np.float64)
        p2 = np.array(node_coords[nC], dtype=np.float64)
        v1 = p1 - p0
        v2 = p2 - p0
        normal = np.cross(v1, v2)
        ln = np.linalg.norm(normal)
        if ln < 1e-12:
            continue
        normal /= ln

        tri_pts = np.vstack([p0, p1, p2])

        # remove from each fracture
        for frac in fractures:
            if "points" not in frac or "spacing" not in frac:
                continue
            if len(frac["points"]) == 0:
                continue
            dr = frac["spacing"]
            pts_arr = frac["points"]
            dvals = dist_point_to_plane(p0, normal, pts_arr)
            close_idx = np.where(np.abs(dvals) <= 0.5 * dr)[0]
            if len(close_idx) == 0:
                continue

            # project & check inside
            pproj = pts_arr[close_idx, 0:3].copy()
            for i_i, idxv in enumerate(close_idx):
                pproj[i_i] -= dvals[idxv] * normal

            inside_mask = np.zeros(len(close_idx), dtype=bool)
            for i_i, pt3d in enumerate(pproj):
                inside_mask[i_i] = point_in_triangle(pt3d, tri_pts)

            if not np.any(inside_mask):
                continue

            used_mask = np.ones(len(pts_arr), dtype=bool)
            used_mask[close_idx[inside_mask]] = False
            frac["points"] = frac["points"][used_mask]

    # done removing near-facet points

    # figure out last node index
    last_point_index = total_nodes
    if len(node_coords):
        max_existing_id = max(node_coords.keys())
        last_point_index = max(last_point_index, max_existing_id)

    # now deduplicate final fracture points, build new node lines
    new_points = []
    point_map = {}
    for frac in fractures:
        pts_arr = frac["points"]
        for row in pts_arr:
            xval, yval, zval, flg = row
            # Round the values to keep things neat
            xval = float(f"{xval:.6E}")
            yval = float(f"{yval:.5E}")
            zval = float(f"{zval:.5E}")
            key = (xval, yval, zval, flg)
            last_point_index += 1
            point_map[key] = last_point_index
            line_str = (f"      {last_point_index:<5d}"
                        f"   {xval:.6E}"
                        f"   {yval:.5E}"
                        f"   {zval:.5E}"
                        f"         1         {int(flg)}")
            new_points.append(line_str)

    total_nodes = last_point_index

    # now rewrite the entire .poly:
    #   line 0 => updated node count (with commas)
    #   node_lines => unchanged
    #   new_points => appended
    #   facet_and_beyond => unchanged
    with open(poly_file, 'w') as f:
        # update the first line with new node count (with commas)
        f.write(f"{total_nodes}, 3, 1, 1,\n")
        # old node lines
        f.writelines(node_lines)
        # new node lines
        if new_points:
            f.write("\n".join(new_points) + "\n")
        # everything else => the same lines
        f.writelines(facet_and_beyond)

    print(f"\nUpdated .poly file '{poly_file}' with {len(new_points)} new fracture points, "
          f"after removing near-facet points.")

###############################################################################
# 4) The main routine
###############################################################################

def main():
    """
    1) Either use default fractures or prompt user for parameters.
    2) build_fracture_points (Tim's code)
    3) add_fracture_points_to_poly => parse facets, remove near-facet points,
       then write them, leaving facet lines verbatim
    4) Optionally run TetGen
    """
    default_fractures = [
        {
            "center": [1245.0, -884.0, 334.58],
            "radius": 25.0,
            "normal": [1.0, 0.37, 0.18],
            "spacing": 0.2,
            "flag": 1001
        },
        {
            "center": [1226.5, -878.5, 344.735],
            "radius": 25.0,
            "normal": [-0.5, 1.0, 0.55],
            "spacing": 0.2,
            "flag": 1002
        },
        {
            "center": [1252.0, -883.0, 337.0],
            "radius": 15.0,
            "normal": [-0.2, 0.85, -0.47],
            "spacing": 0.2,
            "flag": 1003
        },
        {
            "center": [1230.0, -890.0, 345.0],
            "radius": 30.0,
            "normal": [0.7, -0.3, 0.65],
            "spacing": 0.2,
            "flag": 1004
        },
        {
            "center": [1248.0, -875.0, 342.0],
            "radius": 20.0,
            "normal": [0.2, 0.1, 1.0],
            "spacing": 0.2,
            "flag": 1005
        },
    ]

    use_defaults = input("\nDo you want to use the default fractures (5 fracs)? (Y/N): ").strip().upper()
    fractures = []

    if use_defaults == "Y":
        for i, frac in enumerate(default_fractures, start=1):
            print(f"\nUsing predefined parameters for fracture {i}:")
            print(f"  Center: {frac['center']}")
            print(f"  Radius: {frac['radius']}")
            print(f"  Normal: {frac['normal']}")
            print(f"  Spacing: {frac['spacing']}")
            print(f"  Flag: {frac['flag']}")
            pts = build_fracture_points(
                center=frac["center"],
                radius=frac["radius"],
                normal=frac["normal"],
                spacing=frac["spacing"],
                flag=frac["flag"]
            )
            frac["points"] = pts
            fractures.append(frac)
    else:
        num_fract = int(input("Enter number of fractures: ").strip())
        for i in range(num_fract):
            print(f"\nDefining fracture {i+1}:")
            cx, cy, cz = map(np.float64, input("Enter center (x y z): ").split())
            r = np.float64(input("Enter radius: ").strip())
            nx, ny, nz = map(np.float64, input("Enter normal (x y z): ").split())
            spacing = np.float64(input("Enter spacing: ").strip())
            flag = int(input("Enter boundary marker flag: ").strip())

            pts = build_fracture_points(
                center=[cx, cy, cz],
                radius=r,
                normal=[nx, ny, nz],
                spacing=spacing,
                flag=flag
            )
            fractures.append({
                "center": [cx, cy, cz],
                "radius": r,
                "normal": [nx, ny, nz],
                "spacing": spacing,
                "flag": flag,
                "points": pts
            })

    # -------------------------------------------------------------------------
    # SAVE FRACTURE DEFINITIONS TO JSON (to avoid re-prompting in the other script)
    # -------------------------------------------------------------------------
    def to_python_list(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (list, dict)):
            return obj
        return obj

    clean_fractures = []
    for frac in fractures:
        clean_fract = {
            "center": to_python_list(frac["center"]),
            "normal": to_python_list(frac["normal"]),
            "radius": float(frac["radius"]),
            "spacing": float(frac["spacing"]),
            "flag": int(frac["flag"])
        }
        clean_fractures.append(clean_fract)

    with open("fractures.json", "w") as f:
        json.dump(clean_fractures, f, indent=2)
    # -------------------------------------------------------------------------

    # Prompt for the .poly filename
    filename = poly_filename

    # Add fracture points => parse facets, remove near-facet points, then append
    add_fracture_points_to_poly(filename, fractures)
    print("------------------------------------------------------")

    # Run TetGen
    print("\n--> Running TetGen...")
    cmd = f"{tetgen_exe} -pnq1.3a1e12aAA {filename}"
    subprocess.run(cmd, shell=True)
    print("--> Done TetGen!")

if __name__ == "__main__":
    main()
#export PATH="/home/yao/Desktop/tetgen/tetgen1.6.0:$PATH"
