#!/usr/bin/env python3

"""
Replicate Tim's codes, and generalize it for the new workflow.
"""

import os
import numpy as np
from scipy.spatial import Voronoi
import vtk
from vtkmodules.vtkCommonCore import vtkIdList, vtkPoints
from vtkmodules.vtkCommonDataModel import VTK_POLYHEDRON
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridWriter

print(f"\nGenerating Voronoi visualization files...")

# Ask the user for the prefix
pre = input("\nEnter the prefix for the mesh files ('surf2' in this case): ").strip()

# Read the node file (Voronoi cell centers)
tnodes = np.genfromtxt(pre + '.1.node', skip_header=1, usecols=(1,2,3))

# Read material IDs from _materials.txt
matids = np.genfromtxt(pre + '_materials.txt', dtype=int)

# Define material names mapping
material_names = {
    2: "AMU", 3: "AML", 4: "DMU", 5: "DML",
    6: "TL", 7: "TS", 8: "TU", 9: "TN", 10: "TC", 
    11: "well_interval"
}

# Identify fracture materials (>= 11) and name them dynamically
fracture_ids = np.unique(matids[matids >= 12])
fracture_map = {fid: f"Fracture{i+1}" for i, fid in enumerate(sorted(fracture_ids))}

# Merge fracture map into material names
material_names.update(fracture_map)

# Skip host rock (ID=1)
materials_to_process = sorted(set(matids) - {1})

# Build the Voronoi diagram
vor = Voronoi(tnodes)

def map_used_vertices(m_ids, m_keep, vor):
    """
    Maps Voronoi vertices used in visualization.
    
    - `m_ids`: Material IDs for each node.
    - `m_keep`: Material IDs to keep in visualization.
    - `vor`: Voronoi object from SciPy.
    
    Returns:
        - `v_map`: Mapping of Voronoi vertices to new indices (or -999 if unused).
    """
    v_keep = np.zeros(len(vor.vertices), dtype=bool)
    fc_keep = np.zeros(len(vor.ridge_vertices), dtype=bool)
    
    for i, m in enumerate(m_ids):
        if m in m_keep:
            n1 = np.where(vor.ridge_points[:, 0] == i)[0]
            n2 = np.where(vor.ridge_points[:, 1] == i)[0]
            for f in n1:
                fc_keep[f] = True
                for v in vor.ridge_vertices[f]:
                    v_keep[v] = True
            for f in n2:
                fc_keep[f] = True
                for v in vor.ridge_vertices[f]:
                    v_keep[v] = True
    
    v_map = np.full(len(vor.vertices), -999, dtype=int)
    count = 0
    for i, v in enumerate(v_keep):
        if v:
            v_map[i] = count
            count += 1
    
    return v_map

# Follow the same Tim's original logic for processing
for mt in materials_to_process:
    if mt not in material_names:
        material_name = f"Material{mt}"
    else:
        material_name = material_names[mt]

    print(f"\nProcessing Material {mt}: {material_name}")

    # Map vertices used for the current material
    vmap = map_used_vertices(matids, [mt], vor)

    # Create VTK Points object
    points = vtkPoints()
    for i, p in enumerate(vor.vertices):
        if vmap[i] >= 0:
            points.InsertNextPoint(p[0], p[1], p[2])

    # Initialize the VTK unstructured grid and add points
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)

    # Add Voronoi cells to the grid
    for i, flg in enumerate(matids):
        if flg == mt:
            f1 = np.where(vor.ridge_points[:, 0] == i)[0]
            f2 = np.where(vor.ridge_points[:, 1] == i)[0]
            faces = []
            for f in f1:
                faces.append(vmap[vor.ridge_vertices[f]])
            for f in f2:
                faces.append(vmap[vor.ridge_vertices[f]])

            faceId = vtkIdList()
            faceId.InsertNextId(len(faces))  # Number of faces
            for face in faces:
                faceId.InsertNextId(len(face))  # Number of points per face
                [faceId.InsertNextId(j) for j in face]
                
            grid.InsertNextCell(VTK_POLYHEDRON, faceId)

    # Add material ID as cell data
    cellData = vtk.vtkDoubleArray()
    cellData.SetName("Mat_IDS")
    for i, flg in enumerate(matids):
        if flg == mt:
            cellData.InsertNextValue(flg)
    grid.GetCellData().AddArray(cellData)

    # Add elevation data
    cellData = vtk.vtkDoubleArray()
    cellData.SetName("Elev")
    for i, flg in enumerate(matids):
        if flg == mt:
            cellData.InsertNextValue(tnodes[i, 2])
    grid.GetCellData().AddArray(cellData)

    # Write the .vtu file with the appropriate name
    vtu_filename = f"{pre}_{material_name}.vtu"
    writer = vtkXMLUnstructuredGridWriter()
    writer.SetInputData(grid)
    writer.SetFileName(vtu_filename)
    writer.Update()
    writer.Write()
    
    print(f"\nFinished writing {vtu_filename}")

print("\nAll .vtu files (excluding host rock) have been generated.")
