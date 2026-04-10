import sys

# Check if an argument (mesh_name) is provided
if len(sys.argv) < 2:
    raise ValueError("Missing mesh_name argument")

mesh_name = sys.argv[1]

# Define file paths using the mesh_name
materials_file = f"{mesh_name}_materials_updated.txt"
boundaries_file = f"{mesh_name}_boundaries.txt"
output_file = f"{mesh_name}_boundaries_updated.txt"

# Read the materials file
with open(materials_file, 'r') as mf:
    materials = mf.readlines()

# Read the boundaries file
with open(boundaries_file, 'r') as bf:
    boundaries = bf.readlines()

# Ensure the two files have the same number of lines
if len(materials) != len(boundaries):
    raise ValueError("The number of lines in the materials file does not match the number of lines in the boundaries file.")

# Update the boundaries file based on the conditions
updated_boundaries = []
for material_line, boundary_line in zip(materials, boundaries):
    material_value = material_line.strip()
    if material_value == '11':  # interval 
        updated_boundaries.append("13\n")  # Replace the boundary value with 8
    elif material_value == '6':  # TL
        updated_boundaries.append("9\n")
    elif material_value == '7':  # TS
        updated_boundaries.append("10\n")
    elif material_value == '8':  # TU
        updated_boundaries.append("11\n")
    elif material_value == '9':  # TN
        updated_boundaries.append("12\n")
    elif material_value == '10':  # TC
        updated_boundaries.append("8\n")
    elif material_value == '22':  # TL_CELL
        updated_boundaries.append("14\n")
    elif material_value == '23':  # TS_CELL
        updated_boundaries.append("15\n")
    elif material_value == '24':  # TU_CELL
        updated_boundaries.append("16\n")
    elif material_value == '25':  # TN_CELL
        updated_boundaries.append("17\n")
    elif material_value == '26':  # TC_CELL
        updated_boundaries.append("18\n")
    else:
        updated_boundaries.append(boundary_line)  # Keep the original boundary value

# Write the updated boundaries to a new file
with open(output_file, 'w') as of:
    of.writelines(updated_boundaries)

print(f"Updated boundaries file saved to {output_file}")
