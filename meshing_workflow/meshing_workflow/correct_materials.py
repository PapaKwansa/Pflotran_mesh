import sys

# Check if an argument (mesh_name) is provided
if len(sys.argv) < 2:
    raise ValueError("Missing mesh_name argument")

mesh_name = sys.argv[1]

# Define input and output file paths
materials_file = f"{mesh_name}_materials.txt"
updated_materials_file = f"{mesh_name}_materials_updated.txt"

# Read the materials file
with open(materials_file, 'r') as mf:
    materials = mf.readlines()

# Define a dictionary to track replacements
replacements = {
    "6": "22",
    "7": "23",
    "8": "24",
    "9": "25",
    "10": "26"
}

# Keep track of which replacements have been made
replaced_flags = {key: False for key in replacements.keys()}

# Update the materials file
updated_materials = []
for line in materials:
    material_value = line.strip()
    if material_value in replacements and not replaced_flags[material_value]:
        updated_materials.append(replacements[material_value] + "\n")
        replaced_flags[material_value] = True  # Mark this value as replaced
    else:
        updated_materials.append(line)  # Keep the original value

# Write the updated materials to a new file
with open(updated_materials_file, 'w') as umf:
    umf.writelines(updated_materials)

print(f"Updated materials file saved to {updated_materials_file}")
