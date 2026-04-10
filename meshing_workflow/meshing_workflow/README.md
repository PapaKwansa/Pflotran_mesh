### Info ###

Meshing Workflow for 4100_Level Testbed Simulation

Authors: Chuanyao Zhong, Tim Johnson, Satish Karra. Pacific Northwest National Laboratory.

Thrust 2: Field-Scale Assimilation and Validation, Center for Understanding Subsurface Signals and Permeability (CUSSP)


### Main Features ###

1. This automated workflow generates high-resolution and real-world scale tetrahedral and Voronoi meshes, along 
	with all necessary input files (e.g., mesh files, boundary files, and material ID files) for PFLOTRAN simulations

2. Generate the mesh for any number of user-defined arbitrarily distributed fracture planes, and assign material IDs

3. Identify and correctly assign material IDs to cells located at fracture intersections

4. The location and size of borehole intervals (for simulating injection or pressurization operations) can be specified, 
	with corresponding region and material ID files generated
	
5. Identify and correctly assign material IDs to cells located at fracture-borehole intersections

6. Generate separate visualization files for each material in the Voronoi mesh


### Usage ###

1. Have TetGen installed. Please see: https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1

2. Have VORONOI installed. Please see: https://gitlab.pnnl.gov/pflotran/cussp/voronoi

3. Update your local paths for tetgen.exe and voronoi.exe in workflow.py

4. Run workflow.py

5. Define fractures (if desired) and borehole interval parameters by entering responses to the prompts

6. Copy and paste generated files into the designated folder for PFLOTRAN simulations