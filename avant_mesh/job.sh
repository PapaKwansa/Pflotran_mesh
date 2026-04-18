#!/bin/bash
#SBATCH --job-name=pflotran
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=250G
#SBATCH --time=48:00:00
#SBATCH --output=pflotran_%j.out
#SBATCH --error=pflotran_%j.err

set -euo pipefail

cd /home/harhin/Pflotran_mesh/avant_mesh
srun -n 64 /home/harhin/PFLOTRAN/petsc/pflotran/src/pflotran/pflotran -input_prefix layers4
