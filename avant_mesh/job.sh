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

# Stage 1: injection
/home/harhin/PFLOTRAN/petsc/arch-linux-c-opt/bin/mpirun -np 64 \
  /home/harhin/PFLOTRAN/petsc/pflotran/src/pflotran/pflotran \
  -input_prefix layers4_inj

# Stage 2: recovery
/home/harhin/PFLOTRAN/petsc/arch-linux-c-opt/bin/mpirun -np 64 \
  /home/harhin/PFLOTRAN/petsc/pflotran/src/pflotran/pflotran \
  -input_prefix layers4_recovery