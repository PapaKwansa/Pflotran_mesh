#!/bin/bash
#SBATCH --job-name=pflotran
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=250G
#SBATCH --time=48:00:00
#SBATCH --output=pflotran_%j.out
#SBATCH --error=pflotran_%j.err
# If Palmetto requires them, uncomment and fill these in:
##SBATCH --account=YOUR_ACCOUNT
##SBATCH --partition=YOUR_PARTITION

set -euo pipefail

module load pflotran

cd /home/harhin/Pflotran_mesh/avant_mesh
git pull

srun -n "$SLURM_NTASKS" /home/harhin/PFLOTRAN/petsc/pflotran/src/pflotran/pflotran -input_prefix layers4

python postprocess.py