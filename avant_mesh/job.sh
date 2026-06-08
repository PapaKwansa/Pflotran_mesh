#!/bin/bash
#SBATCH --job-name=pflotran_geomech
#SBATCH --partition=work1
#SBATCH --nodes=2
#SBATCH --tasks-per-node=64
#SBATCH --mem=250G
#SBATCH --time=72:00:00
#SBATCH --output=pflotran_%j.out
#SBATCH --error=pflotran_%j.err

set -euo pipefail

module purge
module load spack
spack load pflotran@5.0.0

cd /scratch/$USER/avant_mesh_test

mpiexec -n $SLURM_NTASKS \
    pflotran \
    -input_prefix layers4_geomech

OUTDIR=$HOME/pflotran_results/${SLURM_JOB_ID}

mkdir -p $OUTDIR

rsync -av layers4_geomech* $OUTDIR/

echo "Results copied to:"
echo "$OUTDIR"