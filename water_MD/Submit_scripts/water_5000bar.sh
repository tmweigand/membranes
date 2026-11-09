#!/bin/bash
#SBATCH --job-name=water_5000bar
#SBATCH --output=water_5000bar_%j.out
#SBATCH --error=water_5000bar_%j.err
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=02:00:00

module load lammps/20240929-4

echo "Job started: $(date)"
echo "Running on: $(hostname)"
which lmp

srun lmp -in in.water_5000bar

echo "Job finished: $(date)"
