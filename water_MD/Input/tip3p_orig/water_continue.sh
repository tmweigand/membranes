#!/bin/bash
#SBATCH --job-name=water_continue
#SBATCH --output=/hickory/users/r/g/rghanta/water_data/tip3p_orig/100bar/water_continue_%j.out
#SBATCH --error=/hickory/users/r/g/rghanta/water_data/tip3p_orig/100bar/water_continue_%j.err
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=02:00:00

module load lammps/20240929-4

echo "Job started: $(date)"
echo "Running on: $(hostname)"
which lmp

srun lmp -in in.water_continue

echo "Job finished: $(date)"
