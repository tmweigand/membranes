#!/bin/bash
#SBATCH --job-name=water_test
#SBATCH --output=/hickory/users/r/g/rghanta/water_data/tip3p_orig/100bar/water_test_%j.out
#SBATCH --error=/hickory/users/r/g/rghanta/water_data/tip3p_orig/100bar/water_test_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:15:00
#SBATCH --partition=small

module load lammps/20240929-4

echo "Job started: $(date)"
echo "Node: $(hostname)"
which lmp

# Run LAMMPS using Slurm's task launcher
srun lmp -in in.water

echo "Job finished: $(date)"

