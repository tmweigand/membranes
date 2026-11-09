#!/bin/bash
#SBATCH --job-name=water100_prod
#SBATCH --output=water100_prod_%j.out
#SBATCH --error=water100_prod_%j.err
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=04:00:00

module load lammps/20240929-4

srun lmp -in in.water_100bar_prod
