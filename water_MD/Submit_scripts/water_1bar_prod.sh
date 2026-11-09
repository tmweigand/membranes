#!/bin/bash
#SBATCH --job-name=water1bar_prod
#SBATCH --output=water1bar_prod_%j.out
#SBATCH --error=water1bar_prod_%j.err
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=04:00:00

module load lammps/20240929-4

srun lmp -in Input/riley/in.water_1bar_prod
