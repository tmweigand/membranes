#!/bin/bash
#SBATCH --job-name=water5000_prod
#SBATCH --output=water5000_prod_%j.out
#SBATCH --error=water5000_prod_%j.err
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=04:00:00

module load lammps/20240929-4

srun lmp -in in.water_5000bar_prod

