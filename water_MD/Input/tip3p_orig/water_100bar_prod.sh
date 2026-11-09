#!/bin/bash
#SBATCH --job-name=water100_prod
#SBATCH --output=/hickory/users/r/g/rghanta/water_data/tip3p_orig/100bar/water100_prod_%j.out
#SBATCH --error=/hickory/users/r/g/rghanta/water_data/tip3p_orig/100bar/water100_prod_%j.err
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=04:00:00

module load lammps/20240929-4

srun lmp -in in.water_100bar_prod
