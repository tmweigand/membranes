#!/bin/bash
#SBATCH --job-name=water
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=04:00:00
#SBATCH --output=water_%j.out
#SBATCH --error=water_%j.err

set -euo pipefail

FORCEFIELD="riley"     # Give either "riley" or "tip3p_orig" as forcefield choice
PRESSURE=100           # Pressure in bar, other parameters can be modified in ../Input/water.py

module load lammps/20240929-4

python3 ../Input/water.py --forcefield "$FORCEFIELD" --pressure "$PRESSURE" --output "../Output/${FORCEFIELD}/${PRESSURE}bar"

echo "Simulation completed"
