#!/bin/bash
#SBATCH --job-name=membrane_persistent
#SBATCH --output=membrane_persistent.out
#SBATCH --error=membrane_persistent.err
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=64
#SBATCH --time=48:00:00
#SBATCH --partition=RM

module load anaconda3/2024.10-1
conda activate /jet/home/tweigand/my_env

mpirun -np 320 python rdf/persistent_connected_paths.py
