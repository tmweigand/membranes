#!/bin/bash
#SBATCH --job-name=membrane_connected_profiling
#SBATCH --output=membrane_connected_profiling.out
#SBATCH --error=membrane_connected_profiling.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --time=10:00:00
#SBATCH --partition=RM

module load anaconda3/2024.10-1
conda activate /jet/home/tweigand/my_env

mpirun -np 8 python rdf/generate_rdf.py
mpirun -np 27 python rdf/generate_rdf.py
mpirun -np 64 python rdf/generate_rdf.py
mpirun -np 125 python rdf/generate_rdf.py
