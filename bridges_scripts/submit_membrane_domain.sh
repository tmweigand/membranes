#!/bin/bash
#SBATCH --job-name=membrane_connected
#SBATCH --output=membrane_connected.out
#SBATCH --error=membrane_connected.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=64
#SBATCH --time=48:00:00
#SBATCH --partition=RM

module load anaconda3/2024.10-1
conda activate /jet/home/tweigand/my_env

mpirun -np 216 python rdf/membrane_domain.py
