#!/bin/bash
#SBATCH --job-name=scaling_poly_16
#SBATCH --output=scaling_poly_16_output.out
#SBATCH --error=scaling_poly_16_error.err
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=128
#SBATCH --time=02:00:00
#SBATCH --partition=RM

module load openmpi/4.1.1-gcc8.3.1
mkdir logs
mkdir restarts_polym

mpirun -np 2048 /jet/home/ajotcham/lammps-29Aug2024/build/lmp -in ../../domain_generation/rec_polymerize_1_of_2.in -var mult 64 -var scale 1.0 -var rand 5 -var xlink 0.90