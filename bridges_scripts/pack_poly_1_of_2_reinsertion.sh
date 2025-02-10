#!/bin/bash
#SBATCH --job-name=64x_1.5_0.95_1_of_2_reinsertion
#SBATCH --output=poly_1_of_2_output_%j.out
#SBATCH --error=poly_1_of_2_error_%j.err
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=128
#SBATCH --time=48:00:00
#SBATCH --partition=RM

module load openmpi/4.1.1-gcc8.3.1
#mkdir logs
#mkdir restarts_polym

mpirun -np 1024 /jet/home/ajotcham/lammps-29Aug2024/build/lmp -in ../../domain_generation/rec_polymerize_1_of_2_reinsertion.in -var mult 64 -var scale 1.5 -var rand 3 -var xlink 0.95
