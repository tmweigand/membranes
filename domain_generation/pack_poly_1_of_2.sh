#!/bin/bash
#SBATCH --job-name=64x_1_of_2
#SBATCH --output=lammps_output.out
#SBATCH --error=lammps_error.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --time=48:00:00
#SBATCH --partition=RM

module load openmpi/4.1.1-gcc8.3.1
mkdir logs
mkdir restarts_polym

mpirun -np 512 /jet/home/ajotcham/lammps-29Aug2024/build/lmp -in ../../domain_generation/rec_polymerize_1_of_2.in -var mult 64 -var scale 1.0 -var rand # -var xlink #
