#!/bin/bash
#SBATCH --job-name=scaling_hydr_1
#SBATCH --output=scaling_hydr_1_output.out
#SBATCH --error=scaling_hydr_1_error.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=02:00:00
#SBATCH --partition=RM

module load openmpi/4.1.1-gcc8.3.1
# mkdir logs
# mkdir restarts_polym

mpirun -np 128 /jet/home/ajotcham/lammps-29Aug2024/build/lmp -in ../../domain_generation/hydrate.in -var mult 64 -var rand 2  -var FEEDP 1