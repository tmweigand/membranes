#!/bin/bash
#SBATCH --job-name=scaling_hydr_8
#SBATCH --output=scaling_hydr_8_output.out
#SBATCH --error=scaling_hydr_8_error.err
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=128
#SBATCH --time=02:00:00
#SBATCH --partition=RM

module load openmpi/4.1.1-gcc8.3.1
# mkdir logs
# mkdir restarts_polym

mpirun -np 1024 /jet/home/ajotcham/lammps-29Aug2024/build/lmp -in ../../domain_generation/hydrate.in -var mult 64 -var rand 2  -var FEEDP 1