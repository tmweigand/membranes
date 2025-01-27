#!/bin/bash
#SBATCH --job-name=hydration_restart
#SBATCH --output=hydration_restart_output_%j.out
#SBATCH --error=hydration_restart_error_%j.err
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=128
#SBATCH --time=48:00:00
#SBATCH --partition=RM

module load openmpi/4.1.1-gcc8.3.1
# mkdir hydr
# mkdir hydr/pressure
# mkdir hydr/system
# mkdir hydr/membrane
# mkdir hydr/restarts

mpirun -np 1024 /jet/home/ajotcham/lammps-29Aug2024/build/lmp -in ../../domain_generation/hydrate_restart.in -var mult 64 -var rand 3 -var FEEDP 1