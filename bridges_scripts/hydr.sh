#!/bin/bash
#SBATCH --job-name=hydration
#SBATCH --output=hydration_output_%j.out
#SBATCH --error=hydration_error_%j.err
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=128
#SBATCH --time=48:00:00
#SBATCH --partition=RM

module load openmpi/4.1.1-gcc8.3.1
mkdir hydr
mkdir hydr/pressure
mkdir hydr/system
mkdir hydr/membrane
mkdir hydr/restarts

mpirun -np 1024 /jet/home/ajotcham/lammps-29Aug2024/build/lmp -in ../../domain_generation/hydrate.in -var mult 64 -var rand 5  -var FEEDP 1