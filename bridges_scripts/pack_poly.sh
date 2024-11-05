#!/bin/bash
#SBATCH --job-name=lammps_job
#SBATCH --output=lammps_output.out
#SBATCH --error=lammps_error.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --time=12:00:00
#SBATCH --partition=RM

module load openmpi/4.1.1-gcc8.3.1
mkdir logs
mkdir restarts_polym

pwd
mpirun -np 128 /jet/home/ajotcham/lammps-29Aug2024/build/lmp -in  ../../domain_generation/pack_polymerize.in -var mult 1 -var rand 1 -var xlink 0.93
