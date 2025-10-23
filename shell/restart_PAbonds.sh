#!/bin/bash
#SBATCH --job-name=poly_rectangle
#SBATCH --output=lammps_output.out
#SBATCH --error=lammps_error.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --time=12:00:00
#SBATCH --partition=RM

module load openmpi/4.1.1-gcc8.3.1

mpirun -np 128 /jet/home/ajotcham/lammps-29Aug2024/build/lmp -in ../../domain_generation/restart_PA_bond.in -v mult 0.1 -v scale 1.0 -v xlink 0.90 -v rand 1 -v init_bonds 100
