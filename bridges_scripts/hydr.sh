#!/bin/bash
#SBATCH --job-name=lammps_job
#SBATCH --output=lammps_output.out
#SBATCH --error=lammps_error.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --partition=RM

module load openmpi/4.1.1-gcc8.3.1
mkdir restarts_hydr

mpirun -np 128 /jet/home/ajotcham/lammps-29Aug2024/build/lmp -in  ../../domain_generation/hydrate.in ###Need to add variable initializations here
