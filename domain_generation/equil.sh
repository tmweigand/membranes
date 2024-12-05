#!/bin/bash
#SBATCH --job-name=equilibration
#SBATCH --output=equil_output.out
#SBATCH --error=equil_error.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --time=48:00:00
#SBATCH --partition=RM

module load openmpi/4.1.1-gcc8.3.1
# mkdir logs
# mkdir restarts_polym

mpirun -np 512 /jet/home/ajotcham/lammps-29Aug2024/build/lmp -in ../../domain_generation/equilibrate.in