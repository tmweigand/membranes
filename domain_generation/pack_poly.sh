#!/bin/bash
#SBATCH --job-name=lammps_job
#SBATCH --output=lammps_output.out
#SBATCH --error=lammps_error.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --time=12:00:00
#SBATCH --partition=RM

module load intelmpi/2021.3.0-intel2021.3.0 gcc/10.2.0 cuda/11.7.1 LAMMPS/23Jun22-intel

mpirun -np 128 lmp -in test_pack_polymerize.in -var mult 1 -var rand 1 -var xlink 0.93
