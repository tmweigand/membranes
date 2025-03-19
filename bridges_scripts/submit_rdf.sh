#!/bin/bash
#SBATCH --job-name=rdf_job
#SBATCH --output=rdf_output.out
#SBATCH --error=rdf_error.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=30:00:00
#SBATCH --partition=RM

mpirun -np 8 python  ../rdf/generate_rdf.py
