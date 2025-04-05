#!/bin/bash
#SBATCH --job-name=psd_job
#SBATCH --output=psd_output.out
#SBATCH --error=psd_error.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=30:00:00
#SBATCH --partition=RM

module load anaconda3/2024.10-1
conda activate /jet/home/ajotcham/my_env/

mpirun -np 8 python pore_size_distributions/generate_psd.py
