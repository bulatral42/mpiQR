#!/bin/bash


#SBATCH --job-name=qr-rot-1-1
#SBATCH --output=./results/log_qr_1-1.txt
#SBATCH --partition=x20core
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

module load gnu openblas openmpi

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun ../qr_seq
