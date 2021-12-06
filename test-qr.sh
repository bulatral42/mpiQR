#!/bin/bash


#SBATCH --job-name=qr-rot
#SBATCH --output=./log_qr.txt
#SBATCH --partition=x20core
#SBATCH --time=01:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2

module load mkl openmpi

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun ./qr_par
