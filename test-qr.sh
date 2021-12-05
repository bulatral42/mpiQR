#!/bin/bash

#SBATCH --job-name=qr-rot
#SBATCH --output=./log_qr_rot.txt
#SBATCH --partition=test
#SBATCH --time=00:15:00

####### MPI ranks

#SBATCH --nodes=8

####### 6 OMP threads per MPI rank
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=14


module load mkl openmpi


## set number of OMP thread to the value of --cpus-per-task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpirun ./qr_par 2048 64 16
