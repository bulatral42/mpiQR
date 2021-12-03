#!/bin/bash

#SBATCH --job-name=qr-rot
#SBATCH --output=log_qr_rot.txt
#SBATCH --partition=queue
#SBATCH --time=05:00:00

####### MPI ranks

#SBATCH --nodes=16

####### 6 OMP threads per MPI rank
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20


###########module load gnu openmpi scalapack


## set number of OMP thread to the value of --cpus-per-task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpirun ./qr_par 8192 128 32
