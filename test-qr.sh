#!/bin/bash

module load mkl openmpi

#SBATCH --job-name=qr-rot-seq
#SBATCH --output=./log_qr_seq.txt
#SBATCH --partition=x20core
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun ./qr_seq


#SBATCH --job-name=qr-rot-par-1
#SBATCH --output=./log_qr_par_1.txt
#SBATCH --partition=x20core
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=2
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun ./test_qr_par

#SBATCH --job-name=qr-rot-par-2
#SBATCH --output=./log_qr_par_2.txt
#SBATCH --partition=x20core
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=4
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun ./test_qr_par

#SBATCH --job-name=qr-rot-par-3
#SBATCH --output=./log_qr_par_3.txt
#SBATCH --partition=x20core
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=2
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun ./test_qr_par
