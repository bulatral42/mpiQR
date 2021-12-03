MPI_CC = mpicxx
FLAGS = -Wall -O2 -g -fopenmp
LIBS = -llapack -lblas



all: run

run: runner
	OMP_NUM_THREADS=2 mpirun -np 4 ./runner 1024 128 16
#	salloc -p all -N 8 -c 1 /usr/bin/mpirun ./runner

runner: block_qrdec.cpp
	$(MPI_CC) $(FLAGS) $(LIBS)  block_qrdec.cpp -o runner

clean:
	rm *.o ./runner

source:
	echo source /opt/intel/mkl/bin/mklvars.sh intel64

