MPI_CC = mpicxx
#OMPI_CXX=icpc mpicxx
FLAGS = -Wall -O2 -g -fopenmp
LIBS = -llapack -lblas


all: run

run: runner
	OMP_NUM_THREADS=4 mpirun -np 2 ./runner 16 4 2 #1024 128 16
#	salloc -p all -N 8 -c 1 /usr/bin/mpirun ./runner

runner: block_qrdec.cpp
	$(MPI_CC) $(FLAGS) block_qrdec.cpp -o runner $(LIBS)

clean:
	rm *.o ./runner

source:
	echo source /opt/intel/mkl/bin/mklvars.sh intel64

