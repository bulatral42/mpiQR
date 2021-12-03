MPI_CXX = mpicxx
CXX = g++
FLAGS_PAR = -Wall -O2 -fopenmp -std=c++11
FLAGS_SEQ = -Wall -O2 -std=c++11
LIBS = -llapack -lblas



all: qr_seq qr_par


run_seq: qr_seq
	./qr_seq 1024

run_par: qr_par
	OMP_NUM_THREADS=2 mpirun -np 4 ./qr_par 1024 128 16
#	salloc -p all -N 8 -c 1 /usr/bin/mpirun ./runner


qr_par: block_qrdec.cpp
	$(MPI_CXX) $(FLAGS_PAR) $(LIBS)  block_qrdec.cpp -o qr_par

qr_seq: simple_qrdec.cpp
	$(CXX) $(FLAGS_SEQ) $(LIBS)  simple_qrdec.cpp -o qr_seq


clean:
	rm ./qr_mpi ./qr_simp

