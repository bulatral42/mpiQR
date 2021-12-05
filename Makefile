MPI_CXX = mpicxx
FLAGS= -Wall -O2 -fopenmp -std=c++11
LIBS = -llapack -lblas


all: qr_seq qr_par


run_seq: qr_seq
	./qr_seq 1024

run_par: qr_par
	OMP_NUM_THREADS=2 mpirun -np 4 ./qr_par 8192 32 16


qr_par: block_qrdec.cpp qrdec.hpp
	$(MPI_CXX) $(FLAGS) block_qrdec.cpp -o qr_par $(LIBS)

qr_seq: simple_qrdec.cpp qrdec.hpp
	$(MPI_CXX) $(FLAGS) simple_qrdec.cpp -o qr_seq $(LIBS)


clean:
	rm -f ./qr_seq ./qr_par

