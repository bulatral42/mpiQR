#include <mpi.h>
#include <iostream>
#include <cmath>
#include <random>
#include "BLAS.h"
#include "LAPACK.h"

#include <unistd.h>

#include "block_qrdec.hpp"


void print_matrix(int M, int N, double *A) {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << A[i + M * j] << " ";
        }
        std::cout << std::endl;
    }
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    if (argc < 4) {
        std::cout << "Usage: ./binary N k m" << std::endl;
        return -1;
    }

    std::istringstream sN(argv[1]);
    int N{0}; // Full matrix size
    sN >> N;

    std::istringstream sk(argv[2]);
    int block_size{0}; // number of cols/rows in block
    sk >> block_size;

    std::istringstream sm(argv[3]);
    int subblock_size{0}; // number of cols for rotation by 1 thread
    sm >> subblock_size;

    block_qrdec(N, block_size, subblock_size);

    MPI_Finalize();
    return 0;
}

