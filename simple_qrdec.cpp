#include <iostream>
#include <chrono>
#include "BLAS.h"
#include "LAPACK.h"


int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ./binary N" << std::endl;
        return -1;
    }

    std::istringstream sN(argv[1]);
    int N{0}; // Full matrix size
    sN >> N;

    std::cout << "N = " << N <<  std::endl;
    
    double *block = new double[N * N];

    std::srand(13);
    for (int i = 0; i < N * N; ++i) {
        //block[i] = i;
        block[i] = -1.0 + 2 * static_cast<double>(std::rand()) / RAND_MAX;
    }

    // Main algorithm start
    auto t_start = std::chrono::steady_clock::now();

    int rots_cnt = (N * (N - 1)) / 2;
    double *rotations = new double[2 * rots_cnt];
    double *cur_rots = rotations;

    int n_rots = 0;
    // Main cycle
    for (int col = 0; col < N; ++col) {
        for (int row = col + 1; row < N; ++row) {
            double r{0.0};
            dlartg(block[col * N + col], block[col * N + row], 
                   cur_rots[2 * n_rots], cur_rots[2 * n_rots + 1], r);
            drot(N - col - 1, block + (col + 1) * N + col, N, block + (col + 1) * N + row, N, 
                 cur_rots[2 * n_rots], cur_rots[2 * n_rots + 1]);
            block[col * N + col] = r;
            block[col * N + row] = 0.0;
            ++n_rots;
        }
    }
   
    auto t_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> t_delta = t_end - t_start;
    std::cout << "Simple QR time = " << t_delta.count() << " sec" << std::endl;
 
    delete[] rotations;
    delete[] block;
    
    return 0;
}



