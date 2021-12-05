#include <mpi.h>
#include <iostream>
#include <vector>

#include "qrdec.hpp"


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    std::vector<int> Ns{8192, 16384, 32768, 65536, 131072};
    std::vector<int> block_sizes{64, 128, 256};
    std::vector<int> subblock_sizes{16, 32, 64, 128, 256};

    for (int N : Ns) {
        for (int block_size : block_sizes) {
            for (int subblock_size : subblock_sizes) {
                if (subblock_size > block_size) {
                    continue;
                }
                block_qrdec(N, block_size, subblock_size);
            }
        }
    }

    MPI_Finalize();
    return 0;
}

