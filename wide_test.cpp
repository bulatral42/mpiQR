#include <mpi.h>
#include <iostream>
#include <vector>

#include "qrdec.hpp"


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int block_size = 128;
    int subblock_size = 16;
    std::vector<int> Ns{512, 1024, 2048, 4096, 8192, 16384, 32768};

    int proc_id, n_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

    if (n_proc == 2) {
        Ns.pop_back();
    }
    if (n_proc == 64) {
        Ns.push_back(65536);
    }
    const int repeat = 3;
    for (int N : Ns) {
        if (n_proc > N / block_size) {
            continue;
        }
        for (int i = 1; i <= repeat; ++i) {
            if (proc_id == 0) {
                std::cout << "Test #" << i << std::endl;
            }
            block_qrdec(N, block_size, subblock_size);
        }
    }

    MPI_Finalize();
    return 0;
}

