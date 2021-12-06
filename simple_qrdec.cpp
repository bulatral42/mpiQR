#include <iostream>
#include "qrdec.hpp"


int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    std::vector<int> Ns{256, 512, 1024};

    for (int N : Ns) {
        simple_qrdec(N);
    }
    MPI_Finalize();
    return 0;
}



