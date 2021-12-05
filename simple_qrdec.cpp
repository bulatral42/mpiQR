#include <iostream>
#include "qrdec.hpp"


int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ./binary N" << std::endl;
        return -1;
    }

    std::istringstream sN(argv[1]);
    int N{0}; // Full matrix size
    sN >> N;

    std::cout << "N = " << N <<  std::endl;

    simple_qrdec(N);
    return 0;
}



