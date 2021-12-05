#include <iostream>
#include "qrdec.hpp"


int main(int argc, char *argv[]) {

    std::vector<int> Ns{2048, 4096, 8192, 16384, 32768, 65536};

    for (int N : Ns) {
        simple_qrdec(N);
    }
    return 0;
}



