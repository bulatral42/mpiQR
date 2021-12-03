#include <mpi.h>
#include <iostream>
#include <cmath>
#include <random>
#include "BLAS.h"
#include "LAPACK.h"

#include <unistd.h>


double get_elem_A(size_t i, size_t j) {
    return 1.0 / (i + j + 1);
}

double get_elem_B(size_t i, size_t j) {
    return sin(7 * i + 8 * j + 9);
}
double get_elem_C(size_t i, size_t j) {
    return 100 - std::abs((int)i - 9) - std::abs((int)j - 6);
}
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

    double t_start{0}, t_end{0}, mpi_qr_time{0};


    MPI_Comm comm = MPI_COMM_WORLD;
    int n_proc, proc_id;
    MPI_Comm_rank(comm, &proc_id);
    MPI_Comm_size(comm, &n_proc);

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

    int n_block_cols = N / (n_proc * block_size), n_block_rows = N / block_size; // per proc

    if (proc_id == 0) {
        std::cout << "N = " << N << ", k = " << block_size << ", m = " << subblock_size;
        std::cout << ", nbc = " << n_block_cols << ", nbr = " << n_block_rows << std::endl;
    }
    // blocks are stored rowwise, every block inside is stored columnwise
    double *blocks = new double[n_block_rows * n_block_cols * block_size * block_size];

    // TODO: fill matrix
    std::srand(proc_id);
    for (int i = 0; i < n_block_rows * n_block_cols * block_size * block_size; ++i) {
        //blocks[i] = i + proc_id * n_block_rows * n_block_cols * block_size * block_size;
        blocks[i] = -1.0 + 2 * static_cast<double>(std::rand()) / RAND_MAX;
    }

    // Main algorithm start
    t_start = MPI_Wtime();

    int proc_rots_cnt = (n_block_rows - proc_id) * n_block_cols - n_proc * ((n_block_cols - 1) * n_block_cols) / 2;
    proc_rots_cnt *= block_size * block_size;
    proc_rots_cnt -= n_block_cols * (block_size * (block_size + 1)) / 2;
    //std::cout << "id = " << proc_id << ": rots_cnt = " << proc_rots_cnt << std::endl;
    double *my_rotations = new double[2 * proc_rots_cnt];
    double *out_rotations = new double[2 * block_size * block_size];
    double *cur_rots = my_rotations;

    // Main cycle
    for (int b_col = 0; b_col < n_block_cols * n_proc; ++b_col) {
        int sender_id = b_col % n_proc;
        int local_b_col = b_col / n_proc;
        double *main_block = blocks + (block_size * block_size) * (b_col * n_block_cols + local_b_col);
        for (int b_row = b_col; b_row < n_block_rows; ++b_row) {
            double *cur_block = blocks + (block_size * block_size) * (b_row * n_block_cols + local_b_col);
            int n_rots{0};
            bool triang_zero = b_row == b_col; // Zerorize triangular part or full block
            if (proc_id == sender_id) { // Sender part
                for (int j = 0; j < block_size; ++j) {
                    for (int i = triang_zero ? j + 1 : 0; i < block_size; ++i) {
                        double r{0.0};
                        dlartg(main_block[j * block_size + j], cur_block[j * block_size + i],
                               cur_rots[2 * n_rots], cur_rots[2 * n_rots + 1], r);
                        drot(block_size - j - 1, main_block + (j + 1) * block_size + j, block_size, 
                             cur_block + (j + 1) * block_size + i, block_size, 
                             cur_rots[2 * n_rots], cur_rots[2 * n_rots + 1]);
                        main_block[j * block_size + j] = r;
                        cur_block[j * block_size + i] = 0.0;
                        ++n_rots;
                    }
                }
            } else {
                if (triang_zero) {
                    n_rots = ((block_size - 1) * block_size) / 2;
                } else {
                    n_rots = block_size * block_size;
                }
            }
            // Send/receive rotations
            MPI_Request req_bcast;
            double *senrecv_rots = proc_id == sender_id ? cur_rots : out_rotations;
            MPI_Ibcast(senrecv_rots, 2 * n_rots, MPI_DOUBLE, sender_id, comm, &req_bcast);
                
            // Perform rotations to other blocks in row
            double *up_other_blocks = main_block;
            double *dn_other_blocks = cur_block;
            int n_subblocks = (block_size / subblock_size) * (n_block_cols - b_col);
            if (proc_id == sender_id) { // Current block has already been rotated
                up_other_blocks += block_size * block_size;
                dn_other_blocks += block_size * block_size;
                n_subblocks -= block_size / subblock_size;
            } else {
                MPI_Status status;
                MPI_Wait(&req_bcast, &status);
            }
#pragma omp parallel for schedule(static)
            for (int subbl_id = 0; subbl_id < n_subblocks; ++subbl_id) {
                double *up_subblock = up_other_blocks + block_size * subblock_size * subbl_id;
                double *dn_subblock = dn_other_blocks + block_size * subblock_size * subbl_id;
                int up_off{0};
                int dn_off = triang_zero ? up_off + 1 : 0;
                for (int i = 0; i < n_rots; ++i) {
                    drot(subblock_size, up_subblock + up_off, block_size, dn_subblock + dn_off, 
                            block_size, senrecv_rots[2 * i], senrecv_rots[2 * i + 1]);
                    ++dn_off;
                    if (dn_off % block_size == 0) {
                        ++up_off;
                        dn_off = triang_zero ? up_off + 1 : 0;
                    }
                }   
            }
            if (proc_id == sender_id) { // Current block has already been rotated
                MPI_Status status;
                MPI_Wait(&req_bcast, &status);
                cur_rots += 2 * n_rots;
            }
        }
    }
    MPI_Barrier(comm);
    t_end = MPI_Wtime();
    mpi_qr_time = t_end - t_start;
    if (proc_id == 0) {
        std::cout << "MPI QR time = " << mpi_qr_time << " sec" << std::endl;
    }
    //std::cout << "id = " << proc_id << ", blocks = \n";
    //print_matrix(block_size, n_block_cols * n_block_rows * block_size, blocks);
    delete[] my_rotations;
    delete[] out_rotations;
    delete[] blocks;

    MPI_Finalize();
    return 0;
}

