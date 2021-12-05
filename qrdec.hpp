//
// Created by bulat_v on 05.12.2021.
//

#ifndef MPIQR_QRDEC_HPP
#define MPIQR_QRDEC_HPP

#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include "BLAS.h"
#include "LAPACK.h"

#include <unistd.h>


void block_qrdec(int N, int block_size, int subblock_size, MPI_Comm comm = MPI_COMM_WORLD)
{
    double t_start{0}, t_end{0}, mpi_qr_time{0};

    int n_proc, proc_id, n_thr;
    MPI_Comm_rank(comm, &proc_id);
    MPI_Comm_size(comm, &n_proc);
#pragma omp parallel
    {
        n_thr = omp_get_num_threads();
    }

    int n_block_cols = N / (n_proc * block_size), n_block_rows = N / block_size; // per proc

    if (proc_id == 0) {
        std::cout << "Start mpi-omp-QR, n_proc = " << n_proc << ", n_thr = " << n_thr << std::endl;
        std::cout << "N = " << N << ", k = " << block_size << ", m = " << subblock_size;
        std::cout << ", nbc = " << n_block_cols << ", nbr = " << n_block_rows << std::endl;
    }
    // blocks are stored rowwise, every block inside is stored columnwise
    double *blocks = new double[n_block_rows * n_block_cols * block_size * block_size];

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

    delete[] my_rotations;
    delete[] out_rotations;
    delete[] blocks;
}

void simple_qrdec(int N)
{
    std::cout << "Start simple-QR" << std::endl;
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
}

#endif //MPIQR_QRDEC_HPP
