#include <iostream>
#include <vector>
#include <fstream>


int main()
{
    std::vector<int> Npr{2, 4, 8, 16, 32, 64};
    std::vector<int> Nthr{1, 2, 4};

    for (int p : Npr) {
        for (int t : Nthr) {
            int tasks_per_node = std::min(p, 32 / t);
            int nodes = p / tasks_per_node;
            std::string filename("test-qr-");
            filename += std::to_string(p) + "-" + std::to_string(t) + ".sh";
            std::ofstream out(filename);
            out << "#!/bin/bash\n\n\n";
            out << "#SBATCH --job-name=qr-rot-par-" << p << "-" << t << "\n";
            out << "#SBATCH --output=./results/log_qr_par_" << p << "-" << t << ".txt\n";
            out << "#SBATCH --partition=x20core\n";
            out << "#SBATCH --time=10:00:00\n";
            out << "#SBATCH --nodes=" << nodes << "\n";
            out << "#SBATCH --ntasks-per-node=" << tasks_per_node << "\n";
            out << "#SBATCH --cpus-per-task=" << t << "\n";
            out << "export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n";
            out << "mpirun ../test_qr_par\n";
            out.close();
        }
    }
}
