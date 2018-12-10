#include "iterations.hpp"

int main(int argc, char **argv)
{
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    clargs.parse(argc, argv);

    cnode.init(); // rank, size

    std::vector<int> Ns;
    if (clargs.N > 0)
        Ns.push_back(clargs.N);
    if (Ns.empty()) {
        Ns.push_back(128);
        Ns.push_back(256);
        Ns.push_back(512);
    }

    for (uint i = 0; i < Ns.size(); ++i) {
        int N = Ns[i];

        profiler.start();
        Iterations its(N); // iterations parameters, send/recv buffers
        its.prepare();
        its.run();

        MPI_Barrier(MPI_COMM_WORLD);
        profiler.finish();
        if (cnode.mpi.rank == 0) {
            int nthread = 1;
            std::cout << SSTR("###," << cnode.scTag()
                              << ',' << cnode.mpi.procCount
                              << ',' << nthread
                              << ',' << N
                              << ',' << profiler.time() ) << std::endl;
        }
    }

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
