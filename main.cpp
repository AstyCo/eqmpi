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
        times.clear();
        int N = Ns[i];
        Profiler p_finalization;
        {

            Profiler p;
            Iterations its(N); // iterations parameters, send/recv buffers
            its.prepare();
            MY_ASSERT(cnode.mpi.procCount == 1);
            its.seqRun();

            MPI_Barrier(MPI_COMM_WORLD);
            if (cnode.mpi.rank == 0) {
                int nthread = 1;
                std::cout << SSTR("###," << cnode.scTag()
                                  << ',' << cnode.mpi.procCount
                                  << ',' << nthread
                                  << ',' << N
                                  << ',' << p.time() ) << std::endl;
            }
            p_finalization.start();
        }
        get_time(times.program_finalization, p_finalization);
        cnode.print0(times.get_times(N));
    }


    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
