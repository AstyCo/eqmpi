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
#ifdef BGP
        Ns.push_back(128);
        Ns.push_back(256);
        Ns.push_back(512);
#endif
#ifdef POLUS
        Ns.push_back(128);
        Ns.push_back(256);
        Ns.push_back(512);
#endif
    }

    for (uint i = 0; i < Ns.size(); ++i) {
        int N = Ns[i];

        Profiler prf;

        Iterations its(N); // iterations parameters, send/recv buffers
        its.prepare();
        its.run();

        MPI_Barrier(MPI_COMM_WORLD);
        prf.finish();
        if (cnode.mpi.rank == 0) {
            int nthread = 1;
            std::string sc("U");
#ifdef BGP
            sc = "BGP";
#endif
#ifdef POLUS
            sc = "POLUS";
#endif
#ifdef WITH_OMP
            nthread = 3;
#endif
            std::cout << SSTR("###," << sc
                              << ',' << cnode.mpi.procCount
                              << ',' << nthread
                              << ',' << N
                              << ',' << prf.time() ) << std::endl;
        }
//        cnode.print(SSTR("SIZE " << N << " TIME " << prf.time()));
    }

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
