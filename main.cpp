#include "iterations.hpp"

int main(int argc, char **argv)
{
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    clargs.parse(argc, argv);

    cnode.init(); // rank, size

    Profiler prf;
    prf.start();
    Iterations its; // iterations parameters, send/recv buffers

    its.prepare();
    its.run();

    MPI_Barrier(MPI_COMM_WORLD);
    cnode.print("FINALIZE");
    cnode.print(SSTR("TIME " << prf.time() << "s."));
    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
