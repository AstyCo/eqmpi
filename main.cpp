#include "iterations.hpp"

int main(int argc, char **argv)
{
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    clargs.parse(argc, argv);

    cnode.init(); // rank, size

    Iterations its; // iterations parameters, send/recv buffers

    its.prepare();
    its.run();

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
