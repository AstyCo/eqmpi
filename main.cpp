#include "iterations.hpp"

int main(int argc, char **argv)
{
    clargs.parse(argc, argv);

    // ComputeNode RAII MPI resources (MPI_Initialize, MPI_Finalize)
    ComputeNode cnode; // rank, size

    Iterations its(cnode); // iterations parameters, send/recv buffers

    its.prepare();
    its.run();


    return 0;
}
