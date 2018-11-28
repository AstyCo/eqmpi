#ifndef UTILS_HPP
#define UTILS_HPP

#include "globals.hpp"

struct ComputeNode
{
    // MPI
    struct MPI_data
    {
        uint rank;
//        uint commSize;
        uint procCount;
    } mpi;


    GridDimensions gridDimensions;
    // Equation
    int x;
    int y;
    int z;

    ComputeNode();
    ~ComputeNode();

    int neighbor(ConnectionDirection cdir);
    void fillGridDimensions();
    void fillXYZ();
};


#endif // UTILS_HPP