#ifndef UTILS_HPP
#define UTILS_HPP

#include "globals.hpp"

struct ComputeNode
{
    // MPI
    struct MPI_data
    {
        uint rank;
        uint procCount;
    } mpi;


    GridDimensions gridDimensions;
    // Equation
    int x;
    int y;
    int z;

    ComputeNode();
    ~ComputeNode();

    int neighbor(ConnectionDirection cdir) const;
    bool is(ConnectionDirection cdir) const;
    void fillGridDimensions();
    void fillXYZ();

    int toRank(uint i, uint j, uint k) const;
};


#endif // UTILS_HPP