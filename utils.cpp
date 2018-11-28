#include "utils.hpp"

ComputeNode::ComputeNode()
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

	// Get the number of processes`
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    mpi.rank = world_rank;
    mpi.procCount = world_rank;

    fillGridDimensions();
    fillXYZ();
}

ComputeNode::ComputeNode()
{
	// Finalize the MPI environment.
    MPI_Finalize();
}

int ComputeNode::neighbor(ConnectionDirection cdir)
{
	switch (cdir) {
	case DIR_X: return 
	}
}

void ComputeNode::fillGridDimensions()
{
    gridDimensions = getGridDimensions()[mpi.procCount];
}

void ComputeNode::fillXYZ()
{
    x = mpi.rank % gridDimensions.x;
    uint yzRank = mpi.rank / gridDimensions.x;
    y = yzRank % gridDimensions.y;
    z = yzRank / gridDimensions.y;
}