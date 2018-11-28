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

int ComputeNode::neighbor(ConnectionDirection cdir) const
{
	switch (cdir) {
	case DIR_X: return toRank(x + 1, y, z);
	case DIR_MINUS_X: return toRank(x - 1, y, z);
	case DIR_Y: return toRank(x, y + 1, z);
	case DIR_MINUS_Y: return toRank(x, y - 1, z);
	case DIR_Z: return toRank(x, y, z + 1);
	case DIR_MINUS_Z: return toRank(x, y, z - 1);
	default:
		MY_ASSERT(false);
		return -1;
	}
}

bool is(ConnectionDirection cdir) const
{
	return neighbor(cdir) != -1;
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

int ComputeNode::toRank(uint i, uint j, uint k) const
{
	if (i >= gridDimensions.x
		|| j >= gridDimensions.y
		|| k >= gridDimensions.z) {
		return -1;
	}
	return (i * gridDimensions.y + j) * gridDimensions.z + k;

}
