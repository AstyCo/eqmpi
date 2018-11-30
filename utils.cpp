#include "utils.hpp"

#include <iostream>
#include <limits>

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>

CommandLineArgs clargs;

void initialize_seed()
{
    srand(time(NULL));
}

double randomize(double min, double max)
{
    return min + (static_cast<double>(rand()) / RAND_MAX) * (max - min);
}

ComputeNode::ComputeNode()
{
    std::cout << "DEBUG: ComputeNode\n";
    // Initialize the MPI environment
    MPI_Init(&clargs.argc, &clargs.argv);

#ifdef WITH_OMP
    omp_set_num_threads(3); // set the number of threads for this programm
    omp_set_dynamic(0); // allways use maximum number of threads (not less)
#endif

	// Get the number of processes`
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    mpi.rank = world_rank;
    mpi.procCount = world_size;

    fillGridDimensions();
    fillXYZ();

    std::cout << "DEBUG: ComputeNode +\n";
}

ComputeNode::~ComputeNode()
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

bool ComputeNode::is(ConnectionDirection cdir) const
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

class ProfilerPrivate
{
public:
    ProfilerPrivate(Profiler::Options opts)
        : _started(false),
          _wall_clock_elapsed(0), _opts(opts)
    {

    }

    ~ProfilerPrivate()
    {
        if (_opts & Profiler::PrintOnDestructor) {
            if (_started)
                finish();
            print();
        }
    }

    void start()
    {
        if (_started)
            clear();
        _started = true;

        _wstart = MPI_Wtime();
    }

    void clear()
    {
        _started = false;
    }

    void finish()
    {
        MY_ASSERT(_started);

        _wall_clock_elapsed = MPI_Wtime() - _wstart;

        clear();
    }

    void print() const
    {
        char buff[256];
        snprintf(buff, sizeof(buff),
                 "WC: %lf s.\n",
                 _wall_clock_elapsed);

        std::cout << buff << std::endl;
    }

    double time() const
    {
        return _wall_clock_elapsed;
    }

private:
    bool _started;

    double _wstart;

    double _wall_clock_elapsed;

    Profiler::Options _opts;
};

Profiler::Profiler(Profiler::Options opts)
    : _impl(new ProfilerPrivate(opts))
{

}

void Profiler::start()
{
    _impl->start();
}

void Profiler::finish()
{
    _impl->finish();
}

double Profiler::time() const
{
    return _impl->time();
}

void Profiler::print() const
{
    _impl->print();
}

void CommandLineArgs::parse(int argc_, char **argv_)
{
	argc = argc_;
	argv = argv_;
    for (int i = 1; i < argc; ++i)
        parseArg(argv[i]);
}

void CommandLineArgs::parseArg(char arg[])
{
    const char s_K [] = "K=";
    const char s_N [] = "N=";

    std::string sarg(arg);
    if (sarg.find(s_K) != std::string::npos) {
        long tmp = strtol(sarg.c_str() + sizeof(s_K) - 1,
                         NULL, 10);
        K = tmp;
        return;
    }
    if (sarg.find(s_N) != std::string::npos) {
        long tmp = strtol(sarg.c_str() + sizeof(s_N) - 1,
                         NULL, 10);
        N = tmp;
        return;
    }

    std::cout << "unrecognized argument " << sarg << std::endl;
}

