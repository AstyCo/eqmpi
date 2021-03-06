#include "utils.hpp"

#include <iostream>
#include <limits>

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>

Profiler profiler;
CommandLineArgs clargs;
ComputeNode cnode;

void initialize_seed()
{
    srand(time(NULL));
}

double randomize(double min, double max)
{
    return min + (static_cast<double>(rand()) / RAND_MAX) * (max - min);
}

std::string ComputeNode::scTag() const
{
    switch (sc) {
    case SCPolus: return "POL";
    case SCBluegeneP: return "BGP";
    }
    MY_ASSERT(false);
    return std::string();
}

void ComputeNode::init()
{
#ifdef WITH_OMP
    omp_set_num_threads(3); // set the number of threads for this programm (BlueGene/P)
#endif

#ifdef BGP
    sc = SCBluegeneP;
#endif
#ifdef POLUS
    sc = SCPolus;
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
}

int ComputeNode::neighbor(ConnectionDirection cdir) const
{
	switch (cdir) {
	case DIR_X: return toRank(x + 1, y, z);
	case DIR_MINUS_X: return toRank(x - 1, y, z);
	case DIR_Y: return toRank(x, y + 1, z);
	case DIR_MINUS_Y: return toRank(x, y - 1, z);
    case DIR_Y_PERIOD_FIRST: return toRank(x, y + 1 - gridDimensions.y, z);
    case DIR_Y_PERIOD_LAST: return toRank(x, y + gridDimensions.y - 1, z);
	case DIR_Z: return toRank(x, y, z + 1);
	case DIR_MINUS_Z: return toRank(x, y, z - 1);
	default:
		MY_ASSERT(false);
		return -1;
	}
}

bool ComputeNode::hasNeighbor(ConnectionDirection cdir) const
{
	return neighbor(cdir) != -1;
}


void ComputeNode::fillGridDimensions()
{
    gridDimensions = getGridDimensions()[mpi.procCount];
}

void ComputeNode::fillXYZ()
{
    z = mpi.rank % gridDimensions.z;
    uint yxRank = mpi.rank / gridDimensions.z;
    y = yxRank % gridDimensions.y;
    x = yxRank / gridDimensions.y;
}

int ComputeNode::toRank(int i, int j, int k) const
{
    if (i < 0 || j < 0 || k < 0
            || i >= gridDimensions.x
            || j >= gridDimensions.y
            || k >= gridDimensions.z) {
		return -1;
	}
    return (i * gridDimensions.y + j) * gridDimensions.z + k;

}

void ComputeNode::print(const std::string &str) const
{
    std::cout << titledStr(str) << std::endl;
}

void ComputeNode::error(const std::string &err) const
{
    std::cerr << titledStr(err) << std::endl;
}

std::string ComputeNode::titledStr(const std::string &str) const
{
    return SSTR(mpi.rank << ": " << str);
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

    void step()
    {
        MY_ASSERT(_started);

        _wall_clock_elapsed = MPI_Wtime() - _wstart;
    }

    void finish()
    {
        step();

        clear();
    }

    void print() const
    {
        char buff[256];
        snprintf(buff, sizeof(buff),
                 "WC: %lf s.\n",
                 _wall_clock_elapsed);

        cnode.print(buff);
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

Profiler::~Profiler()
{
    delete _impl;
}

void Profiler::start()
{
    _impl->start();
}

void Profiler::step()
{
    _impl->step();
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
    const char s_Deviation [] = "deviation";

    std::string sarg(arg);
    if (sarg.find(s_Deviation) != std::string::npos) {
        deviation = true;
        return;
    }
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

    cnode.print(std::string("unrecognized argument ") + sarg);
}

