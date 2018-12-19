#include "utils.hpp"

#include <iostream>
#include <limits>

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>

#ifdef CUDA
#include "cuda_runtime.h"
#endif

DetailedTimes times;
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
    sc = SCBluegeneP;

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

void ComputeNode::print0(const std::string &str) const
{
    if (cnode.mpi.rank == 0)
        print(str);
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

#ifdef CUDA
        cudaEventCreate(&_start);
        cudaEventCreate(&_stop);

        cudaEventRecord(_start, 0);
#else
        _wstart = MPI_Wtime();
#endif
    }

    void clear()
    {
        _started = false;
    }

    void step()
    {
        MY_ASSERT(_started);

#ifdef CUDA
        cudaEventRecord(_stop, 0);
        cudaEventSynchronize (_stop);

        float elapsed;
        cudaEventElapsedTime(&elapsed, _start, _stop);
        _wall_clock_elapsed = elapsed / 1000;
        cudaEventDestroy(_start);
        cudaEventDestroy(_stop);
#else
        _wall_clock_elapsed = MPI_Wtime() - _wstart;
#endif
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

    double _wall_clock_elapsed;

#ifdef CUDA
    cudaEvent_t _start, _stop;
#else
    double _wstart;
#endif

    Profiler::Options _opts;
};

Profiler::Profiler(Profiler::Options opts)
    : _impl(new ProfilerPrivate(opts))
{
    if (opts & StartOnConstructor)
        start();
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

double Profiler::time()
{
    step();
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

DetailedTimes::DetailedTimes()
{
    clear();
}

void DetailedTimes::clear()
{
    program_initialization = 0;
    allocations = 0;
    program_finalization = 0;
    parallel_cycles = 0;
    host_device_exchange = 0;
    mpi_send_recv = 0;
    shift_arrays = 0;
    total = 0;
}

std::string DetailedTimes::get_times(long N)
{
    return SSTR("###"
                << ',' << cnode.mpi.procCount
                << ',' << N
                << ',' << host_device_exchange
                << ',' << mpi_send_recv
                << ',' << parallel_cycles
                << ',' << program_initialization
                << ',' << program_finalization
                << ',' << shift_arrays
                << ',' << sum_time()
                << ',' << total);
}

void DetailedTimes::reduce()
{
    reduce(program_initialization);
    reduce(allocations);
    reduce(program_finalization);
    reduce(parallel_cycles);
    reduce(host_device_exchange);
    reduce(mpi_send_recv);
    reduce(shift_arrays);
    reduce(total);
}

void DetailedTimes::reduce(double &var)
{
    double global = 0;
    MPI_Reduce(&var, &global, 1, MPI_DOUBLE,
               MPI_SUM, 0, MPI_COMM_WORLD);
    var = global / cnode.mpi.procCount;
}

//void DetailedTimes::rel()
//{
//    sum = sum_time();
//    double coef = total / sum_time();

//    program_initialization *= coef;
//    allocations *= coef;
//    program_finalization *= coef;
//    parallel_cycles *= coef;
//    host_device_exchange *= coef;
//    mpi_send_recv *= coef;
//    shift_arrays *= coef;
//}

double DetailedTimes::sum_time() const
{
    return program_initialization
            + program_finalization
            + parallel_cycles
            + host_device_exchange
            + mpi_send_recv
            + shift_arrays;
}

void get_time(double &dest, double &local)
{
    dest += local;
}

void get_time(double &dest, Profiler &p)
{
    double time = p.time();
    get_time(dest, time);

    p.start();
}
