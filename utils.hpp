#ifndef UTILS_HPP
#define UTILS_HPP

#include "globals.hpp"

void initialize_seed();
double randomize(double min, double max);

class ProfilerPrivate;
class Profiler
{
public:
    enum Options
    {
        Default = 0x0,
        PrintOnDestructor = 0x1
    };

    explicit Profiler(Options opts = Default);

    void start();
    void finish();

    double time() const;
    void print() const;

private:
    ProfilerPrivate *_impl;
};

struct CommandLineArgs
{
    int argc;
    char **argv;

    uint K;
    uint N;

    CommandLineArgs()
    {
        // default
        K = 20;
        N = 50;
    }

    void parse(int argc_, char **argv_);
    void parseArg(char arg[]);
};

extern CommandLineArgs clargs;

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
