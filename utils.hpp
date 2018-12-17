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
        StartOnConstructor = 0x1,
        PrintOnDestructor = 0x2
    };

    explicit Profiler(Options opts = StartOnConstructor);
    ~Profiler();

    void start();
    void step();
    void finish();

    double time();
    void print() const;

private:
    ProfilerPrivate *_impl;
};

struct CommandLineArgs
{
    int argc;
    char **argv;

    int K;
    int N;
    bool deviation;

    CommandLineArgs()
    {
        // default
        K = 20;
        N = -1;
        deviation = false;
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

    enum Supercomputer
    {
        SCPolus,
        SCBluegeneP
    } sc;
    std::string scTag() const;


    GridDimensions gridDimensions;
    // Equation
    int x;
    int y;
    int z;

    void init();

    int neighbor(ConnectionDirection cdir) const;
    bool hasNeighbor(ConnectionDirection cdir) const;
    void fillGridDimensions();
    void fillXYZ();

    int toRank(int i, int j, int k) const;

    void print(const std::string &str) const;
    void print0(const std::string &str) const;
    void error(const std::string &err) const;
private:
    std::string titledStr(const std::string &str) const;
};

extern ComputeNode cnode;

struct DetailedTimes
{
    double program_initialization;
    double program_finalization;
    double parallel_cycles;
    double host_device_exchange;
    double mpi_send_recv;

    double shift_arrays;

    DetailedTimes();
    void clear();
    std::string get_times();
};

void get_time(double &dest, double &local);
void get_time(double &dest, Profiler &p);

extern DetailedTimes times;


#endif // UTILS_HPP
