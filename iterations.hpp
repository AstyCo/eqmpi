#ifndef ITERATIONS_HPP
#define ITERATIONS_HPP

#include "utils.hpp"

#include <vector>

struct Iterations
{
    typedef std::vector<real> RealVector;
    typedef void (Iterations::*CopyMFuncPtr)(RealVector &, RealVector &,
                                             int, int, int, uint);
    typedef void (Iterations::*IndexesMFuncPtr)(uint, uint, uint);

    enum MPI_OP
    {
        MPI_OP_RECV,
        MPI_OP_SEND
    };

    struct Requests
    {
        struct Info
        {
            ConnectionDirection dir;
        };

        std::vector<MPI_Request> v;
        std::vector<RealVector> buffs;
        std::vector<Info> iv;
        std::vector<MPI_Status> statuses;

        uint size() const { return v.size();}
        void append(Info info, uint sz);

        Requests()
        {
            buffs.reserve(DIR_SIZE);
        }
    };

    uint N;

    uint i0;
    uint j0;
    uint k0;

    uint ic;  // counts
    uint jc;  // counts
    uint kc;  // counts
    uint bigsize;

    real hx;
    real hy;
    real hz;

    real ht; // delta t

    RealVector array;
    RealVector arrayP;
    RealVector arrayPP;

    RealVector sendX;
    RealVector sendXm;
    RealVector sendY;
    RealVector sendYm;
    RealVector sendZ;
    RealVector sendZm;

    RealVector recvX;
    RealVector recvXm;
    RealVector recvY;
    RealVector recvYm;
    RealVector recvZ;
    RealVector recvZm;

    RealVector analyticalSolution;

    std::vector<ConnectionDirection> no_neighbour_edges;
    Requests recv_requests;
    Requests send_requests;

    int next_step;

    Iterations(uint N_);

    void prepare();
    void run();

    void prepareSolution(uint n);
    void printDeviation(uint i, uint j, uint k, uint n);

    void printDeviationPrivate(const RealVector &arr, uint i, uint j, uint k, uint n);
    real getDeviation(const RealVector &arr, uint i, uint j, uint k, uint n) const;

    void printDeviations(uint n);
    void printDeviationsPrivate(const RealVector &arr, uint n);

    uint dir_size(ConnectionDirection cdir);

    void fill(const ComputeNode &n);

    void copy_recv(RealVector &v, RealVector &a, int i, int j, int k, uint offset);
    void copy_send(RealVector &v, RealVector &a, int i, int j, int k, uint offset);

    void copy_data(Requests &requests, uint id, MPI_OP type);
    void calculate(uint i, uint j, uint k);
    void calculate(ConnectionDirection cdir);
    void calculate_edge_values();
   	void it_for_each(IndexesMFuncPtr func);
    void shift_arrays();

    uint get_index(uint i, uint j, uint k) const;

    real x(uint i) const { return (i0 + i) * hx;}
    real y(uint j) const { return (j0 + j) * hy;}
    real z(uint k) const { return (k0 + k) * hz;}
    real time(uint n) const { return n * ht;} // n - iter

    void set_0th(uint i, uint j, uint k);
    void step0();
    void set_1th(uint i, uint j, uint k);
    void step1();

    int edgeId(ConnectionDirection cdir, MPI_OP op_type);
    int sendEdgeId(ConnectionDirection cdir) const;
    int recvEdgeId(ConnectionDirection cdir) const;
    int sendrecvEdgeId(ConnectionDirection cdir) const;
};

#endif // ITERATIONS_HPP
