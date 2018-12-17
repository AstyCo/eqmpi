#ifndef ITERATIONS_HPP
#define ITERATIONS_HPP

#include "utils.hpp"
#include "cuda.hpp"

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

    long N;

    uint i0;
    uint j0;
    uint k0;

    long ic;  // counts
    long jc;  // counts
    long kc;  // counts

    long bigsize;
    long totalEdgeSize;

    real hx;
    real hy;
    real hz;

    real ht; // delta t

    RealDVector dArray;
    RealDVector dArrayP;
    RealDVector dArrayPP;

    RealDVector dEdgeArray;
    RealHVector hEdgeArray;

    LongHVector hEdgeIndices;
    LongDVector dEdgeIndices;

    RealDVector dDeviationsArray;

    RealVector hArrayBuff;

    uint edgeI, edgeIL, edgeJ, edgeJL, edgeK, edgeKL;

    Requests recv_requests;
    Requests send_requests;

    int next_step;

    Iterations(uint N_);

    void prepare();
    void prepareEdgeIndiceArray();
    void prepareEdgeIndices();

    void run();

    void copy_edges_to_h();
    void copy_edges_to_d();

    void async_send_all();
    void async_recv_all();

    void printDeviations(uint n);

    uint dir_size(ConnectionDirection cdir);

    void fill(const ComputeNode &n);

    void copy_recv(RealVector &v, RealVector &a, int i, int j, int k, uint offset);
    void copy_send(RealVector &v, RealVector &a, int i, int j, int k, uint offset);

    void copy_data(Requests &requests, uint id, MPI_OP type);

    long get_index(uint i, uint j, uint k) const;

    real x(uint i) const { return (i0 + i) * hx;}
    real y(uint j) const { return (j0 + j) * hy;}
    real z(uint k) const { return (k0 + k) * hz;}
    real time(uint n) const { return n * ht;} // n - iter

    void step0();
    void step1();

    int edgeId(ConnectionDirection cdir, MPI_OP op_type);
    int sendEdgeId(ConnectionDirection cdir) const;
    int recvEdgeId(ConnectionDirection cdir) const;
    int sendrecvEdgeId(ConnectionDirection cdir) const;
};

#endif // ITERATIONS_HPP
