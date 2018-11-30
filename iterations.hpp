#ifndef ITERATIONS_HPP
#define ITERATIONS_HPP

#include "utils.hpp"

#include <vector>

struct Iterations
{
    typedef std::vector<real> RealVector;
    typedef void (Iterations::*CopyMFuncPtr)(RealVector &, uint, uint, uint);
    typedef void (Iterations::*IndexesMFuncPtr)(uint, uint, uint);

    struct Requests
    {
        struct Info
        {
            ConnectionDirection dir;
        };

        std::vector<MPI_Request> v;
        std::vector<RealVector> buff;
        std::vector<Info> iv;

        uint size() const { return v.size();}
        void append(ConnectionDirection cdir, uint sz);
    };


    static uint K; // time step count
    static real T;
    static real Lx;
    static real Ly;
    static real Lz;

    uint i0;
    uint j0;
    uint k0;

    uint ic;  // counts
    uint jc;  // counts
    uint kc;  // counts
    uint size;
    uint bigsize;

    uint imax;
    uint jmax;
    uint kmax;

    real hx;
    real hy;
    real hz;

    uint ht; // delta t

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

    Requests recv_requests;
    Requests send_requests;

    int step;

    const ComputeNode &cnode;

    Iterations(const ComputeNode &n);

    void prepare();
    void run();


    uint dir_size(ConnectionDirection cdir);

    void fill(const ComputeNode &n);

    void copy_recv(RealVector &v, uint i, uint j, uint k);
    void copy_send(RealVector &v, uint i, uint j, uint k);

    void copy_data(Requests &requests, uint id, CopyMFuncPtr f);
    void calculate(uint i, uint j, uint k);
    void calculate(ConnectionDirection cdir);
    void calculate_edge_values();
   	void it_for_each(IndexesMFuncPtr func);
    void next_step();

    uint get_index(uint i, uint j, uint k) const;

    uint get_p_index(uint i, uint j, uint k) const;

    real x(uint i) const { return (i0 + i) * hx;}
    real y(uint j) const { return (j0 + j) * hy;}
    real z(uint k) const { return (k0 + k) * hz;}

    void set_0th(uint i, uint j, uint k);
    void step0();
    void set_1th(uint i, uint j, uint k);
    void step1();

};

#endif // ITERATIONS_HPP
