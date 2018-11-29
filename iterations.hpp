#ifndef ITERATIONS_HPP
#define ITERATIONS_HPP

#include "utils.hpp"

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

        void append(ConnectionDirection cdir, uint sz);
    };


    static uint K = 20; // time step count
    static real T = 10;
    static real Lx = 10;
    static real Ly = 20;
    static real Lz = 30;

    const uint i0;
    const uint j0;
    const uint k0;

    const uint ic;  // counts
    const uint jc;  // counts
    const uint kc;  // counts
    const uint size;
    const uint bigsize;

    const uint imax;
    const uint jmax;
    const uint kmax;

    const real hx;
    const real hy;
    const real hz;

    const uint ht; // delta t

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


    Iterations(const ComputeNode &n);

    void prepare();
    void run();


    uint dir_size(ConnectionDirection cdir);

    void fill(const ComputeNode &n);

    void copy_recv(RealVector &v, uint i, uint j, uint k);
    void copy_send(RealVector &v, uint i, uint j, uint k);

    void copy_data(uint id, CopyMFuncPtr f);
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