#include <unordered_map>
#include <vector>

#include <cstring>
#include <cmath>

#define SH(x) sinh(x) // different floating point precision *hf, *h, *hl
#define CH(x) cosh(x) // different floating point precision *hf, *h, *hl

typedef unsigned int uint;

#ifdef FLOAT_P
typedef float real;
#   ifndef MPI_TYPE_REAL
#       define MPI_TYPE_REAL MPI_FLOAT
#   endif
#else
typedef double real;
#   ifndef MPI_TYPE_REAL
#       define MPI_TYPE_REAL MPI_FLOAT
#   endif
#endif;

enum MPI_SENDRECV_TAGS
{
    TAG_BOUNDARY_ELEMENTS
};

real phi(real x, real y, real /*z*/)
{
    return SH(x) + CH(y);
}

// Δ = div grad
// Δ phi = (d/dx^2 + d/dy^2 + d/dz^2) phi
real div_grad_phi(real x, real y, real /*z*/)
{
    return SH(x) + CH(y); // sh(x)' = ch(x), ch(x)' = sh(x)
}

struct GridDimensions
{
    uint x;
    uint y;
    uint z;

    GridDimensions(uint x_ = 0, uint y_ = 0, uint z_ = 0)
        : x(x_), y(y_), z(z_)
    {}
};

typedef std::unordered_map<uint, GridDimensions> MapGridDimensions;

MapGridDimensions getGridDimensions()
{
    MapGridDimensions mts;

    mts.insert(128, GridDimensions(4, 4, 8));
    mts.insert(256, GridDimensions(8, 4, 8));
    mts.insert(512, GridDimensions(8, 8, 8));

    return mts;
}

struct ComputeNode
{
    // MPI
    struct MPI_data
    {
        uint rank;
//        uint commSize;
        uint procCount;
    } mpi;


    GridDimensions gridDimensions;
    // Equation
    int x;
    int y;
    int z;


    void fillGridDimensions()
    {
        gridDimensions = getGridDimensions()[mpi.procCount];
    }

    void fillXYZ()
    {
        x = mpi.rank % gridDimensions.x;
        uint yzRank = mpi.rank / gridDimensions.x;
        y = yzRank % gridDimensions.y;
        z = yzRank / gridDimensions.y;
    }
};

enum ConnectionDirection
{
    DIR_X,
    DIR_MINUS_X,
    DIR_Y,
    DIR_MINUS_Y,
    DIR_Z,
    DIR_MINUS_Z,

    DIR_SIZE
};

ConnectionDirection toCD(int i)
{
    return static_cast<ConnectionDirection>(i);
}

struct Iterations
{
    typedef std::vector<real> RealVector;
    typedef void (Iterations::*CopyMFuncPtr)(RealVector &, uint, uint, uint);
    typedef void (Iterations::*IndexesMFuncPtr)(uint, uint, uint);


    static uint K; // time step count
    static real T;
    static real Lx;
    static real Ly;
    static real Lz;

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

    uint dir_size(ConnectionDirection cdir)
    {
        switch (cdir) {
        case DIR_X:
        case DIR_MINUS_X:
            return jc * kc;
        case DIR_Y:
        case DIR_MINUS_Y:
            return jc * kc;
        case DIR_Z:
        case DIR_MINUS_Z:
            return jc * kc;
        default:
            MY_ASSERT(false);
            return 0;
        }
    }

    struct Requests
    {
        struct Info
        {
            ConnectionDirection dir;
        };

        std::vector<MPI_Request> v;
        std::vector<RealVector> buff;
        std::vector<Info> iv;

        void append(ConnectionDirection cdir, uint sz)
        {
            Info i;
            i.dir = cdir;

            iv.push_back(i);

            v.push_back(MPI_Request);

            buff.push_back(RealVector());
            buff.back().reserve(sz);
        }
    };

    Requests recv_requests;
    Requests send_requests;


    Iterations(const ComputeNode &n)
    {
        fill(n);
    }

    void fill(const ComputeNode &n)
    {
        MY_ASSERT(0 == n.mpi.procCount % n.gridDimensions.x);
        MY_ASSERT(0 == n.mpi.procCount % n.gridDimensions.y);
        MY_ASSERT(0 == n.mpi.procCount % n.gridDimensions.z);

        ic = n.mpi.procCount / n.gridDimensions.x;
        jc = n.mpi.procCount / n.gridDimensions.y;
        kc = n.mpi.procCount / n.gridDimensions.z;

        i0 = n.x * ic;
        j0 = n.y * jc;
        k0 = n.z * kc;

        imax = i0 + ic;
        jmax = j0 + jc;
        kmax = k0 + kc;

        hx = Lx / ic;
        hy = Ly / jc;
        hz = Lz / kc;

        ht = T / K;

        size = ic * jc * kc;
        bigsize = (ic + 2) * (jc + 2) * (kc + 2);

        // optimization (allocations may throw std::bad_alloc if no enough memory)
        array.reserve(bigsize);
        arrayP.reserve(bigsize);
        arrayPP.reserve(bigsize);

        for (int i = 0; i < DIR_SIZE; ++i) {
            ConnectionDirection cdir = toCD(i);
            if (n.is(cdir))
                requests.append(cdir, dir_size(cdir));
        }
    }

    void copy_recv(RealVector &v, uint i, uint j, uint k)
    {
        arrayP[get_p_index(i, j, k)] = v[j * kc + k];
    }

    void copy_send(RealVector &v, uint i, uint j, uint k)
    {
        v[j * kc + k] = arrayP[get_p_index(i, j, k)];
    }

    void copy_data(uint id, CopyMFuncPtr f)
    {
        ConnectionDirection = requests.iv[id].dir;
        RealVector &v = requests.v[id];
        switch (cdir) {
        case DIR_X:
        case DIR_MINUS_X:
            uint i = ((cdir == DIR_X)
                        ? ic - 1
                        : 0);
            for (uint j = 0; j < jc; ++j) {
                for (uint k = 0; k < kc; ++k) {
                    this->*f(v, i, j, k);
                }

            }
            break;
        }
    }

    void calculate(ConnectionDirection cdir)
    {
        switch (cdir) {
        case DIR_X:
        case DIR_MINUS_X:
            uint i = ((cdir == DIR_X) ? ic - 1
                                      : 0);
            for (uint j = 1; j < jc - 1; ++j) {
                for (uint k = 1; k < kc - 1; ++k)
                    calculate(i, j, k);
            }
            break;
        }
    }

    void calculate(uint i, uint j, uint k)
    {
        uint index = get_index(i, j, k);
        uint p_index = get_p_index(i, j, k);
        uint pp_index = index;
        array[index] = 2 * arrayP[p_index] - arrayPP[pp_index]
                + ht * ht * (
                    (array[get_index(i-1,j,k)]
                    - 2 * arrayP[p_index]
                    + array[get_index(i+1,j,k)]) / hx / hx
                + (array[get_index(i,j-1,k)]
                - 2 * arrayP[p_index]
                + array[get_index(i,j+1,k)]) / hy / hy
                + (array[get_index(i,j,k-1)]
                - 2 * arrayP[p_index]
                + array[get_index(i,j,k+1)]) / hz / hz
                );
    }
   void it_for_each(IndexesMFuncPtr func)
   {
       for (uint i = 0; i < ic; ++i) {
           for (uint j = 0; j < jc; ++j) {
               for (uint k = 0; k < kc; ++k) {
                   this->*func(i, j, k);
               }
           }
       }
   }

    void next_step()
    {
        // array -> arrayP, arrayP -> arrayPP
        memcpy(arrayPP.data(), arrayP.data(), arrayP.size() * sizeof(real));
        memcpy(arrayPP.data(), arrayP.data(), arrayP.size() * sizeof(real));
    }

    uint get_index(uint i, uint j, uint k) const
    {
        return get_p_index(i,j,k);
//        return (i * jc + j) * kc + k;
    }

    uint get_p_index(uint i, uint j, uint k) const
    {
        return ((i + 1) * (jc + 2) + (j + 1)) * (kc + 2) + (k + 1);
    }

    real x(uint i) const { return (i0 + i) * hx;}
    real y(uint j) const { return (j0 + j) * hy;}
    real z(uint k) const { return (k0 + k) * hz;}

    void set_0th(uint i, uint j, uint k)
    {
        uint index = get_index(i, j, k);
        arrayPP[index] = phi(x(i), y(j), z(k));
    }


    void step0()
    {
        it_for_each(set_0th);
    }

    void set_1th(uint i, uint j, uint k)
    {
        uint index = get_index(i, j, k);
        arrayP[index] = arrayPP[index] + ht * ht / 2 * div_grad_phi(x(i), y(j), z(k));
    }

    void step1()
    {
        it_for_each(set_1th)
    }

};


int main(int argc, char *argv[])
{
    ComputeNode cnode;
    // fill cnode
    // TODO
    Iterations its(cnode);

    its.step0();
    its.step1();

    // STEPS
    for (uint step = 0; step < its.K; ++step) {
        // async receive prev
#       pragma omp parallel for // parallel slices
        for (uint i = 1; i < its.ic - 1; ++i) {
            for (uint j = 1; j < its.jc - 1; ++j) {
                for (uint k = 1; k < its.kc -1; ++k) {
                    its.calculate(i, j, k);
                }
            }
        }
        // sync recv prev
        int request_index;
        MPI_Status status;
#       pragma omp parallel // parallel recv
        {
            const Iterations::Requests &requests = its.recv_requests;
            while (MPI_UNDEFINED != MPI_Waitany(requests.v.data(),
                                                requests.v.size(),
                                                &request_index, &status)) {
                const Iterations::Requests::Info &info
                        = requests.iv[request_index];
                its.recv_data(request_index);
                its.calculate(info.dir);
            }
        }
        its.calculate_angle_values();
        its.next_step(); // update arrays
        // async send prev
        {
            const Iterations::Requests &requests = its.send_requests;
            MPI_Waitall(requests.v.data(),
                        requests.v.size(),
                        &status); // wait for every asynchronous send (so we can use send buffers again)
            // asynchronous send to every neighbor processor
            for (uint i = 0; i < requests.size(); ++i) {
                its.copy_data(i, &Iterations::Requests::copy_send)
                MPI_Isend(requests.buff[i].data(),
                          requests.buff[i].size(),
                          MPI_TYPE_REAL,
                          cnode.neighbor(requests.iv[i].dir),
                          TAG_BOUNDARY_ELEMENTS,
                          MPI_COMM_WORLD,
                          &requests.v[i]);
            }
        }
    } // ENDS STEPS
    return 0;
}
