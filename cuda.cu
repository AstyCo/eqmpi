#include "cuda.hpp"

#include <math.h>

#define LOAD_CONST_MEM(dest, src)\
    cudaMemcpyToSymbol(dest, &src, sizeof(dest), 0, cudaMemcpyHostToDevice)

__constant__ static int N; /// size of 3D grid

__constant__ static int i0;  /// first index of 3D grid
__constant__ static int j0;  /// first index of 3D grid
__constant__ static int k0;  /// first index of 3D grid
__constant__ static long long bigsize; /// size of extended grid

__constant__ static long ic;  /// count in each dimension of 3D grid
__constant__ static long jc;  /// count in each dimension of 3D grid
__constant__ static long kc;  /// count in each dimension of 3D grid

__constant__ static int ei;  /// first edge id in each dimension
__constant__ static int ej;  /// first edge id in each dimension
__constant__ static int ek;  /// first edge id in each dimension

__constant__ static int eic; /// size of edge in each dimension
__constant__ static int ejc; /// size of edge in each dimension
__constant__ static int ekc; /// size of edge in each dimension

__constant__ static real hi; /// step in each dimension
__constant__ static real hj; /// step in each dimension
__constant__ static real hk; /// step in each dimension

__constant__ static real ht; /// step time

__constant__ static RealDVector *array;
__constant__ static RealDVector *arrayP;
__constant__ static RealDVector *arrayPP;

__constant__ static RealDVector *edgeArray;

__constant__ static LongDVector *edgeIndices;

__device__
long get_index(uint i, uint j, uint k) const
{
    return (long(i + 1) * (jc + 2) + (j + 1)) * (kc + 2) + (k + 1);
}

struct ExactIndex
{
    int i, j, k;
    ExactIndex(long id)
    {
        k = id % kc;
        long ij = id / kc;
        j = ij % jc;
        i = ij / jc;
    }
};

__device__
real x(int i) const { return (i0 + i) * hi;}

__device__
real y(int j) const { return (j0 + j) * hj;}

__device__
real z(int k) const { return (k0 + k) * hk;}

struct ExactIndex
{
    int i, j, k;
    ExactIndex(long id)
    {
        k = id % kc;
        long ij = id / kc;
        j = ij % jc;
        i = ij / jc;
    }
};

struct Index
{
    int i, j, k;

    Index(long id)
    {
        k = id % (kc + 2);
        long ij = id / (kc + 2);
        j = ij % (jc + 2);
        i = ij / (jc + 2);
    }
};

__device__
void calculate(long offset)
{
    int i = id.i, j = id.j, k = id.k;
    array[offset] = 2 * arrayP[offset] - arrayPP[offset]
            + ht * ht * (
                (arrayP[get_index(i-1,j,k)]
                - 2 * arrayP[offset]
                + arrayP[get_index(i+1,j,k)]) / hi / hi
            + (arrayP[get_index(i,j-1,k)]
            - 2 * arrayP[offset]
            + arrayP[get_index(i,j+1,k)]) / hj / hj
            + (arrayP[get_index(i,j,k-1)]
            - 2 * arrayP[offset]
            + arrayP[get_index(i,j,k+1)]) / hk / hk
            );
}


struct FStep0
{
    __device__
    real operator(int offset) {
        Index id(offset);
        return phi(x(id.i), y(id.j), z(id.k));
    }
};

void cuda_step_0()
{
    counting_iterator<int> first(0);
    thrust::transform(first, first + bigsize,
                      arrayPP.begin(), FStep0()); // install PHI

    thrust::copy(arrayPP.begin(), arrayPP.end(),
                 array.begin()); // install 0-type boundary conditions
}

struct FStep1
{
    __device__
    real operator(int offset) {
        Index id(offset);
        return arrayPP[offset]
                + ht * ht / 2 * div_grad_phi(x(id.i), y(id.j), z(id.k));
    }
};

void cuda_step_1()
{
    counting_iterator<int> first(0);
    thrust::transform(first, first + bigsize,
                      arrayP.begin(), FStep1()); // install LAPLACIAN PHI
}

struct FCalculateInner
{
    __device__
    void operator(int offset) {
        Index id(offset);
        if (id.i < 2 || id.j < 2 || id.k < 2
                || id.i > ic - 1
                || id.j > jc - 1
                || id.k > kc - 1) {
            // don't change
        }
        else {
            int i = id.i, j = id.j, k = id.k;
            array[offset] = 2 * arrayP[offset] - arrayPP[offset]
                    + ht * ht * (
                        (arrayP[get_index(i-1,j,k)]
                        - 2 * arrayP[offset]
                        + arrayP[get_index(i+1,j,k)]) / hi / hi
                    + (arrayP[get_index(i,j-1,k)]
                    - 2 * arrayP[offset]
                    + arrayP[get_index(i,j+1,k)]) / hj / hj
                    + (arrayP[get_index(i,j,k-1)]
                    - 2 * arrayP[offset]
                    + arrayP[get_index(i,j,k+1)]) / hk / hk
                    );
        }
    }
};

void cuda_calculate_inner()
{
    counting_iterator<int> first(0);
    thrust::for_each(first, first + bigsize,
                     FCalculateInner()); // calculate inner val
}

struct FCalculateEdge
{
    __device__
    real operator(int offset) {
        return calculate(offset);
    }
};

void cuda_calculate_edges()
{
    thrust::transform(edgeIndices.begin(), edgeIndices.end(),
                      edgeArray.begin(), FCalculateEdge());
}

void cuda_load_const_mem(int N_,
                         int i0_,int j0_, int k0_,
                         long long bigsize_,
                         long ic_, long jc_, long kc_,
                         int ei_, int ej_, int ek_,
                         int eic_, int ejc_, int ekc_,
                         real hi_, real hj_, real hk_,
                         real ht_,
                         RealDVector *array_,
                         RealDVector *arrayP_,
                         RealDVector *arrayPP_,
                         RealDVector *edgeArray_,
                         LongDVector *edgeIndices_)
{
    LOAD_CONST_MEM(N, N_);

    LOAD_CONST_MEM(i0, i0_);
    LOAD_CONST_MEM(j0, j0_);
    LOAD_CONST_MEM(k0, k0_);

    LOAD_CONST_MEM(bigsize, bigsize_);

    LOAD_CONST_MEM(ic, ic_);
    LOAD_CONST_MEM(jc, jc_);
    LOAD_CONST_MEM(kc, kc_);

    LOAD_CONST_MEM(ei, ei_);
    LOAD_CONST_MEM(ej, ej_);
    LOAD_CONST_MEM(ek, ek_);

    LOAD_CONST_MEM(eic, eic_);
    LOAD_CONST_MEM(ejc, ejc_);
    LOAD_CONST_MEM(ekc, ekc_);

    LOAD_CONST_MEM(hi, hi_);
    LOAD_CONST_MEM(hj, hj_);
    LOAD_CONST_MEM(hk, hk_);

    LOAD_CONST_MEM(ht, ht_);

    LOAD_CONST_MEM(array, array_);
    LOAD_CONST_MEM(arrayP, arrayP_);
    LOAD_CONST_MEM(arrayPP, arrayPP_);

    LOAD_CONST_MEM(edgeArray, edgeArray_);

    LOAD_CONST_MEM(edgeIndices, edgeIndices_);
}


