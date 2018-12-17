#include "cuda.hpp"
#include "utils.hpp"

#include <thrust/iterator/counting_iterator.h>

#include <math.h>

#define LOAD_CONST_MEM(dest, src)\
    cudaMemcpyToSymbol(dest, &src, sizeof(dest), 0, cudaMemcpyHostToDevice)

__constant__ static int N; /// size of 3D grid

__constant__ static int i_first;  /// first index of 3D grid
__constant__ static int j_first;  /// first index of 3D grid
__constant__ static int k_first;  /// first index of 3D grid
__constant__ static long long bigsize; /// size of extended grid

__constant__ static long i_count;  /// count in each dimension of 3D grid
__constant__ static long j_count;  /// count in each dimension of 3D grid
__constant__ static long k_count;  /// count in each dimension of 3D grid

__constant__ static int edge_i_first;  /// first edge id in each dimension
__constant__ static int edge_j_first;  /// first edge id in each dimension
__constant__ static int edge_k_first;  /// first edge id in each dimension

__constant__ static int edge_i_count; /// size of edge in each dimension
__constant__ static int edge_j_count; /// size of edge in each dimension
__constant__ static int edge_k_count; /// size of edge in each dimension

__constant__ static real step_i; /// step in each dimension
__constant__ static real step_j; /// step in each dimension
__constant__ static real step_k; /// step in each dimension

__constant__ static real step_t; /// step time

__constant__ static RealDVector *array;
__constant__ static RealDVector *arrayP;
__constant__ static RealDVector *arrayPP;

__constant__ static real *arr;
__constant__ static real *arrP;
__constant__ static real *arrPP;



__constant__ static RealDVector *edgeArray;

__constant__ static LongDVector *edgeIndices;


__device__
static real phi(real x, real y, real z)
{
    return sin(x) * cos(y - VAL_LY/2) * sin(z);
}

// Δ = div grad
// Δ phi = (d/dx^2 + d/dy^2 + d/dz^2) phi
__device__
static real div_grad_phi(real x, real y, real z)
{
    return -3.0 * phi(x, y, z);
}

__device__
static long get_index(uint i, uint j, uint k)
{
    return (long(i + 1) * (j_count + 2) + (j + 1)) * (k_count + 2) + (k + 1);
}

__device__
static real x_val(int i) { return (i_first + i) * step_i;}

__device__
static real y_val(int j) { return (j_first + j) * step_j;}

__device__
static real z_val(int k) { return (k_first + k) * step_k;}

struct Index
{
    int i, j, k;

    __device__
    Index(long id)
    {
        k = id % (k_count + 2);
        long ij = id / (k_count + 2);
        j = ij % (j_count + 2);
        i = ij / (j_count + 2);
    }
};

__device__
void calculate(long offset)
{
    Index id(offset);
    int i = id.i, j = id.j, k = id.k;

    arr[offset] = 2 * arrP[offset] - arrPP[offset]
            + step_t * step_t * (
                (arrP[get_index(i-1,j,k)]
                - 2 * arrP[offset]
                + arrP[get_index(i+1,j,k)]) / step_i / step_i
            + (arrP[get_index(i,j-1,k)]
            - 2 * arrP[offset]
            + arrP[get_index(i,j+1,k)]) / step_j / step_j
            + (arrP[get_index(i,j,k-1)]
            - 2 * arrP[offset]
            + arrP[get_index(i,j,k+1)]) / step_k / step_k
            );
}




struct FStep0
{
    __device__
    real operator()(int offset) {
        Index id(offset);
        return phi(x_val(id.i), y_val(id.j), z_val(id.k));
    }
};

void cuda_resize(RealDVector &dArray,
                 RealDVector &dArrayP,
                 RealDVector &dArrayPP,
                 RealDVector &dEdgeArray,
                 RealHVector &hEdgeArray,
                 RealDVector &analyticalSolution,
                 long totalEdgeSize,
                 long bigsize)
{
    dArray.resize(bigsize);
    dArrayP.resize(bigsize);
    dArrayPP.resize(bigsize);

    dEdgeArray.resize(totalEdgeSize);
    hEdgeArray.resize(totalEdgeSize);

    if (clargs.deviation)
        analyticalSolution.resize(bigsize);
}

void cuda_step_0(RealDVector &dArray, RealDVector &dArrayPP)
{
    thrust::counting_iterator<int> first(0);
    thrust::transform(first, first + dArrayPP.size(),
                      dArrayPP.begin(), FStep0()); // install PHI

    thrust::copy(dArrayPP.begin(), dArrayPP.end(),
                 dArray.begin()); // install 0-type boundary conditions
}

struct FStep1
{
    __device__
    real operator()(int offset) {
        Index id(offset);
        return arrP[offset]
                + step_t * step_t / 2 * div_grad_phi(x_val(id.i), y_val(id.j), z_val(id.k));
    }
};

void cuda_step_1(RealDVector &dArrayP)
{
    thrust::counting_iterator<int> first(0);
    thrust::transform(first, first + dArrayP.size(),
                      dArrayP.begin(), FStep1()); // install LAPLACIAN PHI
}

struct FCalculateInner
{
    __device__
    void operator()(int offset) {
        Index id(offset);
        if (id.i < 2 || id.j < 2 || id.k < 2
                || id.i > i_count - 1
                || id.j > j_count - 1
                || id.k > k_count - 1) {
            // don't change
        }
        else {
            calculate(offset);
        }
    }
};

void cuda_calculate_inner(long bigsize)
{
    thrust::counting_iterator<int> first(0);
    thrust::for_each_n(first, bigsize,
                     FCalculateInner()); // calculate inner val
}

struct FCalculateEdge
{
    __device__
    real operator()(int offset) {
        calculate(offset);
        return arr[offset];
    }
};

void cuda_calculate_edges(LongDVector &dEdgeIndices, RealDVector &dEdgeArray)
{
    thrust::transform(dEdgeIndices.begin(), dEdgeIndices.end(),
                      dEdgeArray.begin(), FCalculateEdge());
}

void cuda_shift_arrays(RealDVector &dArray,
                       RealDVector &dArrayP,
                       RealDVector &dArrayPP)
{
    dArrayP.swap(dArrayPP);
    dArrayP = dArray;
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

    LOAD_CONST_MEM(i_first, i0_);
    LOAD_CONST_MEM(j_first, j0_);
    LOAD_CONST_MEM(k_first, k0_);

    LOAD_CONST_MEM(bigsize, bigsize_);

    LOAD_CONST_MEM(i_count, ic_);
    LOAD_CONST_MEM(j_count, jc_);
    LOAD_CONST_MEM(k_count, kc_);

    LOAD_CONST_MEM(edge_i_first, ei_);
    LOAD_CONST_MEM(edge_j_first, ej_);
    LOAD_CONST_MEM(edge_k_first, ek_);

    LOAD_CONST_MEM(edge_i_count, eic_);
    LOAD_CONST_MEM(edge_j_count, ejc_);
    LOAD_CONST_MEM(edge_k_count, ekc_);

    LOAD_CONST_MEM(step_i, hi_);
    LOAD_CONST_MEM(step_j, hj_);
    LOAD_CONST_MEM(step_k, hk_);

    LOAD_CONST_MEM(step_t, ht_);

    LOAD_CONST_MEM(array, array_);
    LOAD_CONST_MEM(arrayP, arrayP_);
    LOAD_CONST_MEM(arrayPP, arrayPP_);

    const real *arr_ = array_->data().get();
    const real *arrP_ = arrayP_->data().get();
    const real *arrPP_ = arrayPP_->data().get();
    LOAD_CONST_MEM(arr, arr_);
    LOAD_CONST_MEM(arrP, arrP_);
    LOAD_CONST_MEM(arrPP, arrPP_);


    LOAD_CONST_MEM(edgeArray, edgeArray_);

    LOAD_CONST_MEM(edgeIndices, edgeIndices_);
}


