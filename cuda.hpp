#include "globals.hpp"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

typedef thrust::device_vector<real> RealDVector;
typedef thrust::host_vector<real> RealHVector;

typedef thrust::host_vector<long> LongHVector;
typedef thrust::device_vector<long> LongDVector;

void cuda_resize(RealDVector &dArray, RealDVector &dArrayP, RealDVector &dArrayPP, RealDVector &dEdgeArray, RealHVector &hEdgeArray, RealDVector &analyticalSolution, long totalEdgeSize, long bigsize);

void cuda_step_0(RealDVector &darray, RealDVector &dArrayPP);
void cuda_step_1(RealDVector &dArrayP);

void cuda_calculate_inner(long bigsize);

void cuda_calculate_edges(LongDVector &dEdgeIndices, RealDVector &dEdgeArray);

void cuda_shift_arrays(RealDVector &dArray,
                       RealDVector &dArrayP,
                       RealDVector &dArrayPP);

void cuda_load_const_mem(int N_,
                         int i0_, int j0_, int k0_,
                         long long bigsize_,
                         long ic_, long jc_, long kc_,
                         int ei_, int ej_, int ek_,
                         int eic_, int ejc_, int ekc_,
                         real hi_, real hj_, real hk_,
                         real ht_, RealDVector *array_,
                         RealDVector *arrayP_,
                         RealDVector *arrayPP_,
                         RealDVector *edgeArray_, LongDVector *edgeIndices_);

template <typename THVec, typename TDVec>
void copy_h_to_d(THVec &hv, TDVec &dv)
{
    dv = hv;
}

template <typename THVec, typename TDVec>
void copy_d_to_h(TDVec &dv, THVec &hv)
{
    hv = dv;
}

