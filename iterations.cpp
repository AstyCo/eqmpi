#include "iterations.hpp"


void Requests::append(ConnectionDirection cdir, uint sz)
{
    Info i;
    i.dir = cdir;

    iv.push_back(i);

    v.push_back(MPI_Request);

    buff.push_back(RealVector());
    buff.back().reserve(sz);
}

Iterations::Iterations(const ComputeNode &n)
{
    fill(n);
}

uint Iterations::dir_size(ConnectionDirection cdir)
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

void Iterations::fill(const ComputeNode &n)
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

void Iterations::copy_recv(RealVector &v, uint i, uint j, uint k)
{
    arrayP[get_p_index(i, j, k)] = v[j * kc + k];
}

void Iterations::copy_send(RealVector &v, uint i, uint j, uint k)
{
    v[j * kc + k] = arrayP[get_p_index(i, j, k)];
}

void Iterations::copy_data(uint id, CopyMFuncPtr f)
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

void Iterations::calculate(ConnectionDirection cdir)
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

void Iterations::calculate(uint i, uint j, uint k)
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
void Iterations::it_for_each(IndexesMFuncPtr func)
{
   for (uint i = 0; i < ic; ++i) {
       for (uint j = 0; j < jc; ++j) {
           for (uint k = 0; k < kc; ++k) {
               this->*func(i, j, k);
           }
       }
   }
}

void Iterations::next_step()
{
    // array -> arrayP, arrayP -> arrayPP
    memcpy(arrayPP.data(), arrayP.data(), arrayP.size() * sizeof(real));
    memcpy(arrayPP.data(), arrayP.data(), arrayP.size() * sizeof(real));
}

uint Iterations::get_index(uint i, uint j, uint k) const
{
    return get_p_index(i,j,k);
//        return (i * jc + j) * kc + k;
}

uint Iterations::get_p_index(uint i, uint j, uint k) const
{
    return ((i + 1) * (jc + 2) + (j + 1)) * (kc + 2) + (k + 1);
}

void Iterations::set_0th(uint i, uint j, uint k)
{
    uint index = get_index(i, j, k);
    arrayPP[index] = phi(x(i), y(j), z(k));
}


void Iterations::step0()
{
    it_for_each(set_0th);
}

void Iterations::set_1th(uint i, uint j, uint k)
	{
	uint index = get_index(i, j, k);
	arrayP[index] = arrayPP[index] + ht * ht / 2 * div_grad_phi(x(i), y(j), z(k));
	}

void Iterations::step1()
{
    it_for_each(set_1th)
}
