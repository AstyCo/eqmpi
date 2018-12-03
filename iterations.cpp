#include "iterations.hpp"

real Iterations::T = 10;
real Iterations::Lx = 2 * M_PI;
real Iterations::Ly = 2 * M_PI;
real Iterations::Lz = 2 * M_PI;

void Iterations::Requests::append(Iterations::Requests::Info info, uint sz)
{
    iv.push_back(info);

    v.push_back(MPI_REQUEST_NULL);

    buffs.push_back(RealVector());
    buffs.back().resize(sz);

    statuses.push_back(MPI_Status());
}

Iterations::Iterations(uint N_)
    : N(N_), next_step(0)
{
    fill(cnode);
}

void Iterations::prepare()
{
	step0();
	step1();
}


void Iterations::run()
{
    MY_ASSERT(next_step == 2);
	// STEPS
    bool first = true;
    for (; next_step < clargs.K + 1; ++next_step) {
//        cnode.print(SSTR("STEP " << next_step));
        // async receive prev
        if (!first) {
            // asynchronous recv from every neighbor processor
            for (uint i = 0; i < recv_requests.size(); ++i) {
//                cnode.print(SSTR("RECV FROM " << cnode.neighbor(recv_requests.iv[i].dir)));
                MPI_Irecv(recv_requests.buffs[i].data(),
                          recv_requests.buffs[i].size(),
                          MPI_TYPE_REAL,
                          cnode.neighbor(recv_requests.iv[i].dir),
                          TAG_BOUNDARY_ELEMENTS,
                          MPI_COMM_WORLD,
                          &recv_requests.v[i]);
            }
        }
#ifdef WITH_OMP
#       pragma omp parallel for // parallel slices
#endif
        for (int i = 1; i < ic - 1; ++i) {
            for (uint j = 1; j < jc - 1; ++j) {
                for (uint k = 1; k < kc - 1; ++k) {
                    calculate(i, j, k);
                }
            }
        }
        // sync recv prev
//#       ifdef WITH_OMP
//#           pragma omp parallel // parallel recv
//#       endif
        {
            int request_index;
            MPI_Status status;
//            cnode.print("MPI_Waitany Irecv");
            for(;;) {
                MPI_Waitany(recv_requests.v.size(),
                            recv_requests.v.data(),
                            &request_index,
                            &status);
                if (request_index == MPI_UNDEFINED) {
//                    cnode.print(SSTR("UNDEFINED"));
                    break;
                }
                CHECK_INDEX(request_index, 0, recv_requests.iv.size());
//                cnode.print(
//                    SSTR("MESSAGE FROM " <<
//                         cnode.neighbor(recv_requests.iv[request_index].dir)));
                const Iterations::Requests::Info &info
                        = recv_requests.iv[request_index];
//                cnode.print(SSTR("info " << info.dir <<',' << request_index));
                copy_data(recv_requests, request_index, MPI_OP_RECV);
                calculate(info.dir);
            }
        }
        calculate_edge_values();
        shift_arrays(); // update arrays
        // async send prev
        if (next_step < clargs.K) {
//            cnode.print("MPI_Waitall");
            MPI_Waitall(send_requests.v.size(),
                        send_requests.v.data(),
                        send_requests.statuses.data()); // wait for every asynchronous send (so we can use send buffers again)
            // asynchronous send to every neighbor processor
//#           ifdef WITH_OMP
//#               pragma omp parallel for // parallel recv
//#           endif
            for (int i = 0; i < send_requests.size(); ++i) {
//                cnode.print(SSTR("SEND TO " << cnode.neighbor(send_requests.iv[i].dir)));
                copy_data(send_requests, i, MPI_OP_SEND);

                MPI_Isend(send_requests.buffs[i].data(),
                          send_requests.buffs[i].size(),
                          MPI_TYPE_REAL,
                          cnode.neighbor(send_requests.iv[i].dir),
                          TAG_BOUNDARY_ELEMENTS,
                          MPI_COMM_WORLD,
                          &send_requests.v[i]);
            }
        }
        first = false;
    } // ENDS STEPS
}


uint Iterations::dir_size(ConnectionDirection cdir)
{
    switch (cdir) {
    case DIR_X:
    case DIR_MINUS_X:
        return jc * kc;
    case DIR_Y:
    case DIR_MINUS_Y:
    case DIR_Y_PERIOD_FIRST:
    case DIR_Y_PERIOD_LAST:
        return ic * kc;
    case DIR_Z:
    case DIR_MINUS_Z:
        return ic * jc;
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

    ic = N / n.gridDimensions.x;
    jc = N / n.gridDimensions.y;
    kc = N / n.gridDimensions.z;

    // TODO
    MY_ASSERT(0 == N % n.gridDimensions.x);
    MY_ASSERT(0 == N % n.gridDimensions.y);
    MY_ASSERT(0 == N % n.gridDimensions.z);

    i0 = n.x * ic;
    j0 = n.y * jc;
    k0 = n.z * kc;

    imax = i0 + ic;
    jmax = j0 + jc;
    kmax = k0 + kc;

    hx = Lx / ic;
    hy = Ly / jc;
    hz = Lz / kc;

    ht = T / clargs.K;

    bigsize = (ic + 2) * (jc + 2) * (kc + 2);

    // optimization (allocations may throw std::bad_alloc if no enough memory)
    try {
    array.resize(bigsize);
    arrayP.resize(bigsize);
    arrayPP.resize(bigsize);

    for (int i = 0; i < DIR_SIZE; ++i) {
        ConnectionDirection cdir = toCD(i);
        if (n.is(cdir)) {
            Requests::Info info;
            info.dir = cdir;

            send_requests.append(info, dir_size(cdir));
            recv_requests.append(info, dir_size(cdir));
        }
    }
//    cnode.print(SSTR("(x, y, z): (" << n.x << ',' << n.y << ',' << n.z << ')'
//                     << " i0, j0, k0 (ic, jc, kc): "
//                     << i0 << ',' << j0 << ',' << k0
//                     << " (" << ic  << ',' << jc << ',' << kc << ")"
//                     << std::endl));
    }
    catch(...) {
        MY_ASSERT(false);
    }
}

void Iterations::copy_recv(RealVector &v, RealVector &a,
                           int i, int j, int k, uint offset)
{
    CHECK_INDEX(offset, 0, v.size());
    CHECK_INDEX(get_p_index(i, j, k), 0, a.size());

    a[get_p_index(i, j, k)] = v[offset];
//    cnode.print(SSTR("copy_recv " << i << ',' << j << ',' << k
//                     << '(' << get_p_index(i, j, k) << ')'
//                     << '=' << offset));
}

void Iterations::copy_send(RealVector &v, RealVector &a,
                           int i, int j, int k, uint offset)
{
    CHECK_INDEX(offset, 0, v.size());
    CHECK_INDEX(get_p_index(i, j, k), 0, a.size());

    v[offset] = a[get_p_index(i, j, k)];
//    cnode.print(SSTR("copy_send " << offset
//                     << '=' << i << ',' << j << ',' << k
//                     << '(' << get_p_index(i, j, k) << ')'));
}

void Iterations::copy_data(Requests &requests, uint id, MPI_OP type)
{
    CHECK_INDEX(id, 0, requests.iv.size());
    CopyMFuncPtr copy_func = ((type == MPI_OP_SEND) ? &Iterations::copy_send
                                        : &Iterations::copy_recv);
    ConnectionDirection cdir= requests.iv[id].dir;
    RealVector &v = requests.buffs[id];
    uint offset = 0;
    switch (cdir) {
    case DIR_X:
    case DIR_MINUS_X:
    {
        int i = edgeId(cdir, type);
        for (uint j = 0; j < jc; ++j) {
            for (uint k = 0; k < kc; ++k) {
                (this->*copy_func)(v, arrayP, i, j, k, offset++);
            }
        }
        break;
    }
    case DIR_Y:
    case DIR_MINUS_Y:
    case DIR_Y_PERIOD_FIRST:
    case DIR_Y_PERIOD_LAST:
    {
        RealVector &a = (((cdir == DIR_Y_PERIOD_LAST)
                         && (type == MPI_OP_RECV))
                         ? array
                         : arrayP);
        int j = edgeId(cdir, type);
        for (uint i = 0; i < ic; ++i) {
            for (uint k = 0; k < kc; ++k) {
                (this->*copy_func)(v, a, i, j, k, offset++);
            }
        }
        break;
    }
    case DIR_Z:
    case DIR_MINUS_Z:
    {
        int k = edgeId(cdir, type);
        for (uint i = 0; i < ic; ++i) {
            for (uint j = 0; j < jc; ++j) {
                (this->*copy_func)(v, arrayP, i, j, k, offset++);
            }
        }
        break;
    }
    default: MY_ASSERT(false);
    }
}

void Iterations::calculate(ConnectionDirection cdir)
{
    switch (cdir) {
    case DIR_X:
    case DIR_MINUS_X:
    {
        uint i = sendEdgeId(cdir);
        for (uint j = 1; j < jc - 1; ++j) {
            for (uint k = 1; k < kc - 1; ++k)
                calculate(i, j, k);
        }
        break;
    }
    case DIR_Y:
    case DIR_MINUS_Y:
    {
        uint j = sendEdgeId(cdir);
        for (uint i = 1; i < ic - 1; ++i) {
            for (uint k = 1; k < kc - 1; ++k)
                calculate(i, j, k);
        }
        break;
    }
    case DIR_Z:
    case DIR_MINUS_Z:
    {
        uint k = sendEdgeId(cdir);
        for (uint i = 1; i < ic - 1; ++i) {
            for (uint j = 1; j < jc - 1; ++j)
                calculate(i, j, k);
        }
        break;
    }
    case DIR_Y_PERIOD_FIRST:
    case DIR_Y_PERIOD_LAST:
        // just copy
        break;
    default: MY_ASSERT(false);
    }
}

void Iterations::calculate_edge_values()
{
    for (uint i = 0; i < ic; ++i) {
        calculate(i, 0, 0);
        calculate(i, jc - 1, 0);
        calculate(i, 0, kc - 1);
        calculate(i, jc - 1, kc - 1);
    }
    for (uint j = 1; j < jc - 1; ++j) {
        calculate(0, j, 0);
        calculate(ic - 1, j, 0);
        calculate(0, j, kc - 1);
        calculate(ic - 1, j, kc - 1);
    }
    for (uint k = 1; k < kc - 1; ++k) {
        calculate(0, 0, k);
        calculate(ic - 1, 0, k);
        calculate(0, jc - 1, k);
        calculate(ic - 1, jc - 1, k);
    }
}

void Iterations::calculate(uint i, uint j, uint k)
{
    uint index = get_index(i, j, k);
    uint p_index = get_p_index(i, j, k);
    uint pp_index = index;

//    cnode.print(SSTR("calculate " << i << ',' << j << ',' << k));

    CHECK_INDEX(index, 0, array.size());
    CHECK_INDEX(p_index, 0, array.size());
    CHECK_INDEX(get_index(i-1,j,k), 0, array.size());
    CHECK_INDEX(get_index(i+1,j,k), 0, array.size());
    CHECK_INDEX(get_index(i,j-1,k), 0, array.size());
    CHECK_INDEX(get_index(i,j+1,k), 0, array.size());
    CHECK_INDEX(get_index(i,j,k-1), 0, array.size());
    CHECK_INDEX(get_index(i,j,k+1), 0, array.size());


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
#ifdef WITH_OMP
#   pragma omp parallel for
#endif
   for (int i = 0; i < ic; ++i) {
       for (uint j = 0; j < jc; ++j) {
           for (uint k = 0; k < kc; ++k) {
               (this->*func)(i, j, k);
           }
       }
   }
}

void Iterations::shift_arrays()
{
    uint byteSize = bigsize * sizeof(real);

//    // array -> arrayP, arrayP -> arrayPP
    memcpy(arrayPP.data(), arrayP.data(), byteSize);
    memcpy(arrayP.data(), array.data(), byteSize);
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
    CHECK_INDEX(index, 0, arrayPP.size());
    arrayPP[index] = phi(x(i), y(j), z(k));
}


void Iterations::step0()
{
    MY_ASSERT(next_step == 0);
    it_for_each(&Iterations::set_0th);
    next_step = 1;
}

void Iterations::set_1th(uint i, uint j, uint k)
{
    uint index = get_index(i, j, k);
    CHECK_INDEX(index, 0, arrayP.size());
    arrayP[index] = arrayPP[index] + ht * ht / 2 * div_grad_phi(x(i), y(j), z(k));
}

void Iterations::step1()
{
    MY_ASSERT(next_step == 1);
    it_for_each(&Iterations::set_1th);
    next_step = 2;
}

int Iterations::edgeId(ConnectionDirection cdir, Iterations::MPI_OP op_type)
{
    if (op_type == MPI_OP_RECV)
        return recvEdgeId(cdir);
    else
        return sendEdgeId(cdir);
}

int Iterations::sendEdgeId(ConnectionDirection cdir) const
{
    switch (cdir) {
    case DIR_Y_PERIOD_FIRST: return jc - 1;
    case DIR_Y_PERIOD_LAST: return 1;
    case DIR_X: return ic - 1;
    case DIR_MINUS_X: return 0;
    case DIR_Y: return jc - 1;
    case DIR_MINUS_Y:return 0;
    case DIR_Z: return kc - 1;
    case DIR_MINUS_Z: return 0;
    default: MY_ASSERT(false); return 0;
    }
}

int Iterations::recvEdgeId(ConnectionDirection cdir) const
{
    switch (cdir) {
    case DIR_Y_PERIOD_FIRST:
        return jc;
    case DIR_Y_PERIOD_LAST:
        return 0;
    case DIR_X:
    case DIR_Y:
    case DIR_Z:
        return sendEdgeId(cdir) + 1;
    case DIR_MINUS_X:
    case DIR_MINUS_Y:
    case DIR_MINUS_Z:
        return sendEdgeId(cdir) - 1;
    default: MY_ASSERT(false); return 0;
    }
}

int Iterations::sendrecvEdgeId(ConnectionDirection cdir) const
{
    switch (cdir) {
    case DIR_X: return ic - 1;
    case DIR_MINUS_X: return 0;
    case DIR_Y: return jc - 1;
    case DIR_MINUS_Y:return 0;
    case DIR_Z: return kc - 1;
    case DIR_MINUS_Z: return 0;
    default: MY_ASSERT(false); return 0;
    }
}
