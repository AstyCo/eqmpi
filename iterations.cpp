#include "iterations.hpp"

#include <algorithm>

#include <string.h>
//#include <mem.h>

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

void Iterations::prepareEdgeIndiceArray()
{
    MY_ASSERT(hEdgeIndices.empty());
    for (int i = 0; i < DIR_SIZE; ++i) {
        ConnectionDirection cdir = static_cast<ConnectionDirection>(i);
        if (!cnode.hasNeighbor(cdir))
            continue;
        switch (cdir) {
        case DIR_X:
        case DIR_MINUS_X:
        {
            uint i = sendEdgeId(cdir);
            for (uint j = edgeJ; j < edgeJL; ++j) {
                for (uint k = edgeK; k < edgeKL; ++k)
                    hEdgeIndices.push_back(get_index(i, j, k));
            }
            break;
        }
        case DIR_Y:
        case DIR_MINUS_Y:
        {
            uint j = sendEdgeId(cdir);
            for (uint i = edgeI; i < edgeIL; ++i) {
                for (uint k = edgeK; k < edgeKL; ++k)
                    hEdgeIndices.push_back(get_index(i, j, k));
            }
            break;
        }
        case DIR_Z:
        case DIR_MINUS_Z:
        {
            uint k = sendEdgeId(cdir);
            for (uint i = edgeI; i < edgeIL; ++i) {
                for (uint j = edgeJ; j < edgeJL; ++j)
                    hEdgeIndices.push_back(get_index(i, j, k));
            }
            break;
        }
        case DIR_Y_PERIOD_FIRST:
        {
            uint j = recvEdgeId(cdir);
            for (uint i = edgeI; i < edgeIL; ++i) {
                for (uint k = edgeK; k < edgeKL; ++k)
                    hEdgeIndices.push_back(get_index(i, j, k));
            }
            break;
        }
        case DIR_Y_PERIOD_LAST:
            // do nothing
            break;
        default: MY_ASSERT(false);
        }
    }
    totalEdgeSize = hEdgeIndices.size();

    copy_h_to_d(dEdgeIndices, hEdgeIndices);
}

void Iterations::prepareEdgeIndices()
{
    edgeI = edgeJ = edgeK = 0;
    edgeIL = ic;
    edgeJL = jc;
    edgeKL = kc;
    for (uint i = 0; i <= DIR_MINUS_Z; ++i) {
        ConnectionDirection cdir = static_cast<ConnectionDirection>(i);
        if (cnode.hasNeighbor(cdir))
            continue;
        switch (cdir) {
        case DIR_X:
            --edgeIL;
            break;
        case DIR_MINUS_X:
            ++edgeI;
            break;
        case DIR_MINUS_Y:
            ++edgeJ;
            break;
        case DIR_Z:
            --edgeKL;
            break;
        case DIR_MINUS_Z:
            ++edgeK;
            break;
        default:
            return;
        }
    }
}

void Iterations::run()
{
    MY_ASSERT(next_step == 2);
	// STEPS
    for (; next_step < clargs.K + 1; ++next_step) {
        Profiler p;

        if (next_step > 2) {
            async_recv_all();
            async_send_all();
        }
        get_time(times.mpi_send_recv, p);

        cuda_calculate_inner(bigsize);

        get_time(times.parallel_cycles, p);

        // sync recv prev
        MPI_Waitall(recv_requests.v.size(),
                    recv_requests.v.data(),
                    recv_requests.statuses.data());

        for (uint i = 0; i < recv_requests.size(); ++i)
            copy_data(recv_requests, i, MPI_OP_RECV);

        get_time(times.mpi_send_recv, p);

        copy_edges_to_d();

        get_time(times.host_device_exchange, p);

        cuda_calculate_edges(dEdgeIndices, dEdgeArray);

        get_time(times.parallel_cycles, p);

        copy_edges_to_h();

        get_time(times.host_device_exchange, p);

        if (clargs.deviation)
            printDeviations(next_step);

        if (next_step < clargs.K) {
            p.start();
            cuda_shift_arrays(dArray, dArrayP, dArrayPP);
            get_time(times.shift_arrays, p);
        }
    } // ENDS STEPS
}

void Iterations::copy_edges_to_h()
{
    cuda_copy_from_dvector(dEdgeIndices, dEdgeArray);
    copy_d_to_h(dEdgeArray, hEdgeArray);
    for (uint i = 0; i < hEdgeIndices.size(); ++i) {
        long offset = hEdgeIndices[i];
        hArrayBuff[offset] = hEdgeArray[i];
    }
}

void Iterations::copy_edges_to_d()
{
    for (uint i = 0; i < hEdgeIndices.size(); ++i) {
        long offset = hEdgeIndices[i];
        hEdgeArray[i] = hArrayBuff[offset];
    }

    copy_h_to_d(hEdgeArray, dEdgeArray);
    cuda_copy_to_dvector(dEdgeIndices, dEdgeArray);
}

void Iterations::async_send_all()
{
    // asynchronous send to every neighbor processor
    for (int i = 0; i < send_requests.size(); ++i) {
        copy_data(send_requests, i, MPI_OP_SEND);

        MPI_Isend(send_requests.buffs[i].data(),
                  send_requests.buffs[i].size(),
                  MPI_TYPE_REAL,
                  cnode.neighbor(send_requests.iv[i].dir),
                  send_requests.iv[i].dir,
                  MPI_COMM_WORLD,
                  &send_requests.v[i]);
    }
}

void Iterations::async_recv_all()
{
    // asynchronous recv from every neighbor processor
    for (uint i = 0; i < recv_requests.size(); ++i) {
        MPI_Irecv(recv_requests.buffs[i].data(),
                  recv_requests.buffs[i].size(),
                  MPI_TYPE_REAL,
                  cnode.neighbor(recv_requests.iv[i].dir),
                  CDPair(recv_requests.iv[i].dir),
                  MPI_COMM_WORLD,
                  &recv_requests.v[i]);
    }
}

void Iterations::printDeviations(uint n)
{
    real avgDeviation = cuda_get_local_avg_deviation(bigsize,
                                                     ic * jc * kc,
                                                     time(n),
                                                     dDeviationsArray);
    real globalDeviation = 0;
    MPI_Reduce(&avgDeviation, &globalDeviation, 1, MPI_TYPE_REAL,
               MPI_SUM, 0, MPI_COMM_WORLD);

    globalDeviation /= cnode.mpi.procCount;

    if (cnode.mpi.rank == 0) {
        std::cout << SSTR("%%%," << cnode.scTag()
                         << ',' << cnode.mpi.procCount
                         << ',' << N
                         << ',' << next_step
                         << ',' << globalDeviation) << std::endl;
    }

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

    int iMissedItemCount = N % n.gridDimensions.x;
    int jMissedItemCount = N % n.gridDimensions.y;
    int kMissedItemCount = N % n.gridDimensions.z;

    i0 = MIN(n.x, iMissedItemCount) * (ic + 1) + MAX(n.x - iMissedItemCount, 0) * ic;
    j0 = MIN(n.y, jMissedItemCount) * (jc + 1) + MAX(n.y - jMissedItemCount, 0) * jc;
    k0 = MIN(n.z, kMissedItemCount) * (kc + 1) + MAX(n.z - kMissedItemCount, 0) * kc;

    if (cnode.x < iMissedItemCount)
        ++ic;
    if (cnode.y < jMissedItemCount)
        ++jc;
    if (cnode.z < kMissedItemCount)
        ++kc;

    hx = VAL_LX / ic;
    hy = VAL_LY / jc;
    hz = VAL_LZ / kc;


    ht = VAL_T / clargs.K;

    bigsize = (ic + 2) * (jc + 2) * (kc + 2);

    prepareEdgeIndices();
    prepareEdgeIndiceArray();

    // optimization (allocations may throw std::bad_alloc if no enough memory)
    try {
        hArrayBuff.resize(bigsize);
        MPI_Barrier(MPI_COMM_WORLD);
        Profiler p_allocations;
        cuda_resize(dArray, dArrayP, dArrayPP,
                    dEdgeArray, hEdgeArray,
                    dDeviationsArray,
                    totalEdgeSize, bigsize);
        get_time(times.allocations, p_allocations);

        for (int i = 0; i < DIR_SIZE; ++i) {
            ConnectionDirection cdir = toCD(i);
            if (n.hasNeighbor(cdir)) {
                Requests::Info info;
                info.dir = cdir;

                send_requests.append(info, dir_size(cdir));
                recv_requests.append(info, dir_size(cdir));
            }
        }
    }
    catch(...) {
        MY_ASSERT(false);
    }

    Profiler p_allocations;
    cuda_load_const_mem(N,
                        i0, j0, k0,
                        ic, jc, kc,
                        hx, hy, hz,
                        ht,
                        &dArray, &dArrayP, &dArrayPP);
    get_time(times.allocations, p_allocations);
}

void Iterations::copy_recv(RealVector &v, RealVector &a,
                           int i, int j, int k, long offset)
{
    CHECK_INDEX(offset, 0, v.size());
    CHECK_INDEX(get_index(i, j, k), 0, a.size());

    a[get_index(i, j, k)] = v[offset];
}

void Iterations::copy_send(RealVector &v, RealVector &a,
                           int i, int j, int k, long offset)
{
    CHECK_INDEX(offset, 0, v.size());
    CHECK_INDEX(get_index(i, j, k), 0, a.size());

    v[offset] = a[get_index(i, j, k)];
}

void Iterations::copy_data(Requests &requests, uint id, MPI_OP type)
{
    CHECK_INDEX(id, 0, requests.iv.size());
    CopyMFuncPtr copy_func = ((type == MPI_OP_SEND) ? &Iterations::copy_send
                                        : &Iterations::copy_recv);
    ConnectionDirection cdir= requests.iv[id].dir;
    RealVector &v = requests.buffs[id];

    long offset = 0;
    switch (cdir) {
    case DIR_X:
    case DIR_MINUS_X:
    {
        int i = edgeId(cdir, type);
        for (uint j = 0; j < jc; ++j) {
            for (uint k = 0; k < kc; ++k) {
                (this->*copy_func)(v, hArrayBuff, i, j, k, offset++);
            }
        }
        break;
    }
    case DIR_Y:
    case DIR_MINUS_Y:
    case DIR_Y_PERIOD_FIRST:
    case DIR_Y_PERIOD_LAST:
    {
        int j = edgeId(cdir, type);
        for (uint i = 0; i < ic; ++i) {
            for (uint k = 0; k < kc; ++k) {
                (this->*copy_func)(v, hArrayBuff, i, j, k, offset++);
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
                (this->*copy_func)(v, hArrayBuff, i, j, k, offset++);
            }
        }
        break;
    }
    default: MY_ASSERT(false);
    }
}

long Iterations::get_index(uint i, uint j, uint k) const
{
    return (long(i + 1) * (jc + 2) + (j + 1)) * (kc + 2) + (k + 1);
}

void Iterations::step0()
{
    MY_ASSERT(next_step == 0);
    cuda_step_0(dArray, dArrayPP);
    next_step = 1;
}

void Iterations::step1()
{
    MY_ASSERT(next_step == 1);
    cuda_step_1(dArrayP);
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
