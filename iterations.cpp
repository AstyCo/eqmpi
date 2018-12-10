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
//    printDeviationsPrivate(arrayPP, 0);
    step1();
//    printDeviationsPrivate(arrayP, 1);
}


void Iterations::run()
{
    MY_ASSERT(next_step == 2);
	// STEPS
    for (; next_step < clargs.K + 1; ++next_step) {
        profiler.step();
        cnode.print(SSTR("ITER " << next_step << ',' << profiler.time()));
        async_recv_all();
        async_send_all();

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
        {
            int request_index;
            MPI_Status status;

            for(;;) {
                MPI_Waitany(recv_requests.v.size(),
                            recv_requests.v.data(),
                            &request_index,
                            &status);
                if (request_index == MPI_UNDEFINED) {
                    break;
                }
                CHECK_INDEX(request_index, 0, recv_requests.iv.size());
                const Iterations::Requests::Info &info
                        = recv_requests.iv[request_index];
                copy_data(recv_requests, request_index, MPI_OP_RECV);
                calculate(info.dir);
            }
        }
        calculate_edge_values();

        if (clargs.deviation) {
            prepareSolution(next_step);
            printDeviations(next_step);
        }
        if (next_step < clargs.K)
            shift_arrays(); // update arrays
    } // ENDS STEPS
}

void Iterations::async_send_all()
{
    MPI_Waitall(send_requests.v.size(),
                send_requests.v.data(),
                send_requests.statuses.data()); // wait for every asynchronous send (so we can use send buffers again)
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

void Iterations::prepareSolution(uint n)
{
#ifdef WITH_OMP
#   pragma omp parallel for
#endif
   for (int i = 0; i < ic; ++i) {
       for (uint j = 0; j < jc; ++j) {
           for (uint k = 0; k < kc; ++k) {
               long id = get_index(i, j, k);
               analyticalSolution[id] = u(x(i), y(j), z(k), time(n));
           }
       }
   }
}

real Iterations::getDeviation(const Iterations::RealVector &arr, uint i, uint j, uint k, uint n) const
{
    long long id = get_index(i, j, k);
    real correctAnswer = analyticalSolution[id];
    real approxAnswer = arr[id];

    return ABS(correctAnswer - approxAnswer);
}

void Iterations::printDeviations(uint n)
{
    printDeviationsPrivate(array, n);
}

void Iterations::printDeviationsPrivate(const Iterations::RealVector &arr, uint n)
{
    prepareSolution(n);

    real avgDeviation = 0;
    real corr = 0, appr = 0;

#ifdef WITH_OMP
#   pragma omp parallel for
#endif
    for (int i = 0; i < ic; ++i) {
        for (uint j = 0; j < jc; ++j) {
            for (uint k = 0; k < kc; ++k) {
                long long id = get_index(i, j, k);
                real correctAnswer = analyticalSolution[id];
                real approxAnswer = arr[id];

                corr += correctAnswer;
                appr += approxAnswer;

                avgDeviation += getDeviation(arr, i, j, k, n);
            }
        }
    }
    real globalDeviation=0;
    MPI_Reduce(&avgDeviation, &globalDeviation, 1, MPI_TYPE_REAL,
               MPI_SUM, 0, MPI_COMM_WORLD);
    globalDeviation /= N*N*N;

//    long size = ic * jc * kc;
//    avgDeviation /= size;
//    if (avgDeviation > 0.01) {
//        cnode.print(SSTR("local delta for step " << n << " equals "
//                         << avgDeviation << " app " << appr << " corr " << corr));
//    }

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

    // optimization (allocations may throw std::bad_alloc if no enough memory)
    try {
        prepareEdgeIndices();
        array.resize(bigsize);
        arrayP.resize(bigsize);
        arrayPP.resize(bigsize);
        if (clargs.deviation)
            analyticalSolution.resize(bigsize);

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
}

void Iterations::copy_recv(RealVector &v, RealVector &a,
                           int i, int j, int k, uint offset)
{
    CHECK_INDEX(offset, 0, v.size());
    CHECK_INDEX(get_index(i, j, k), 0, a.size());

    a[get_index(i, j, k)] = v[offset];
}

void Iterations::copy_send(RealVector &v, RealVector &a,
                           int i, int j, int k, uint offset)
{
    CHECK_INDEX(offset, 0, v.size());
    CHECK_INDEX(get_index(i, j, k), 0, a.size());

    v[offset] = a[get_index(i, j, k)];
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
    {
        uint j = recvEdgeId(cdir);
        for (uint i = 1; i < ic - 1; ++i) {
            for (uint k = 1; k < kc - 1; ++k)
                calculate(i, j, k);
        }
        break;
    }
    case DIR_Y_PERIOD_LAST:
        // just copy
        break;
    default: MY_ASSERT(false);
    }
}

void Iterations::calculate_edge_values()
{
    for (long i = 0; i < edgeIndeces.size(); ++i) {
        const Indice &ind = edgeIndeces[i];
        calculate(ind.i, ind.j, ind.k);
    }
}

void Iterations::calculate(uint i, uint j, uint k)
{
    long index = get_index(i, j, k);

//    array[index]++;
//    cnode.print(SSTR("calculate " << i << ',' << j << ',' << k));
//    return;

    CHECK_INDEX(index, 0, array.size());
    CHECK_INDEX(index, 0, array.size());
    CHECK_INDEX(get_index(i-1,j,k), 0, array.size());
    CHECK_INDEX(get_index(i+1,j,k), 0, array.size());
    CHECK_INDEX(get_index(i,j-1,k), 0, array.size());
    CHECK_INDEX(get_index(i,j+1,k), 0, array.size());
    CHECK_INDEX(get_index(i,j,k-1), 0, array.size());
    CHECK_INDEX(get_index(i,j,k+1), 0, array.size());


    array[index] = 2 * arrayP[index] - arrayPP[index]
            + ht * ht * (
                (arrayP[get_index(i-1,j,k)]
                - 2 * arrayP[index]
                + arrayP[get_index(i+1,j,k)]) / hx / hx
            + (arrayP[get_index(i,j-1,k)]
            - 2 * arrayP[index]
            + arrayP[get_index(i,j+1,k)]) / hy / hy
            + (arrayP[get_index(i,j,k-1)]
            - 2 * arrayP[index]
            + arrayP[get_index(i,j,k+1)]) / hz / hz
            );
}

void Iterations::it_for_each(IndexesMFuncPtr func)
{
#ifdef WITH_OMP
#   pragma omp parallel for
#endif
   for (int i = 0; i < ic; ++i) {
       for (int j = 0; j < jc; ++j) {
           for (int k = 0; k < kc; ++k) {
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


void Iterations::prepareEdgeIndices()
{
    std::vector<int> tmp;
    tmp.resize(ic * jc * kc);
    for (long i = 0; i < tmp.size(); ++i)
        tmp[i] = 0;

    for (uint i = 0; i < ic; ++i) {
        tmp[get_exact_index(i, 0, 0)] = 1;
        tmp[get_exact_index(i, 0, kc - 1)] = 1;
        tmp[get_exact_index(i, jc - 1, 0)] = 1;
        tmp[get_exact_index(i, jc - 1, kc - 1)] = 1;
    }
    for (uint j = 0; j < jc; ++j) {
        tmp[get_exact_index(0, j, 0)] = 1;
        tmp[get_exact_index(0, j, kc - 1)] = 1;
        tmp[get_exact_index(ic - 1, j, 0)] = 1;
        tmp[get_exact_index(ic - 1, j, kc - 1)] = 1;
    }
    for (uint k = 0; k < kc; ++k) {
        tmp[get_exact_index(0, 0, k)] = 1;
        tmp[get_exact_index(0, jc - 1, k)] = 1;
        tmp[get_exact_index(ic - 1, 0, k)] = 1;
        tmp[get_exact_index(ic - 1, jc - 1, k)] = 1;
    }

    for (uint i = 0; i < DIR_Y_PERIOD_FIRST; ++i) {
        ConnectionDirection cdir = static_cast<ConnectionDirection>(i);
        if (cnode.hasNeighbor(cdir))
            continue;
        switch (cdir) {
        case DIR_X:
            for (uint j = 0; j < jc; ++j) {
                tmp[get_exact_index(ic - 1, j, 0)] = 0;
                tmp[get_exact_index(ic - 1, j, kc - 1)] = 0;
            }
            for (uint k = 0; k < kc; ++k) {
                tmp[get_exact_index(ic - 1, 0, k)] = 0;
                tmp[get_exact_index(ic - 1, jc - 1, k)] = 0;
            }
            break;
        case DIR_MINUS_X:
            for (uint j = 0; j < jc; ++j) {
                tmp[get_exact_index(0, j, 0)] = 0;
                tmp[get_exact_index(0, j, kc - 1)] = 0;
            }
            for (uint k = 0; k < kc; ++k) {
                tmp[get_exact_index(0, 0, k)] = 0;
                tmp[get_exact_index(0, jc - 1, k)] = 0;
            }
            break;
        case DIR_Y:
            if (cnode.hasNeighbor(DIR_Y_PERIOD_FIRST))
                continue;
            for (uint i = 0; i < ic; ++i) {
                tmp[get_exact_index(i, jc - 1, 0)] = 0;
                tmp[get_exact_index(i, jc - 1, kc - 1)] = 0;
            }
            for (uint k = 0; k < kc; ++k) {
                tmp[get_exact_index(0, jc - 1, k)] = 0;
                tmp[get_exact_index(ic - 1, jc - 1, k)] = 0;
            }
            break;
        case DIR_MINUS_Y:
            for (uint i = 0; i < ic; ++i) {
                tmp[get_exact_index(i, 0, 0)] = 0;
                tmp[get_exact_index(i, 0, kc - 1)] = 0;
            }
            for (uint k = 0; k < kc; ++k) {
                tmp[get_exact_index(0, 0, k)] = 0;
                tmp[get_exact_index(ic - 1, 0, k)] = 0;
            }
            break;
        case DIR_Z:
            for (uint i = 0; i < ic; ++i) {
                tmp[get_exact_index(i, 0, kc - 1)] = 0;
                tmp[get_exact_index(i, jc - 1, kc - 1)] = 0;
            }
            for (uint j = 0; j < jc; ++j) {
                tmp[get_exact_index(0, j, kc - 1)] = 0;
                tmp[get_exact_index(ic - 1, j, kc - 1)] = 0;
            }
            break;
        case DIR_MINUS_Z:
            for (uint i = 0; i < ic; ++i) {
                tmp[get_exact_index(i, 0, 0)] = 0;
                tmp[get_exact_index(i, jc - 1, 0)] = 0;
            }
            for (uint j = 0; j < jc; ++j) {
                tmp[get_exact_index(0, j, 0)] = 0;
                tmp[get_exact_index(ic - 1, j, 0)] = 0;
            }
            break;
        default:
            MY_ASSERT(false);
            break;
        }
    }

    for (long i = 0; i < tmp.size(); ++i) {
        if (!tmp[i])
            continue;
        long vz = i % kc;
        long vxy = i / kc;
        long vy = vxy % jc;
        long vx = vxy / jc;

        edgeIndeces.push_back(Indice(vx, vy, vz));
    }
}

long Iterations::get_index(uint i, uint j, uint k) const
{
    return (long(i + 1) * (jc + 2) + (j + 1)) * (kc + 2) + (k + 1);
}

long Iterations::get_exact_index(uint i, uint j, uint k) const
{
    return (long(i)*jc + j) * kc + k;
}

void Iterations::set_0th(uint i, uint j, uint k)
{
    long index = get_index(i, j, k);
    CHECK_INDEX(index, 0, arrayPP.size());
    arrayPP[index] = phi(x(i), y(j), z(k));
}


void Iterations::step0()
{
    MY_ASSERT(next_step == 0);
    it_for_each(&Iterations::set_0th);
    memcpy(array.data(), arrayPP.data(), bigsize * sizeof(real));
    next_step = 1;
}

void Iterations::set_1th(uint i, uint j, uint k)
{
    long index = get_index(i, j, k);
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

void Iterations::printArrayDebug()
{
    if (cnode.mpi.rank != cnode.mpi.procCount - 1)
        return;
    cnode.print(SSTR("STEP " << next_step));
    for (uint i = ic - 4; i < ic; ++i) {
        for (uint j = jc - 4; j < jc; ++j) {
            for (uint k = kc - 4; k < kc; ++k) {
                cnode.print(SSTR('(' << i << ',' << j << ',' << k << ')'
                                 << ' ' << array[get_index(i,j,k)]));
            }
        }
    }
}
