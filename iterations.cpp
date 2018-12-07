#include "iterations.hpp"

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
    printDeviationsPrivate(arrayPP, 0);
	step1();
    printDeviationsPrivate(arrayP, 1);
}


void Iterations::run()
{
    MY_ASSERT(next_step == 2);
	// STEPS
    for (; next_step < clargs.K + 1; ++next_step) {
        // async receive prev
        if (next_step > 2) {
            // asynchronous recv from every neighbor processor
            for (uint i = 0; i < recv_requests.size(); ++i) {
                MPI_Irecv(recv_requests.buffs[i].data(),
                          recv_requests.buffs[i].size(),
                          MPI_TYPE_REAL,
                          cnode.neighbor(recv_requests.iv[i].dir),
                          next_step,
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
        {
            int request_index;
            MPI_Status status;
            for (uint i = 0; i < no_neighbour_edges.size(); ++i)
                calculate(no_neighbour_edges[i]);

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
        prepareSolution(next_step);
        printDeviations(next_step);
        // async send prev
        if (next_step < clargs.K) {
            shift_arrays(); // update arrays
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
                          next_step + 1,
                          MPI_COMM_WORLD,
                          &send_requests.v[i]);
            }
        }
    } // ENDS STEPS
}

void Iterations::prepareSolution(uint n)
{
#ifdef WITH_OMP
#   pragma omp parallel for
#endif
   for (int i = 0; i < ic; ++i) {
       for (uint j = 0; j < jc; ++j) {
           for (uint k = 0; k < kc; ++k) {
               uint id = get_index(i, j, k);
               analyticalSolution[id] = u(x(i), y(j), z(k), time(n));
           }
       }
   }
}

void Iterations::printDeviation(uint i, uint j, uint k, uint n)
{
    printDeviationPrivate(array, i, j, k, n);
}

void Iterations::printDeviationPrivate(const Iterations::RealVector &arr, uint i, uint j, uint k, uint n)
{
    uint id = get_index(i, j, k);
    real correctAnswer = analyticalSolution[id];
    real approxAnswer = arr[id];

    cnode.print(SSTR("print " << id << " approx " << approxAnswer));
    cnode.print(SSTR("print " << id << " u " << correctAnswer));

    cnode.print(SSTR("n=" << n << " (" << i << ',' << j << ',' << k
                     <<") DEV, CORR, APPROX: "
                     << ABS(correctAnswer - approxAnswer)
                     << ',' << ' ' << correctAnswer
                         << ',' << ' ' << approxAnswer));
}

real Iterations::getDeviation(const Iterations::RealVector &arr, uint i, uint j, uint k, uint n) const
{
    uint id = get_index(i, j, k);
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
#ifdef WITH_OMP
#   pragma omp parallel for
#endif
    real avgDeviation = 0;

//    long long totalDevDiv = 0;
    for (int i = 0; i < ic; ++i) {
        for (uint j = 0; j < jc; ++j) {
            for (uint k = 0; k < kc; ++k) {
                avgDeviation += getDeviation(arr, i, j, k, n);
//                uint index = get_index(i, j, k);
//                avgDeviation += ABS(array[index]-arrayPP[index]);
            }
        }
    }
    real globalDeviation=0;
    MPI_Reduce(&avgDeviation, &globalDeviation, 1, MPI_TYPE_REAL,
               MPI_SUM, 0, MPI_COMM_WORLD);
    if (cnode.mpi.rank == 0) {
        cnode.print(SSTR("global deviation for step " << n << " equals "
                         << avgDeviation / (double(N) * N * N)));
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
    j0 = MIN(n.y, iMissedItemCount) * (jc + 1) + MAX(n.y - iMissedItemCount, 0) * jc;
    k0 = MIN(n.z, iMissedItemCount) * (kc + 1) + MAX(n.z - iMissedItemCount, 0) * kc;

    if (cnode.x < iMissedItemCount)
        ++ic;
    if (cnode.y < jMissedItemCount)
        ++jc;
    if (cnode.z < kMissedItemCount)
        ++kc;

    hx = VAL_LX / ic;
    hy = VAL_LY / jc;
    hz = VAL_LZ / kc;

    MY_ASSERT(ABS(hx - hy) < 0.0001);
    MY_ASSERT(ABS(hx - hz) < 0.0001);


    ht = VAL_T / clargs.K;

    bigsize = (ic + 2) * (jc + 2) * (kc + 2);

    // optimization (allocations may throw std::bad_alloc if no enough memory)
    try {
        array.resize(bigsize);
        arrayP.resize(bigsize);
        arrayPP.resize(bigsize);
        analyticalSolution.resize(bigsize);

        for (int i = 0; i < DIR_SIZE; ++i) {
            ConnectionDirection cdir = toCD(i);
            if (n.is(cdir)) {
                Requests::Info info;
                info.dir = cdir;

                send_requests.append(info, dir_size(cdir));
                recv_requests.append(info, dir_size(cdir));
            }
            else if (cdir != DIR_Y_PERIOD_FIRST &&
                     cdir != DIR_Y_PERIOD_LAST) {
                no_neighbour_edges.push_back(cdir);
            }
        }
            cnode.print(SSTR("(x, y, z): (" << n.x << ',' << n.y << ',' << n.z << ')'
                             << " i0, j0, k0 (ic, jc, kc): "
                             << i0 << ',' << j0 << ',' << k0
                             << " (" << ic  << ',' << jc << ',' << kc << ")"
                             << std::endl));
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
//    cnode.print(SSTR("copy_recv " << i << ',' << j << ',' << k
//                     << '(' << get_p_index(i, j, k) << ')'
//                     << '=' << offset));
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

    switch (type) {
    case MPI_OP_RECV:
        cnode.print(SSTR("RECEIVING FROM " << cnode.neighbor(cdir)));
        break;
    case MPI_OP_SEND:
        cnode.print(SSTR("SENDING TO " << cnode.neighbor(cdir)));
        break;
    }

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

//    cnode.print(SSTR("calculate " << i << ',' << j << ',' << k));

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
