#include "iterations.hpp"


#include <unordered_map>
#include <vector>

#include <cstring>
#include <cmath>

int main(int argc, char *argv[])
{
    // ComputeNode RAII MPI resources (MPI_Initialize, MPI_Finalize)
    ComputeNode cnode; // rank, size

    Iterations its(cnode); // iterations parameters, send/recv buffers

    its.step0(); // u^0
    its.step1(); // u^1

    // STEPS
    for (uint step = 0; step < its.K; ++step) {
        // async receive prev
        {
            const Iterations::Requests &requests = its.recv_requests;
            // asynchronous recv from every neighbor processor
            for (uint i = 0; i < requests.size(); ++i) {
                MPI_Irecv(requests.buff[i].data(),
                          requests.buff[i].size(),
                          MPI_TYPE_REAL,
                          cnode.neighbor(requests.iv[i].dir),
                          TAG_BOUNDARY_ELEMENTS,
                          MPI_COMM_WORLD,
                          &requests.v[i]);
            }
        }
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
                its.copy_data(request_index, &Iterations::Requests::copy_recv)
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
