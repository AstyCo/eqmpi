#include "iterations.hpp"

#ifdef CUDA
#include "cuda_runtime.h"
#endif

#ifdef SC_INFO
#include <cstdio>
#include <unistd.h> // sysconf
int main(int argc, char **argv)
{
    int num_cpu = sysconf(_SC_NPROCESSORS_ONLN);
    int dev_count;
    cudaGetDeviceCount(&dev_count);

    std::cout << "num_cpu " << num_cpu << std::endl
              << "dev_count " << dev_count << std::endl;
    cudaGetDeviceCount(&dev_count);
    std::cout << "float size" << sizeof(float) << std::endl
              << "int size " << sizeof(int) << std::endl
              << "long size " << sizeof(long) << std::endl;

    for (int i = 0; i < dev_count; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        printf("Device Number: %d\n", i);
        printf("  Device name: %s\n", prop.name);
        printf("  Memory Clock Rate (KHz): %d\n",
               prop.memoryClockRate);
        printf("  Memory Bus Width (bits): %d\n",
               prop.memoryBusWidth);
        printf("  Peak Memory Bandwidth (GB/s): %f\n",
               2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
        printf("  sharedMemPerBlock (bytes): %lu\n",
               prop.sharedMemPerBlock);
        printf("  totalGlobalMem (bytes): %lu\n",
               prop.totalGlobalMem);
        printf("  l2CacheSize (bytes): %d\n",
               prop.l2CacheSize);
        printf("  maxThreadsPerBlock: %d\n",
               prop.maxThreadsPerBlock);
        printf("  maxThreadsDim: %d %d %d\n",
               prop.maxThreadsDim[0], prop.maxThreadsDim[1],
               prop.maxThreadsDim[2]);
        printf("  sharedMemPerMultiprocessor (bytes): %lu\n",
               prop.sharedMemPerMultiprocessor);
        printf("  regsPerMultiprocessor (32-bit): %d\n",
               prop.regsPerMultiprocessor);
        printf("  maxGridSize: %d %d %d\n",
               prop.maxGridSize[0], prop.maxGridSize[1],
               prop.maxGridSize[2]);


        printf("\n");
    }
    return 0;
}

#else
int main(int argc, char **argv)
{
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    clargs.parse(argc, argv);

    cnode.init(); // rank, size

    std::vector<int> Ns;
    if (clargs.N > 0)
        Ns.push_back(clargs.N);
    if (Ns.empty()) {
        Ns.push_back(128);
        Ns.push_back(256);
        Ns.push_back(512);
    }

    for (uint i = 0; i < Ns.size(); ++i) {
        times.clear();
        int N = Ns[i];
        MPI_Barrier(MPI_COMM_WORLD);
        Profiler profiler;
        Profiler p_finalization;
        {

            Profiler p_init;
            Iterations its(N); // iterations parameters, send/recv buffers
            get_time(times.program_initialization, p_init);

            {
                Profiler p_s0_s1;
                its.prepare();
                get_time(times.parallel_cycles, p_s0_s1);
            }

            its.run();

            MPI_Barrier(MPI_COMM_WORLD);
            p_finalization.start();
        }
        get_time(times.total, profiler);
        get_time(times.program_finalization, p_finalization);

        times.reduce();
        if (cnode.mpi.rank == 0)
            std::cout << times.get_times(N) << std::endl;
    }


    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
#endif
