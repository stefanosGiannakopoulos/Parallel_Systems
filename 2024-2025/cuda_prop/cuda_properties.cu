#include <stdio.h>
#include <cuda_runtime.h>

int main() {
    int deviceCount;
    cudaError_t err;

    // Get the number of devices
    err = cudaGetDeviceCount(&deviceCount);
    if (err != cudaSuccess) {
        printf("Error getting device count: %s\n", cudaGetErrorString(err));
        return -1;
    }

    printf("Number of CUDA devices: %d\n", deviceCount);

    // Loop through each device
    for (int device = 0; device < deviceCount; ++device) {
        cudaDeviceProp deviceProp;

        // Get device properties
        err = cudaGetDeviceProperties(&deviceProp, device);
        if (err != cudaSuccess) {
            printf("Error getting properties for device %d: %s\n", device, cudaGetErrorString(err));
            continue;
        }

        printf("\nDevice %d: \"%s\"\n", device, deviceProp.name);
        printf("  Compute capability: %d.%d\n", deviceProp.major, deviceProp.minor);
        printf("  Total global memory: %lu bytes\n", deviceProp.totalGlobalMem);
        printf("  Shared memory per block: %lu bytes\n", deviceProp.sharedMemPerBlock);
        printf("  Registers per block: %d\n", deviceProp.regsPerBlock);
        printf("  Warp size: %d\n", deviceProp.warpSize);
        printf("  Max threads per block: %d\n", deviceProp.maxThreadsPerBlock);
        printf("  Max threads dimensions: [%d, %d, %d]\n",
               deviceProp.maxThreadsDim[0],
               deviceProp.maxThreadsDim[1],
               deviceProp.maxThreadsDim[2]);
        printf("  Max grid size: [%d, %d, %d]\n",
               deviceProp.maxGridSize[0],
               deviceProp.maxGridSize[1],
               deviceProp.maxGridSize[2]);
        printf("  Clock rate: %.2f MHz\n", deviceProp.clockRate / 1000.0);
        printf("  Memory clock rate: %.2f MHz\n", deviceProp.memoryClockRate / 1000.0);
        printf("  Memory bus width: %d bits\n", deviceProp.memoryBusWidth);
        printf("  Multiprocessor count: %d\n", deviceProp.multiProcessorCount);
        printf("  L2 cache size: %d bytes\n", deviceProp.l2CacheSize);
        printf("  Concurrent kernels: %s\n", deviceProp.concurrentKernels ? "Yes" : "No");
        printf("  ECC enabled: %s\n", deviceProp.ECCEnabled ? "Yes" : "No");
    }

    return 0;
}
