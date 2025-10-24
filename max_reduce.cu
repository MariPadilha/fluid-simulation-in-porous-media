#include <stdio.h>
#include <cuda_runtime.h>
#include <cub/cub.cuh>

double max_reduce(double* dev_matriz_linearizada, int imax, int jmax){
    int tamanho = imax * jmax;
    double *dev_max;
    cudaMalloc(&dev_max, sizeof(double));
    void *d_temp_storage = NULL;
    size_t temp_storage_bytes = 0;

    
    cub::DeviceReduce::Max(d_temp_storage, temp_storage_bytes, dev_matriz_linearizada, dev_max, tamanho);
    cudaMalloc(&d_temp_storage, temp_storage_bytes);

    cub::DeviceReduce::Max(d_temp_storage, temp_storage_bytes, dev_matriz_linearizada, dev_max, tamanho);

    double h_max;
    cudaMemcpy(&h_max, dev_max, sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(dev_max);
    cudaFree(d_temp_storage);

    return h_max;
}
