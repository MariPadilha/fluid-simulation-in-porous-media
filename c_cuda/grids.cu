#include <cuda_runtime.h>
#include <device_launch_parameters.h>

int grid_1d(int tamanho, int threads){
    return ((tamanho) + threads - 1) / threads;
}