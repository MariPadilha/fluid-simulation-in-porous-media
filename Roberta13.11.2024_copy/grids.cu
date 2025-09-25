#include <cuda_runtime.h>
#include <device_launch_parameters.h>

dim3 grid_2d(int imax, int jmax, dim3 threads){
    return dim3(
        ((imax+1) + threads.x - 1) / threads.x,
        ((jmax+1) + threads.y - 1) / threads.y
    );
}

int grid_1d(int tamanho, int threads){
    return ((tamanho+1) + threads - 1) / threads;
}