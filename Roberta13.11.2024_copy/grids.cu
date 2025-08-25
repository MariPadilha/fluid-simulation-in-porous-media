#include <cuda_runtime.h>
#include <device_launch_parameters.h>


dim3 grid_1d(int n, int threads){
    int blocks = (n + threads - 1) / threads;
    return dim3(blocks, 1, 1);
}

dim3 grid_2d(int imax, int jmax, int threads_x, int threads_y){
    int blocks_x = (imax + threads_x - 1) / threads_x;
    int blocks_y = (jmax + threads_y - 1) / threads_y;
    return dim3(blocks_x, blocks_y, 1);
}