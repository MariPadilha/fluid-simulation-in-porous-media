#include "comum.h"

__global__ void calc_dx(double *dev_dx, double *dev_x, int imax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if(i > imax) return;
    
    dev_dx[i] = dev_x[i]-dev_x[i-1];
}

__global__ void calc_dy(double *dev_dy, double *dev_y, int jmax){
    int j = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if(j > jmax) return;
        
    dev_dy[j] = dev_y[j]-dev_y[j-1];
}

void calcula_dx_dy(){
    int threads = 256;
    dim3 blocks = grid_1d(imax, threads);
    calc_dx<<<blocks, threads>>>(dev_dx, dev_x, imax);

    blocks = grid_1d(jmax, threads);
    calc_dy<<<blocks, threads>>>(dev_dy, dev_y, jmax);
}