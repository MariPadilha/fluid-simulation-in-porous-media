#include "comum.h"

__global__ void xm_borda(double *dev_x, double *dev_xm, int imax){
    dev_xm[1] = -(dev_x[2] + dev_x[1]) * (double)(1.0/2.0);
    dev_xm[imax+1] = (dev_x[imax] + dev_x[imax-1]) * (double)(1.0/2.0) + (dev_x[imax] - dev_x[imax-1]);
}

__global__ void xm_interior(double *dev_x, double *dev_xm, int imax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if(i <= imax)
        dev_xm[i] = (dev_x[i] + dev_x[i-1]) * (double)(1.0/2.0);
}

__global__ void ym_borda(double *dev_y, double *dev_ym, int jmax){
    dev_ym[1] = dev_y[1] - (dev_y[2] - dev_y[1]) * (double)(1.0/2.0);
    dev_ym[jmax+1] = (dev_y[jmax] + dev_y[jmax-1]) * (double)(1.0/2.0) + (dev_y[jmax] - dev_y[jmax-1]);
}

__global__ void ym_interior(double *dev_y, double *dev_ym, int jmax){
    int j = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if(j <= jmax)
        dev_ym[j] = (dev_y[j] + dev_y[j-1]) * 0.5; 
}


void calcula_xm_ym(){
    int threads = 256;
    dim3 blocks = grid_1d(imax-1, threads);

    xm_borda<<<1,1>>>(dev_x, dev_xm, imax);
    xm_interior<<<blocks, threads>>>(dev_x, dev_xm, imax);

    blocks = grid_1d(jmax-1, threads);
    ym_borda<<<1,1>>>(dev_y, dev_ym, jmax);
    ym_interior<<<blocks, threads>>>(dev_y, dev_ym, jmax);
}