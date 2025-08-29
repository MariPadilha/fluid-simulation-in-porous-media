#include "comum.h"

__global__ void calc_vol_u(double *dev_vol_u, double *dev_x, double *dev_ym, int imax, int jmax){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    int n_i = imax - 1;
    int total = jmax * n_i;

    if(idx < total){
        int j = idx / n_i + 1;
        int i = idx % n_i + 2;

        dev_vol_u[j * (imax+2) + i] = (dev_x[i] - dev_x[i-1]) * (dev_ym[j+1] - dev_ym[j]);
    }
}

__global__ void calc_vol_v(double *dev_vol_v, double *dev_y, double *dev_xm, int imax, int jmax){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    int n_i = imax;
    int total = (jmax-1) * n_i;

    if(idx < total){
        int j = idx / n_i + 2;
        int i = idx % n_i + 1;
        dev_vol_v[j * (imax+1) + i] = (dev_y[j]-dev_y[j-1]) * (dev_xm[i+1]-dev_xm[i]); 
    }
}

__global__ void calc_vol_p(double *dev_vol_p, double *dev_ym, double *dev_xm, int imax, int jmax){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    int n_i = imax;
    int total = (jmax-1) * n_i;

    if(idx < total){
        int j = idx / n_i + 2;
        int i = idx % n_i + 1;
        dev_vol_p[j * (imax+1) + i] = (dev_ym[j+1]-dev_ym[j]) * (dev_xm[i+1]-dev_xm[i]); 
    }
}

void calcula_vol_u_v_p(){
    int threads = 256;

    dim3 blocks = grid_1d(jmax * (imax - 1), threads);
    calc_vol_u<<<blocks, threads>>>(dev_vol_u, dev_x, dev_ym, imax, jmax);

    blocks = grid_1d((jmax-1) * imax, threads);
    calc_vol_v<<<blocks, threads>>>(dev_vol_v, dev_y, dev_xm, imax, jmax);
    calc_vol_p<<<blocks, threads>>>(dev_vol_p, dev_ym, dev_xm, imax, jmax);

    cudaDeviceSynchronize();
}