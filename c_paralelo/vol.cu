#include "comum.h"

__global__ void calc_vol_u(double *dev_x, double *dev_ym, double *dev_vol_u, int imax, int jmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if(i <= imax && j <= jmax){
        dev_vol_u[i*(jmax+1)+j] = (dev_x[i]-dev_x[i-1]) * (dev_ym[j+1]-dev_ym[j]); 
    }
}

__global__ void calc_vol_v(double *dev_xm, double *dev_y, double *dev_vol_v, int imax, int jmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;

    if(i <= imax && j <= jmax)
        dev_vol_v[i*(jmax+2)+j] = (dev_y[j]-dev_y[j-1]) * (dev_xm[i+1]-dev_xm[i]); 
}

__global__ void calc_vol_p(double *dev_xm, double *dev_ym, double *dev_vol_p, int imax, int jmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;

    if(i <= imax && j <= jmax)
            dev_vol_p[i*(jmax+1)+j] = (dev_ym[j+1]-dev_ym[j]) * (dev_xm[i+1]-dev_xm[i]); 
}

void calcula_vol(){
    dim3 blockDim(16,16);
    dim3 gridDimVolVP((imax + blockDim.x - 1)/blockDim.x, (jmax-1 + blockDim.y - 1)/blockDim.y);
    dim3 gridDimVolU((imax-1 + blockDim.x - 1)/blockDim.x, (jmax + blockDim.y - 1)/blockDim.y);
    calc_vol_u<<<gridDimVolVP, blockDim>>>(dev_x, dev_ym, dev_vol_u, imax, jmax);
    calc_vol_p<<<gridDimVolVP, blockDim>>>(dev_xm, dev_ym, dev_vol_p, imax, jmax);
    calc_vol_v<<<gridDimVolU, blockDim>>>(dev_xm, dev_y, dev_vol_v, imax, jmax);
}