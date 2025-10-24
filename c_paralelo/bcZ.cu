#include "comum.h"

static __global__ void contorno_superior_inferior(double *dev_zt, int imax, int jmax, double tinf){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if(i >= imax) return;
    
    dev_zt[i*(jmax+1)+1] = tinf;
    dev_zt[i*(jmax+1)+jmax] = dev_zt[i*(jmax+1)+(jmax-1)] + (dev_zt[i*(jmax+1)+(jmax-1)]-dev_zt[i*(jmax+1)+(jmax-2)]);
}

static __global__ void contorno_esquerdo_direito(double *dev_zt, int imax, int jmax){
    int j = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if(j >= jmax) return;
    
    dev_zt[1*(jmax+1)+j] = dev_zt[2*(jmax+1)+j];
    dev_zt[imax*(jmax+1)+j] = dev_zt[(imax-1)*(jmax+1)+j];
}

static __global__ void contorno_cantos(double *dev_zt, int imax, int jmax){
    if(threadIdx.x == 0 && blockIdx.x == 0){
        dev_zt[1*(jmax+1)+1] = dev_zt[2*(jmax+1)+1];
        dev_zt[1*(jmax+1)+jmax] = dev_zt[2*(jmax+1)+jmax];
        dev_zt[imax*(jmax+1)+1] = dev_zt[(imax-1)*(jmax+1)+1];
        dev_zt[imax*(jmax+1)+jmax] = dev_zt[(imax-1)*(jmax+1)+jmax];
    }
}

void bcZ(double *dev_zt){
    int threads = 256;
    int blocks_h = grid_1d(imax-2, threads);
    int blocks_v = grid_1d(jmax-2, threads);
    

    cudaStream_t s1, s2;
    cudaStreamCreate(&s1);
    cudaStreamCreate(&s2);

    contorno_superior_inferior<<<blocks_h, threads, 0, s1>>>(dev_zt, imax, jmax, tinf);
    contorno_esquerdo_direito<<<blocks_v, threads, 0, s2>>>(dev_zt, imax, jmax);

    cudaStreamSynchronize(s1);
    cudaStreamSynchronize(s2);
    
    contorno_cantos<<<1, 1>>>(dev_zt, imax, jmax);
    cudaDeviceSynchronize();

    cudaStreamDestroy(s1);
    cudaStreamDestroy(s2);
}