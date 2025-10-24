#include "comum.h"

static __global__ void contorno_superior_inferior(double *dev_cn, int imax, int jmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if(i >= imax) return;
    
    dev_cn[i*(jmax+1)+1] = 0.0;
    dev_cn[i*(jmax+1)+jmax] = dev_cn[i*(jmax+1)+(jmax-1)] + (dev_cn[i*(jmax+1)+(jmax-1)]-dev_cn[i*(jmax+1)+(jmax-2)]);
}

static __global__ void contorno_esquerdo_direito(double *dev_cn, int imax, int jmax){
    int j = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if(j >= jmax) return;
    
    dev_cn[1*(jmax+1)+j] = dev_cn[2*(jmax+1)+j];
    dev_cn[imax*(jmax+1)+j] = dev_cn[(imax-1)*(jmax+1)+j];
}

static __global__ void contorno_cantos(double *dev_cn, int imax, int jmax){
    if(threadIdx.x == 0 && blockIdx.x == 0){
        dev_cn[1*(jmax+1)+1] = dev_cn[2*(jmax+1)+1];
        dev_cn[1*(jmax+1)+jmax] = dev_cn[2*(jmax+1)+jmax];
        dev_cn[imax*(jmax+1)+1] = dev_cn[imax-1*(jmax+1)+1];
        dev_cn[imax*(jmax+1)+jmax] = dev_cn[imax-1*(jmax+1)+jmax];
    }
}

void bcC(double *dev_cn){
    int threads = 256;
    int blocks_h = grid_1d(imax-2, threads);
    int blocks_v = grid_1d(jmax-2, threads);
    
    cudaStream_t s1, s2;
    cudaStreamCreate(&s1);
    cudaStreamCreate(&s2);

    contorno_superior_inferior<<<blocks_h, threads, 0, s1>>>(dev_cn, imax, jmax);
    contorno_esquerdo_direito<<<blocks_v, threads, 0, s2>>>(dev_cn, imax, jmax);

    cudaStreamSynchronize(s1);
    cudaStreamSynchronize(s2);
    
    contorno_cantos<<<1, 1>>>(dev_cn, imax, jmax);
    cudaDeviceSynchronize();

    cudaStreamDestroy(s1);
    cudaStreamDestroy(s2);
}