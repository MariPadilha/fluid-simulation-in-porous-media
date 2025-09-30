#include "comum.h"

__global__ void pontos_medios_velocidade(double *dev_u, double *dev_v, double *dev_um, double *dev_vm, int imax, int jmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if(i <= imax && j <= jmax){
        dev_u[i*(jmax+1)+j] = (dev_um[(i+1)*(jmax+1)+j]+dev_um[i*(jmax+1)+j])*0.50;
        dev_v[i*(jmax+1)+j] = (dev_vm[i*(jmax+2)+(j+1)]+dev_vm[i*(jmax+2)+j])*0.50;
    }
}

void comp_mean(double *dev_u, double *dev_v, double *dev_um, double *dev_vm){
    dim3 threads(16,16);
    dim3 blocks = grid_2d(imax-1, jmax-1, threads);
    
    pontos_medios_velocidade<<<blocks, threads>>>(dev_u, dev_v, dev_um, dev_vm, imax, jmax);
}