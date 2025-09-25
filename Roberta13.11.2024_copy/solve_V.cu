#include "comum.h"

__global__ void calc_1(double *dev_res_v, double *dev_vm, double *dev_vm_tau, double *dev_rv, double *dev_vi, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_v[i*(jmax+2)+j] = ((dev_vm[i*(jmax+2)+j]-dev_vm_tau[i*(jmax+2)+j]) + dev_rv[i*(jmax+2)+j]*dt) * dtau;
        dev_vi[i*(jmax+2)+j] = (dev_vm_tau[i*(jmax+2)+j] + dev_res_v[i*(jmax+2)+j]);
    }
}

__global__ void calc_2(double *dev_res_v, double *dev_vm, double *dev_vm_tau, double *dev_rv, double *dev_vi, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_v[i*(jmax+2)+j] = ((dev_vm[i*(jmax+2)+j]-dev_vm_tau[i*(jmax+2)+j]) + dev_rv[i*(jmax+2)+j]*dt) * dtau;
        dev_vi[i*(jmax+2)+j] = ((double)(3.0/4.0) * dev_vm_tau[i*(jmax+2)+j] + (double)(1.0/4.0) * (dev_vi[i*(jmax+2)+j] + dev_res_v[i*(jmax+2)+j]));
    }
}

__global__ void calc_3(double *dev_res_v, double *dev_vm, double *dev_vm_tau, double *dev_vm_n_tau, double *dev_rv, double *dev_vi, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_v[i*(jmax+2)+j] =((dev_vm[i*(jmax+2)+j]-dev_vm_tau[i*(jmax+2)+j]) +  dev_rv[i*(jmax+2)+j]*dt) * dtau;
        dev_vm_n_tau[i*(jmax+2)+j] = (double)(1.0 / 3.0) * dev_vm_tau[i*(jmax+2)+j] + (double)(2.0 / 3.0) * (dev_vi[i*(jmax+2)+j] + dev_res_v[i*(jmax+2)+j]);
    }
}

void solve_V(double *dev_um, double *dev_vm, double *dev_vm_n, double *dev_um_tau, double *dev_vm_tau, double *dev_vm_n_tau, double **p, double **T, double *residual_v){
    int i, j;
    double *dev_vi, *dev_rv, *dev_res_v;
    dim3 threads(16,16), blocks;

    cudaMalloc((void**)&dev_rv, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_vi, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_res_v, sizeof(double)*(imax+1)*(jmax+2));

    RESV(dev_um_tau, dev_vm_tau, p, T, dev_rv);

    blocks = grid_2d((imax-1-2), (jmax-1-3), threads);
    calc_1<<<blocks, threads>>>(dev_res_v, dev_vm, dev_vm_tau, dev_rv, dev_vi, jmax, imax, dt, dtau);

    bcUV(dev_um_tau, dev_vi);
    RESV(dev_um_tau, dev_vi, p, T, dev_rv);

    calc_2<<<blocks, threads>>>(dev_res_v, dev_vm, dev_vm_tau, dev_rv, dev_vi, jmax, imax, dt, dtau);

    bcUV(dev_um_tau, dev_vi);
    RESV(dev_um_tau, dev_vi, p, T, dev_rv);

    calc_3<<<blocks, threads>>>(dev_res_v, dev_vm, dev_vm_tau, dev_vm_n_tau, dev_rv, dev_vi, jmax, imax, dt, dtau);

    bcUV(dev_um_tau, dev_vm_n_tau);
    (*residual_v) = maior_valor(res_v, imax+1, jmax+2, 2, 3); 

    //desalocando
    cudaFree(dev_rv);
    cudaFree(dev_vi);
    cudaFree(dev_res_v);
}