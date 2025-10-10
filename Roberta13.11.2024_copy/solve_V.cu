#include "comum.h"
#define idx i*(jmax+2)+j

__global__ void calc_1(double *dev_res_v, double *dev_vm, double *dev_vm_tau, double *dev_rv, double *dev_vi, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_v[idx] = ((dev_vm[idx]-dev_vm_tau[idx]) + dev_rv[idx]*dt) * dtau;
        dev_vi[idx] = (dev_vm_tau[idx] + dev_res_v[idx]);
    }
}

__global__ void calc_2(double *dev_res_v, double *dev_vm, double *dev_vm_tau, double *dev_rv, double *dev_vi, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_v[idx] = ((dev_vm[idx]-dev_vm_tau[idx]) + dev_rv[idx]*dt) * dtau;
        dev_vi[idx] = ((double)(3.0/4.0) * dev_vm_tau[idx] + (double)(1.0/4.0) * (dev_vi[idx] + dev_res_v[idx]));
    }
}

__global__ void calc_3(double *dev_res_v, double *dev_vm, double *dev_vm_tau, double *dev_vm_n_tau, double *dev_rv, double *dev_vi, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_v[idx] =((dev_vm[idx]-dev_vm_tau[idx]) +  dev_rv[idx]*dt) * dtau;
        dev_vm_n_tau[idx] = (double)(1.0 / 3.0) * dev_vm_tau[idx] + (double)(2.0 / 3.0) * (dev_vi[idx] + dev_res_v[idx]);
    }
}

double solve_V(double *dev_um, double *dev_vm, double *dev_vm_n, double *dev_um_tau, double *dev_vm_tau, double *dev_vm_n_tau, double *dev_p, double *dev_t){
    double *dev_vi, *dev_rv, *dev_res_v;
    dim3 blockDim(16,16);
    dim3 gridDim((imax-3 + blockDim.x - 1)/blockDim.x, (jmax-4 + blockDim.y - 1)/blockDim.y);

    cudaMalloc((void**)&dev_rv, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_vi, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_res_v, sizeof(double)*(imax+1)*(jmax+2));

    RESV(dev_um_tau, dev_vm_tau, dev_p, dev_t, dev_rv);

    calc_1<<<gridDim, blockDim>>>(dev_res_v, dev_vm, dev_vm_tau, dev_rv, dev_vi, jmax, imax, dt, dtau);

    bcUV(dev_um_tau, dev_vi);
    RESV(dev_um_tau, dev_vi, dev_p, dev_t, dev_rv);

    calc_2<<<gridDim, blockDim>>>(dev_res_v, dev_vm, dev_vm_tau, dev_rv, dev_vi, jmax, imax, dt, dtau);

    bcUV(dev_um_tau, dev_vi);
    RESV(dev_um_tau, dev_vi, dev_p, dev_t, dev_rv);

    calc_3<<<gridDim, blockDim>>>(dev_res_v, dev_vm, dev_vm_tau, dev_vm_n_tau, dev_rv, dev_vi, jmax, imax, dt, dtau);

    for(int i = 1; i <= imax; i++){
        for(int j = 1; j <= jmax+1; j++){
            printf("[%i][%i] vm_n_tau = %lf\n", i, j, dev_vm_n_tau[i*(jmax+2)+j]);
        }
    }

    bcUV(dev_um_tau, dev_vm_n_tau);

    double residual_v = max_reduce(dev_res_v, imax-1, jmax); 

    cudaFree(dev_rv);
    cudaFree(dev_vi);
    cudaFree(dev_res_v);

    return residual_v;
}