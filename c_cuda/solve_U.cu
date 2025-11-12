#include "comum.h"
#define idx i*(jmax+1)+j

static __global__ void calc_1(double *dev_res_u, double *dev_um, double *dev_um_tau, double *dev_ru, double *dev_ui, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 3;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_u[idx] = ((dev_um[idx]-dev_um_tau[idx]) + dev_ru[idx]*dt) * dtau;
        dev_ui[idx] = dev_um_tau[idx] + dev_res_u[idx];
    }
}

static __global__ void calc_2(double *dev_res_u, double *dev_um, double *dev_um_tau, double *dev_ru, double *dev_ui, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 3;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_u[idx] = ((dev_um[idx]-dev_um_tau[idx]) + dev_ru[idx]*dt) * dtau;
        dev_ui[idx] = 0.75 * dev_um_tau[idx] + 0.25 * (dev_ui[idx]+dev_res_u[idx]);            
    }
}

static __global__ void calc_3(double *dev_res_u, double *dev_um, double *dev_um_tau, double *dev_um_n_tau, double *dev_ru, double *dev_ui, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 3;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_u[idx] = ((dev_um[idx]-dev_um_tau[idx]) + dev_ru[idx]*dt) * dtau;
        dev_um_n_tau[idx] = 1.0 / 3.0 * dev_um_tau[idx] + 2.0 / 3.0 * (dev_ui[idx]+dev_res_u[idx]); 
    }
}

double solve_U(double *dev_um, double *dev_vm, double *dev_um_n, double *dev_um_tau, double *dev_vm_tau, double *dev_um_n_tau, double *dev_p, double *dev_ui, double *dev_ru, double *dev_res_u){
    dim3 blockDim(16, 16);
    dim3 gridDim((imax-4 + blockDim.x - 1)/blockDim.x, (jmax-3 + blockDim.y - 1)/blockDim.y);
    
    
    RESU(dev_um_tau, dev_vm_tau, dev_p, dev_ru);

    calc_1<<<gridDim, blockDim>>>(dev_res_u, dev_um, dev_um_tau, dev_ru, dev_ui, jmax, imax, dt, dtau);
    
    bcUV(dev_ui, dev_vm_tau);

    
    RESU(dev_ui, dev_vm_tau, dev_p, dev_ru);
    
    calc_2<<<gridDim, blockDim>>>(dev_res_u, dev_um, dev_um_tau, dev_ru, dev_ui, jmax, imax, dt, dtau);
    
    bcUV(dev_ui, dev_vm_tau);
    
    
    RESU(dev_ui, dev_vm_tau, dev_p, dev_ru);
    
    
    calc_3<<<gridDim, blockDim>>>(dev_res_u, dev_um, dev_um_tau, dev_um_n_tau, dev_ru, dev_ui, jmax, imax, dt, dtau);
    
    
    bcUV(dev_um_n_tau, dev_vm_tau);
    
    double residual_u =  max_reduce(dev_res_u, imax-1, jmax-2); 
    
    return residual_u;
}