#include "comum.h"

__global__ void calc_1(double *dev_res_u, double *dev_um, double *dev_um_tau, double *dev_ru, double *dev_ui, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 3;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_u[i*(jmax+1)+j] = ((dev_um[i*(jmax+1)+j]-dev_um_tau[i*(jmax+1)+j]) + dev_ru[i*(jmax+1)+j]*dt) * dtau;
        dev_ui[i*(jmax+1)+j] = dev_um_tau[i*(jmax+1)+j] + dev_res_u[i*(jmax+1)+j];
    }
}

__global__ void calc_2(double *dev_res_u, double *dev_um, double *dev_um_tau, double *dev_ru, double *dev_ui, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 3;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_u[i][j] = ((dev_um[i][j]-dev_um_tau[i][j]) + dev_ru[i*(jmax+1)+j]*dt) * dtau;
        dev_ui[i][j] = 0.75 * dev_um_tau[i][j] + 0.25 * (dev_ui[i][j]+dev_res_u[i][j]);            
    }
}

__global__ void calc_3(double *dev_res_u, double *dev_um, double *dev_um_tau, double *dev_um_n_tau, double *dev_ru, double *dev_ui, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 3;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_u[i][j] = ((dev_um[i][j]-dev_um_tau[i][j]) + dev_ru[i*(jmax+1)+j]*dt) * dtau;
        dev_um_n_tau[i][j] = 1.0 / 3.0 * dev_um_tau[i][j] + 2.0 / 3.0 * (dev_ui[i][j]+dev_res_u[i][j]); 
    }
}

void solve_U(double *dev_um, double *dev_vm, double *dev_um_n, double *dev_um_tau, double *dev_vm_tau, double *dev_um_n_tau, double **p, double *residual_u){
    int i, j;
    double *dev_ui, *dev_ru, *dev_res_u;
    dim3 blocks, threads(16, 16);

    cudaMalloc((void**)&dev_ru, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_ui, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_res_u, sizeof(double)*(imax+2)*(jmax+1));
    
    RESU(dev_um_tau, dev_vm_tau,p,dev_ru);

    blocks = grid_2d((imax-1-3), (jmax-1-2), threads);
    calc_1<<<blocks, threads>>>(dev_res_u, dev_um, dev_um_tau, dev_ru, dev_ui, jmax, imax, dt, dtau);

    bcUV(dev_ui, dev_vm_tau);
    RESU(dev_ui, dev_vm_tau, p, dev_ru);

    blocks = grid_2d((imax-1-3), (jmax-1-2), threads);
    calc_2<<<blocks, threads>>>(dev_res_u, dev_um, dev_um_tau, dev_ru, dev_ui, jmax, imax, dt, dtau);
    
    bcUV(dev_ui, dev_vm_tau);
    RESU(dev_ui, dev_vm_tau, p, dev_ru);

    calc_3(dev_res_u, dev_um, dev_um_tau, dev_um_n_tau, dev_ru, dev_ui, jmax, imax, dt, dtau);

    bcUV(dev_um_n_tau, dev_vm_tau);

    (*residual_u) =  maior_valor(res_u, imax+2, jmax+1, 3, 2); 

    //desalocando
    cudaFree(dev_ru);
    cudaFree(dev_ui);
    cudaFree(dev_res_u);
}