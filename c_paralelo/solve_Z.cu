#include "comum.h"
#define idx i*(jmax+1)+j

static __global__ void calc_1(double *dev_res_z, double *dev_zi, double *dev_z, double *dev_z_tau, double *dev_rz, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	
    if(i <= imax-1 && j <= jmax-1){
        dev_res_z[idx] = ((dev_z[idx]-dev_z_tau[idx]) + dev_rz[idx]*dt) * dtau; 
        dev_zi[idx] = dev_z_tau[idx] + dev_res_z[idx];             
    }
}

static __global__ void calc_2(double *dev_res_z, double *dev_zi, double *dev_z, double *dev_z_tau, double *dev_rz, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_z[idx] = ((dev_z[idx]-dev_z_tau[idx]) + dev_rz[idx]*dt) * dtau;
        dev_zi[idx] = 0.75 * dev_z_tau[idx] + 0.25 * (dev_zi[idx]+dev_res_z[idx]);
    }
}

static __global__ void calc_3(double *dev_res_z, double *dev_zi, double *dev_z, double *dev_z_tau, double *dev_z_n_tau, double *dev_rz, int jmax, int imax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_z[idx] = ((dev_z[idx]-dev_z_tau[idx]) + dev_rz[idx]*dt) * dtau; 
        dev_z_n_tau[idx] = 1.0 / 3.0 * dev_z_tau[idx] + 2.0 / 3.0 * (dev_zi[idx]+dev_res_z[idx]);
    }
}

//--- Solve Mixture Fraction ---
void solve_Z(double *dev_um_n, double *dev_vm_n, double *dev_z, double *dev_z_n_tau, double *dev_z_tau, double *dev_rz){
    dim3 blockDim(16,16);
    dim3 gridDim((imax-3 + blockDim.x - 1)/blockDim.x, (jmax-3 + blockDim.y - 1)/blockDim.y);


    // RALSTON'S METHOD (Second Order Runge-Kutta)
    RESZ(dev_um_n, dev_vm_n, dev_z_tau, dev_rz);

    calc_1<<<gridDim, blockDim>>>(dev_res_z, dev_zi, dev_z, dev_z_tau, dev_rz, jmax, imax, dt, dtau);

    bcZ(dev_zi);
    RESZ(dev_um_n, dev_vm_n, dev_zi, dev_rz);

    calc_2<<<gridDim, blockDim>>>(dev_res_z, dev_zi, dev_z, dev_z_tau, dev_rz, jmax, imax, dt, dtau);
    
    bcZ(dev_zi);
    RESZ(dev_um_n, dev_vm_n, dev_zi, dev_rz);

    calc_3<<<gridDim, blockDim>>>(dev_res_z, dev_zi, dev_z, dev_z_tau, dev_z_n_tau, dev_rz, jmax, imax, dt, dtau);

    bcZ(dev_z_n_tau);
}