#include "comum.h"
#define idx i*(jmax+1)+j

static __global__ void calc_1(double *dev_res_c, double *dev_ci, double *dev_c, double *dev_c_tau, double *dev_rc, int jmax, int imax, double dtau, double dt){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	
    if(i <= imax-1 && j <= jmax-1){
        dev_res_c[idx] = ((dev_c[idx]-dev_c_tau[idx]) + dev_rc[idx]*dt) * dtau; 
        dev_ci[idx] = dev_c_tau[idx] + dev_res_c[idx];        
    }
}

static __global__ void calc_2(double *dev_res_c, double *dev_ci, double *dev_c, double *dev_c_tau, double *dev_rc, int jmax, int imax, double dtau, double dt){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	if(i <= imax-1 && j <= jmax-1){    
        dev_res_c[idx] = ((dev_c[idx]-dev_c_tau[idx]) + dev_rc[idx]*dt) * dtau; 
        dev_ci[idx] = 0.75 * dev_c_tau[idx] + 0.25 * (dev_ci[idx] + dev_res_c[idx]);   
        
    }
}

static __global__ void calc_3(double *dev_res_c, double *dev_ci, double *dev_c, double *dev_c_tau, double *dev_c_n_tau, double *dev_rc, int jmax, int imax, double dtau, double dt){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_c[idx] = ((dev_c[idx]-dev_c_tau[idx]) + dev_rc[idx]*dt) * dtau; 
        dev_c_n_tau[idx] = 1.0 / 3.0 * dev_c_tau[idx] + 2.0 / 3.0 * (dev_ci[idx] + dev_res_c[idx]);            
    } 
}

//--- solve_C - concentration ---
void solve_C(double *dev_um_n, double *dev_vm_n, double *dev_c, double *dev_c_n_tau, double *dev_c_tau, double *dev_rc){
    dim3 blockDim(16,16);
    dim3 gridDim((imax-3 + blockDim.x - 1)/blockDim.x, (jmax-3 + blockDim.y - 1)/blockDim.y);
    
    // RALSTON'S METHOD (Second Order Runge-Kutta)
    RESC(dev_um_n, dev_vm_n, dev_c_tau, dev_rc);

    calc_1<<<gridDim, blockDim>>>(dev_res_c, dev_ci, dev_c, dev_c_tau, dev_rc, jmax, imax, dtau, dt);

    bcC(dev_ci);
    RESC(dev_um_n, dev_vm_n, dev_ci, dev_rc);
    
    calc_2<<<gridDim, blockDim>>>(dev_res_c, dev_ci, dev_c, dev_c_tau, dev_rc, jmax, imax, dtau, dt);

    bcC(dev_ci);
    RESC(dev_um_n, dev_vm_n, dev_ci, dev_rc);

    calc_3<<<gridDim, blockDim>>>(dev_res_c, dev_ci, dev_c, dev_c_tau, dev_c_n_tau, dev_rc, jmax, imax, dtau, dt);

    bcC(dev_c_n_tau);
}