#include "comum.h"

__global__ void calc_1(double *dev_res_z, double *dev_zi){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	
    if(i <= imax-1 && j <= jmax-1){
        res_z[i][j] = ((z[i][j]-z_tau[i][j]) + dev_rz[i*(jmax+1)+j]*dt) * dtau; 
        zi[i][j] = z_tau[i][j] + res_z[i][j];             
    }
}

__global__ void calc_2(double *dev_pi, double *dev_p, double *dev_rp, int jmax, int imax, double dtau, double beta){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
	if(i <= imax-1 && j <= jmax-1){
        for(i = 2; i <= imax-1; i++){    
            for(j = 2; j <= jmax-1; j++){
            res_z[i][j] = ((z[i][j]-z_tau[i][j]) + dev_rz[i*(jmax+1)+j]*dt) * dtau;
            zi[i][j] = 0.75 * z_tau[i][j] + 0.25 * (zi[i][j]+res_z[i][j]);
        }
    }
    }
}

__global__ void calc_3(double *dev_res_p, double *dev_pn, double *dev_rp, double *dev_p, int jmax, int imax, double dtau, double beta){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	if(i <= imax-1 && j <= jmax-1){
        for(i = 2; i <= imax-1; i++){
            for(j = 2; j <= jmax-1; j++){
            res_z[i][j] = ((z[i][j]-z_tau[i][j]) + dev_rz[i*(jmax+1)+j]*dt) * dtau; 
            z_n_tau[i][j] = 1.0 / 3.0 * z_tau[i][j] + 2.0 / 3.0 * (zi[i][j]+res_z[i][j]);
        }
    }
    }
}

//--- Solve Mixture Fraction ---
void solve_Z(double *dev_um_n, double *dev_vm_n, double **z, double **z_n_tau, double **z_tau){
    int i, j;
    double *dev_rz;
    dim3 threads(16,16), blocks;

    cudaMalloc((void**)&dev_rz, sizeof(double)*(imax+1)*(jmax+1));

    // RALSTON'S METHOD (Second Order Runge-Kutta)
    RESZ(dev_um_n, dev_vm_n, z_tau, dev_rz);

    blocks = grid_2d((imax-1-2), (jmax-1-2), threads);
    calc_1<<<blocks, threads>>>(dev_dudx, dev_dvdy, dev_rp, dev_pi, dev_um_n, dev_vm_n, dev_areau_e, dev_areav_n, dev_areau_w, dev_areav_s, dev_p, imax, jmax, dtau, iterations.beta);

    bcZ(zi);
    RESZ(um_n,vm_n,zi,dev_rz);

    blocks = grid_2d((imax-1-2), (jmax-1-2), threads);
    calc_2<<<blocks, threads>>>(dev_dudx, dev_dvdy, dev_rp, dev_pi, dev_um_n, dev_vm_n, dev_areau_e, dev_areav_n, dev_areau_w, dev_areav_s, dev_p, imax, jmax, dtau, iterations.beta);
    
    bcZ(zi);
    RESZ(um_n,vm_n,zi,dev_rz);

    blocks = grid_2d((imax-1-2), (jmax-1-2), threads);
    calc_3<<<blocks, threads>>>(dev_dudx, dev_dvdy, dev_rp, dev_pi, dev_um_n, dev_vm_n, dev_areau_e, dev_areav_n, dev_areau_w, dev_areav_s, dev_p, imax, jmax, dtau, iterations.beta);

    bcZ(z_n_tau);

    cudaFree(dev_rz);
}