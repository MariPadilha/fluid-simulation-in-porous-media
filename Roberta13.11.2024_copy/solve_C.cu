#include "comum.h"

__global__ void calc_1(double *dev_dudx, double *dev_dvdy, double *dev_rp, double *dev_pi, double *dev_um_n, double *dev_vm_n, double *dev_areau_e, double *dev_areav_n, double *dev_areau_w, double *dev_areav_s, double *dev_p, int imax, int jmax, double dtau, double beta){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	
    if(i <= imax-1 && j <= jmax-1){
        dev_dudx[i*(jmax+1)+j] = dev_um_n[(i+1)*(jmax+1)+j] * dev_areau_e[j] - dev_um_n[i*(jmax+1)+j] * dev_areau_w[j];
        dev_dvdy[i*(jmax+1)+j] = dev_vm_n[i*(jmax+1)+(j+1)] * dev_areav_n[i] - dev_vm_n[i*(jmax+1)+j] * dev_areav_s[i];
        dev_rp[i*(jmax+1)+j] = - (dev_dudx[i*(jmax+1)+j]+ dev_dvdy[i*(jmax+1)+j]); 
        dev_pi[i*(jmax+1)+j] = dev_p[i*(jmax+1)+j] + dtau * dev_rp[i*(jmax+1)+j] * beta;
    }
}

__global__ void calc_2(double *dev_pi, double *dev_p, double *dev_rp, int jmax, int imax, double dtau, double beta){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
	if(i <= imax-1 && j <= jmax-1){
        dev_pi[i*(jmax+1)+j] = (double)(3.0/4.0) * dev_p[i*(jmax+1)+j] + (double)(1.0/4.0) * (dev_pi[i*(jmax+1)+j] + dtau * dev_rp[i*(jmax+1)+j] * beta);
    }
}

__global__ void calc_3(double *dev_res_p, double *dev_pn, double *dev_rp, double *dev_p, int jmax, int imax, double dtau, double beta){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_p[i*(jmax+1)+j] = dtau * dev_rp[i*(jmax+1)+j] * beta;
        dev_pn[i*(jmax+1)+j] = (double)(1.0 / 3.0) * dev_p[i*(jmax+1)+j] + (double)(2.0 / 3.0) * (dev_pi[i*(jmax+1)+j] + dev_res_p[i*(jmax+1)+j]);
    }
}

void solve_P(double *dev_p, double *dev_um_n, double *dev_vm_n, double *dev_pn, double *residual_p){
    int i, j;
    dim3 threads(16,16), blocks;

    blocks = grid_2d((imax-1-2), (jmax-1-2), threads);
    calc_1<<<blocks, threads>>>(dev_dudx, dev_dvdy, dev_rp, dev_pi, dev_um_n, dev_vm_n, dev_areau_e, dev_areav_n, dev_areau_w, dev_areav_s, dev_p, imax, jmax, dtau, iterations.beta);

    bcP(dev_pi);

    blocks = grid_2d((imax-1-2), (jmax-1-3), threads);
    calc_2<<<blocks, threads>>>(dev_pi, dev_p, dev_rp, jmax, imax, dtau, iterations.beta);

    bcP(dev_pi);

    blocks = grid_2d((imax-1-2), (jmax-1-2), threads);
    calc_3<<<blocks, threads>>>(dev_res_p, dev_pn, dev_rp, dev_p, jmax, imax, dtau, iterations.beta);

    bcP(dev_pn);

    (*residual_p) = maior_valor(dev_res_p, imax+1, jmax+1, 1, 1);
}



//--- solve_C - concentration ---
void solve_C(double **um_n, double **vm_n, double **c, double **c_n_tau, double **c_tau){
    int i, j;
    double *dev_rc;

    cudaMalloc((void**)&dev_rc, sizeof(double)*(imax+1)*(jmax+1));

    // RALSTON'S METHOD (Second Order Runge-Kutta)
    RESC(um_n, vm_n, c_tau, dev_rc);

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_c[i][j] = ((c[i][j]-c_tau[i][j]) + dev_rc[i*(jmax+1)+j]*dt) * dtau; 
            ci[i][j] = c_tau[i][j] + res_c[i][j];
        }
    }

    bcC(ci);
    RESC(um_n,vm_n,ci,dev_rc);
    
    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_c[i][j] = ((c[i][j]-c_tau[i][j]) + dev_rc[i*(jmax+1)+j]*dt) * dtau; 
            ci[i][j] = 0.75 * c_tau[i][j] + 0.25 * (ci[i][j] + res_c[i][j]);   
        }
    }

    bcC(ci);
    RESC(um_n, vm_n, ci, dev_rc);

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_c[i][j] = ((c[i][j]-c_tau[i][j]) + dev_rc[i*(jmax+1)+j]*dt) * dtau; 
            c_n_tau[i][j] = 1.0 / 3.0 * c_tau[i][j] + 2.0 / 3.0 * (ci[i][j] + res_c[i][j]);            
        }
    }

    bcC(c_n_tau);

    cudaFree(dev_rc);
}