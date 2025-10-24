#include "comum.h"
#define idx i*(jmax+1)+j

__global__ void calc_1(double *dev_dudx, double *dev_dvdy, double *dev_rp, double *dev_pi, double *dev_um_n, double *dev_vm_n, double *dev_areau_e, double *dev_areav_n, double *dev_areau_w, double *dev_areav_s, double *dev_p, int imax, int jmax, double dtau, double beta){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	
    if(i <= imax-1 && j <= jmax-1){
        dev_dudx[idx] = dev_um_n[(i+1)*(jmax+1)+j] * dev_areau_e[j] - dev_um_n[idx] * dev_areau_w[j];
        dev_dvdy[idx] = dev_vm_n[i*(jmax+2)+(j+1)] * dev_areav_n[i] - dev_vm_n[i*(jmax+2)+j] * dev_areav_s[i];
        dev_rp[idx] = - (dev_dudx[idx]+ dev_dvdy[idx]); 
        dev_pi[idx] = dev_p[idx] + dtau * dev_rp[idx] * beta;
    }
}

__global__ void calc_2(double *dev_pi, double *dev_p, double *dev_rp, int jmax, int imax, double dtau, double beta){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
	if(i <= imax-1 && j <= jmax-1){
        dev_pi[idx] = (double)(3.0/4.0) * dev_p[idx] + (double)(1.0/4.0) * (dev_pi[idx] + dtau * dev_rp[idx] * beta);
    }
}

__global__ void calc_3(double *dev_res_p, double *dev_pn, double *dev_rp, double *dev_p, double *dev_pi, int jmax, int imax, double dtau, double beta){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
	if(i <= imax-1 && j <= jmax-1){
        dev_res_p[idx] = dtau * dev_rp[idx] * beta;
        dev_pn[idx] = (double)(1.0 / 3.0) * dev_p[idx] + (double)(2.0 / 3.0) * (dev_pi[idx] + dev_res_p[idx]);
    }
}

double solve_P(double *dev_p, double *dev_um_n, double *dev_vm_n, double *dev_pn){
    dim3 blockDim(16,16);
    dim3 gridDimCalc1_3((imax-3 + blockDim.x - 1)/blockDim.x, (jmax-3 + blockDim.y - 1)/blockDim.y);
    dim3 gridDimCalc2((imax-3 + blockDim.x - 1)/blockDim.x, (jmax-4 + blockDim.y - 1)/blockDim.y);

    calc_1<<<gridDimCalc1_3, blockDim>>>(dev_dudx, dev_dvdy, dev_rp, dev_pi, dev_um_n, dev_vm_n, dev_areau_e, dev_areav_n, dev_areau_w, dev_areav_s, dev_p, imax, jmax, dtau, iterations.beta);

    bcP(dev_pi);

    calc_2<<<gridDimCalc2, blockDim>>>(dev_pi, dev_p, dev_rp, jmax, imax, dtau, iterations.beta);

    bcP(dev_pi);

    calc_3<<<gridDimCalc1_3, blockDim>>>(dev_res_p, dev_pn, dev_rp, dev_p, dev_pi, jmax, imax, dtau, iterations.beta);

    bcP(dev_pn);

    return max_reduce(dev_res_p, imax, jmax);
}
