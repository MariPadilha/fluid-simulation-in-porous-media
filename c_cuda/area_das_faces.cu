#include "comum.h"

__global__ void calc_areau_n_s(double *dev_x, double *dev_areau_n, double *dev_areau_s, int imax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if(i <= imax){
        dev_areau_n[i] = dev_x[i] - dev_x[i-1];
        dev_areau_s[i] = dev_x[i] - dev_x[i-1];
    }
}

__global__ void calc_areav_n_s(double *dev_xm, double *dev_areav_n, double *dev_areav_s, int imax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    if(i <= imax){
        dev_areav_n[i] = dev_xm[i+1] - dev_xm[i];
        dev_areav_s[i] = dev_xm[i+1] - dev_xm[i];
    }
}

__global__ void calc_areau_e_w(double *dev_ym, double *dev_areau_e, double *dev_areau_w, int jmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    if(i <= jmax){
        dev_areau_e[i] = dev_ym[i+1] - dev_ym[i];
        dev_areau_w[i] = dev_ym[i+1] - dev_ym[i];
    }
}

__global__ void calc_areav_e_w(double *dev_y, double *dev_areav_e, double *dev_areav_w, int jmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    if(i <= jmax){
        dev_areav_e[i] = dev_y[i] - dev_y[i-1];
        dev_areav_w[i] = dev_y[i] - dev_y[i-1];
    }
}

void calcula_area_das_fases(){
    int threads = 256;
    int blocks = grid_1d(imax-1, threads);
    calc_areau_n_s<<<blocks, threads>>>(dev_x, dev_areau_n, dev_areau_s, imax);

    blocks = grid_1d(jmax, threads);
    calc_areau_e_w<<<blocks, threads>>>(dev_ym, dev_areau_e, dev_areau_w, jmax);
    
    blocks = grid_1d(imax, threads);
    calc_areav_n_s<<<blocks, threads>>>(dev_xm, dev_areav_n, dev_areav_s, imax);

    blocks = grid_1d(jmax-1, threads);
    calc_areav_e_w<<<blocks, threads>>>(dev_y, dev_areav_e, dev_areav_w, jmax);
}