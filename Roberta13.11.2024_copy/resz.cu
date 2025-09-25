#include "comum.h"
#define idx i*(jmax+1)+j


//--- ResZ ---
__global__ void calc_resz(double *dev_dzudx, double *dev_dzvdy, double *dev_de, double *dev_dw, double *dev_dn, double *dev_ds, double *dev_dp, 
    double *dev_xm, double *dev_x, double *dev_y, double *dev_ym, double *dev_areau_e, double *dev_areau_w, double *dev_areav_n,
    double *dev_areav_s, double *dev_epsilon1, double *dev_liga_poros, double pe, int imax, int jmax,
    double *dev_um_n, double *dev_vm_n, double **z, double *dev_rz){

    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
    
    if(i <= imax-1 && j <= jmax-1){
        dev_dzudx[idx] = 0.5 * (z[i+1][j]+z[i][j]) * dev_um_n[(i+1)*(jmax+1)+j] * dev_areau_e[j]
                    - 0.5 * (z[i-1][j]+z[i][j]) * dev_um_n[idx] * dev_areau_w[j];
        dev_dzvdy[idx] = 0.5 * (z[i][j+1]+z[i][j]) * dev_vm_n[i*(jmax+2)+(j+1)] * dev_areav_n[i]
                    - 0.5 * (z[i][j-1]+z[i][j]) * dev_vm_n[i*(jmax+2)+j] * dev_areav_s[i];

        dev_de[idx] = (dev_ym[j+1]-dev_ym[j]) * (1.0/pe) / (dev_x[i+1]-dev_x[i]);  
        dev_dw[idx] = (dev_ym[j+1]-dev_ym[j]) * (1.0/pe) / (dev_x[i]-dev_x[i-1]); 
        dev_dn[idx] = (dev_xm[i+1]-dev_xm[i]) * (1.0/pe) / (dev_y[j+1]-dev_y[j]);  
        dev_ds[idx] = (dev_xm[i+1]-dev_xm[i]) * (1.0/pe) / (dev_y[j]-dev_y[j-1]);  
        dev_dp[idx] = dev_de[idx] + dev_dw[idx] + dev_dn[idx] + dev_ds[idx];

        dev_rz[idx] = 1.0 / (dev_xm[i+1]-dev_xm[i]) / (dev_ym[j+1]-dev_ym[j]) 
                *  (-dev_dp[idx]*z[i][j] + dev_de[idx]*z[i+1][j] 
                +  dev_dw[idx]*z[i-1][j] + dev_dn[idx]*z[i][j+1] 
                +  dev_ds[idx]*z[i][j-1] - (1.0-dev_liga_poros[idx])
                *  (dev_dzudx[idx]+dev_dzvdy[idx])) / (dev_liga_poros[idx]
                *  (dev_epsilon1[idx]-1.0)+1.0);
    }
}

void RESZ(double *dev_um_n, double *dev_vm_n, double **z, double *dev_rz){
    int threads = 256;
    dim3 blocks = grid_1d((imax-1-2)*(jmax-1-2), threads);
    
    calc_resz<<<blocks, threads>>>(dev_dzudx, dev_dzvdy, dev_de, dev_dw, dev_dn, dev_ds, dev_dp, 
    dev_xm, dev_x, dev_y, dev_ym, dev_areau_e, dev_areau_w, dev_areav_n,
    dev_areav_s, dev_epsilon1, dev_liga_poros, pe, imax, jmax, dev_um_n, dev_vm_n, z, dev_rz);
}
