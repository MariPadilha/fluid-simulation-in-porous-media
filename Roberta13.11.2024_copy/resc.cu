#include "comum.h"
#define idx i*(jmax+1)+j

//--- ResC ---
__global__ void calc_resc(double *dev_dcudx, double *dev_dcvdy, double *dev_areau_e, double *dev_areau_w, double *dev_areav_n, double *dev_areav_s, 
    double *dev_de, double *dev_dw, double *dev_dn, double *dev_ds, double *dev_dp, double *dev_xm, double *dev_ym, double *dev_x, double *dev_y, 
    double *dev_liga_poros, double *dev_epsilon1, double re, double sc, int imax,  int jmax, double *dev_um_n, double *dev_vm_n, double **c, double *dev_rc){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;

    if(i <= imax-1 && j <= jmax-1){
        dev_dcudx[idx] = 0.5 * (c[i+1][j]+c[i][j]) * dev_um_n[(i+1)*(jmax+1)+j] * dev_areau_e[j]
                    - 0.5 * (c[i-1][j]+c[i][j]) * dev_um_n[idx] * dev_areau_w[j];
        dev_dcvdy[idx] = 0.5 * (c[i][j+1]+c[i][j]) * dev_vm_n[i*(jmax+2)+(j+1)] * dev_areav_n[i]
                    - 0.5 * (c[i][j-1]+c[i][j]) * dev_vm_n[i*(jmax+2)+j] * dev_areav_s[i];

        dev_de[idx] = (dev_ym[j+1]-dev_ym[j]) * (1.0/re/sc) / (dev_x[i+1]-dev_x[i]);
        dev_dw[idx] = (dev_ym[j+1]-dev_ym[j]) * (1.0/re/sc) / (dev_x[i]-dev_x[i-1]);
        dev_dn[idx] = (dev_xm[i+1]-dev_xm[i]) * (1.0/re/sc) / (dev_y[j+1]-dev_y[j]);
        dev_ds[idx] = (dev_xm[i+1]-dev_xm[i]) * (1.0/re/sc) / (dev_y[j]-dev_y[j-1]);
        dev_dp[idx] = dev_de[idx] + dev_dw[idx] + dev_dn[idx] + dev_ds[idx];

        dev_rc[idx] = 1.0 / (dev_xm[i+1]-dev_xm[i]) / (dev_ym[j+1]-dev_ym[j]) 
                *  (-dev_dp[idx]*c[i][j] + dev_de[idx]*c[i+1][j] 
                +  dev_dw[idx]*c[i-1][j] + dev_dn[idx]*c[i][j+1] 
                +  dev_ds[idx]*c[i][j-1] - (dev_dcudx[idx] + dev_dcvdy[idx]) 
                /  (dev_liga_poros[idx]*(dev_epsilon1[idx]-1.0)+1.0));
    }
}

void RESC(double *dev_um_n, double *dev_vm_n, double **c, double *dev_rc){
    int threads = 256;
    dim3 blocks = grid_1d((imax-1-2)*(jmax-1-2), threads);
    
    calc_resc<<<blocks, threads>>>(dev_dcudx, dev_dcvdy, dev_areau_e, dev_areau_w, dev_areav_n, dev_areav_s, 
    dev_de, dev_dw, dev_dn, dev_ds, dev_dp, dev_xm, dev_ym, dev_x, dev_y, 
    dev_liga_poros, dev_epsilon1, re, sc, imax, jmax, dev_um_n, dev_vm_n, c, dev_rc);
}