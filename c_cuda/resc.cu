#include "comum.h"

//--- ResC ---
__global__ void calc_resc(double *dev_areau_e, double *dev_areau_w, double *dev_areav_n, double *dev_areav_s, double *dev_xm, 
    double *dev_ym, double *dev_x, double *dev_y, double *dev_liga_poros, double *dev_epsilon1, double re, double sc, int imax,  
    int jmax, double *dev_um_n, double *dev_vm_n, double *dev_c, double *dev_rc){
        
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
    int idx = i*(jmax+1)+j;
    double dcudx = 0.0, dcvdy = 0.0;
    double de, dw, dn, ds, dp;
    double aux = 1.0/re/sc;

    if(i <= imax-1 && j <= jmax-1){
        dcudx = 0.5 * (dev_c[(i+1)*(jmax+1)+j]+dev_c[idx]) * dev_um_n[(i+1)*(jmax+1)+j] * dev_areau_e[j]
                    - 0.5 * (dev_c[(i-1)*(jmax+1)+j]+dev_c[idx]) * dev_um_n[idx] * dev_areau_w[j];
        dcvdy = 0.5 * (dev_c[i*(jmax+1)+(j+1)]+dev_c[idx]) * dev_vm_n[i*(jmax+2)+(j+1)] * dev_areav_n[i]
                    - 0.5 * (dev_c[i*(jmax+1)+(j-1)]+dev_c[idx]) * dev_vm_n[i*(jmax+2)+j] * dev_areav_s[i];

        de = (dev_ym[j+1]-dev_ym[j]) * aux / (dev_x[i+1]-dev_x[i]);
        dw = (dev_ym[j+1]-dev_ym[j]) * aux / (dev_x[i]-dev_x[i-1]);
        dn = (dev_xm[i+1]-dev_xm[i]) * aux / (dev_y[j+1]-dev_y[j]);
        ds = (dev_xm[i+1]-dev_xm[i]) * aux / (dev_y[j]-dev_y[j-1]);
        dp = de + dw + dn + ds;

        dev_rc[idx] = 1.0 / (dev_xm[i+1]-dev_xm[i]) / (dev_ym[j+1]-dev_ym[j]) 
                *  (-dp*dev_c[idx] + de*dev_c[(i+1)*(jmax+1)+j] 
                +  dw*dev_c[(i-1)*(jmax+1)+j] + dn*dev_c[i*(jmax+1)+(j+1)] 
                +  ds*dev_c[i*(jmax+1)+(j-1)] - (dcudx + dcvdy) 
                /  (dev_liga_poros[idx]*(dev_epsilon1[idx]-1.0)+1.0));
    }
}

void RESC(double *dev_um_n, double *dev_vm_n, double *dev_c, double *dev_rc){
    dim3 blockDim(16,16);
    dim3 gridDim((imax-3 + blockDim.x - 1)/blockDim.x, (jmax-3 + blockDim.y - 1)/blockDim.y);
    
    calc_resc<<<gridDim, blockDim>>>(dev_areau_e, dev_areau_w, dev_areav_n, dev_areav_s, dev_xm, dev_ym, dev_x, dev_y, 
    dev_liga_poros, dev_epsilon1, re, sc, imax, jmax, dev_um_n, dev_vm_n, dev_c, dev_rc);
}