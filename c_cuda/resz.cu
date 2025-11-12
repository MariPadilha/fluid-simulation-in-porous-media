#include "comum.h"

//--- ResZ ---
__global__ void calc_resz(
    double *dev_xm, double *dev_x, double *dev_y, double *dev_ym, double *dev_areau_e, double *dev_areau_w, double *dev_areav_n,
    double *dev_areav_s, double *dev_epsilon1, double *dev_liga_poros, double pe, int imax, int jmax,
    double *dev_um_n, double *dev_vm_n, double *dev_z, double *dev_rz){

    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
    int idx = i*(jmax+1)+j;
    double de, dw, dn, ds, dp;
    double dzudx = 0.0, dzvdy = 0.0;
    double aux = 1.0/pe;

    if(i <= imax-1 && j <= jmax-1){
        dzudx = 0.5 * (dev_z[(i+1)*(jmax+1)+j]+dev_z[idx]) * dev_um_n[(i+1)*(jmax+1)+j] * dev_areau_e[j]
                    - 0.5 * (dev_z[(i-1)*(jmax+1)+j]+dev_z[idx]) * dev_um_n[idx] * dev_areau_w[j];
        dzvdy = 0.5 * (dev_z[i*(jmax+1)+(j+1)]+dev_z[idx]) * dev_vm_n[i*(jmax+2)+(j+1)] * dev_areav_n[i]
                    - 0.5 * (dev_z[i*(jmax+1)+(j-1)]+dev_z[idx]) * dev_vm_n[i*(jmax+2)+j] * dev_areav_s[i];

        de = (dev_ym[j+1]-dev_ym[j]) * aux / (dev_x[i+1]-dev_x[i]);  
        dw = (dev_ym[j+1]-dev_ym[j]) * aux / (dev_x[i]-dev_x[i-1]); 
        dn = (dev_xm[i+1]-dev_xm[i]) * aux / (dev_y[j+1]-dev_y[j]);  
        ds = (dev_xm[i+1]-dev_xm[i]) * aux / (dev_y[j]-dev_y[j-1]);  
        dp = de + dw + dn + ds;

        dev_rz[idx] = 1.0 / (dev_xm[i+1]-dev_xm[i]) / (dev_ym[j+1]-dev_ym[j]) 
                *  (-dp*dev_z[idx] +   de*dev_z[(i+1)*(jmax+1)+j] 
                +  dw*dev_z[(i-1)*(jmax+1)+j] + dn*dev_z[i*(jmax+1)+(j+1)] 
                +  ds*dev_z[i*(jmax+1)+(j-1)] - (1.0-dev_liga_poros[idx])
                *  (dzudx+dzvdy)) / (dev_liga_poros[idx]
                *  (dev_epsilon1[idx]-1.0)+1.0);
    }
}

void RESZ(double *dev_um_n, double *dev_vm_n, double *dev_z, double *dev_rz){
    dim3 blockDim(16, 16);
    dim3 gridDim((imax-3 + blockDim.x - 1)/blockDim.x, (jmax-3 + blockDim.y - 1)/blockDim.y);
    
    calc_resz<<<gridDim, blockDim>>>( 
    dev_xm, dev_x, dev_y, dev_ym, dev_areau_e, dev_areau_w, dev_areav_n,
    dev_areav_s, dev_epsilon1, dev_liga_poros, pe, imax, jmax, dev_um_n, dev_vm_n, dev_z, dev_rz);
}
