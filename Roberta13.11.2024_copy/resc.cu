#include "comum.h"
#define idx 1

//--- ResC ---
__global__ void calc_resc(double **um_n, double **vm_n, double **c, double **rc){
    int i, j;

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            dcudx[i][j] = 0.5 * (c[i+1][j]+c[i][j]) * um_n[i+1][j] * dev_areau_e[j]
                        - 0.5 * (c[i-1][j]+c[i][j]) * um_n[i][j] * dev_areau_w[j];
            dcvdy[i][j] = 0.5 * (c[i][j+1]+c[i][j]) * vm_n[i][j+1] * dev_areav_n[i]
                        - 0.5 * (c[i][j-1]+c[i][j]) * vm_n[i][j] * dev_areav_s[i];

            dev_de[i][j] = (ym[j+1]-ym[j]) * (1.0/re/sc) / (dev_x[i+1]-dev_x[i]);
            dev_dw[i][j] = (ym[j+1]-ym[j]) * (1.0/re/sc) / (dev_x[i]-dev_x[i-1]);
            dev_dn[i][j] = (xm[i+1]-xm[i]) * (1.0/re/sc) / (dev_y[j+1]-dev_y[j]);
            dev_ds[i][j] = (xm[i+1]-xm[i]) * (1.0/re/sc) / (dev_y[j]-dev_y[j-1]);
            dp[i][j] = de[i][j] + dw[i][j] + dn[i][j] + ds[i][j];

            rc[i][j] = 1.0 / (xm[i+1]-xm[i]) / (ym[j+1]-ym[j]) 
                    *  (-dp[i][j]*c[i][j] + de[i][j]*c[i+1][j] 
                    +  dw[i][j]*c[i-1][j] + dn[i][j]*c[i][j+1] 
                    +  ds[i][j]*c[i][j-1] - (dcudx[i][j] + dcvdy[i][j]) 
                    /  (liga_poros[i][j]*(epsilon1[i][j]-1.0)+1.0));
        }
    }
}

void RESC(double **um_n, double **vm_n, double **c, double **rc){
    calc_resc<<<thread, block>>>
}