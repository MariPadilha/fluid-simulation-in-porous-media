#include "comum.h"
#define idx 1


//--- ResZ ---
__global__ void calc_resz(double **um_n, double **vm_n, double **z, double **rz){
    int i, j;

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            dzudx[i][j] = 0.5 * (z[i+1][j]+z[i][j]) * um_n[i+1][j] * areau_e[j]
                        - 0.5 * (z[i-1][j]+z[i][j]) * um_n[i][j] * areau_w[j];
            dzvdy[i][j] = 0.5 * (z[i][j+1]+z[i][j]) * vm_n[i][j+1] * areav_n[i]
                        - 0.5 * (z[i][j-1]+z[i][j]) * vm_n[i][j] * areav_s[i];

            de[i][j] = (ym[j+1]-ym[j]) * (1.0/pe) / (x[i+1]-x[i]);  
            dw[i][j] = (ym[j+1]-ym[j]) * (1.0/pe) / (x[i]-x[i-1]); 
            dn[i][j] = (xm[i+1]-xm[i]) * (1.0/pe) / (y[j+1]-y[j]);  
            ds[i][j] = (xm[i+1]-xm[i]) * (1.0/pe) / (y[j]-y[j-1]);  
            dp[i][j] = de[i][j] + dw[i][j] + dn[i][j] + ds[i][j];

            rz[i][j] = 1.0 / (xm[i+1]-xm[i]) / (ym[j+1]-ym[j]) 
                    *  (-dp[i][j]*z[i][j] + de[i][j]*z[i+1][j] 
                    +  dw[i][j]*z[i-1][j] + dn[i][j]*z[i][j+1] 
                    +  ds[i][j]*z[i][j-1] - (1.0-liga_poros[i][j])
                    *  (dzudx[i][j]+dzvdy[i][j])) / (liga_poros[i][j]
                    *  (epsilon1[i][j]-1.0)+1.0);
        }
    }
}

void RESZ(double **um_n, double **vm_n, double **z, double **rz){
    calc_resz<<<thread, blocks>>>
}
