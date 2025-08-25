#include "comum.h"

__global__ void calc_vol_u(){
    for(j = 1 ; j <= jmax; j++){
        for(i = 2; i <= imax; i++){
            vol_u[i][j] = (dev_x[i]-dev_x[i-1]) * (dev_ym[j+1]-dev_ym[j]); 
        }
    }
}

__global__ void calc_vol_v(){
    for(j = 2; j <= jmax; j++){
        for(i = 1; i <= imax; i++){
            vol_v[i][j] = (dev_y[j]-dev_y[j-1]) * (dev_xm[i+1]-dev_xm[i]); 
        }
    }
}

__global__ void calc_vol_p(){
    for(j = 2; j <= jmax; j++){
        for(i = 1; i <= imax; i++){
            vol_p[i][j] = (dev_ym[j+1]-dev_ym[j]) * (dev_xm[i+1]-dev_xm[i]); 
        }
    }
}

void calcula_vol_u_v_p(){
    int threads = 256;
    dim3 blocks = grid_2d(imax, threads);

}