#include "comum.h"
#define idx i*(jmax+1)+j

__global__ void calc_upwind_Ui(double *dev_fn, double *dev_fs, double *dev_fe, double *dev_fw, double **vm, double **um, double *dev_areau_n, 
    double *dev_areau_s, double *dev_areau_e, double *dev_areau_w, double *dev_epsilon1, double *dev_df, double *dev_dn, double *dev_ds, 
    double *dev_de, double *dev_dw, double *dev_ym, double *dev_x, double *dev_y, double *dev_aw, double *dev_as, double *dev_ae, double *dev_an, 
    double *dev_ap, double *dev_u_w, double *dev_u_e, double *dev_u_s, double *dev_u_n, double *dev_v_p, double *dev_u_p, double *dev_dudxdx, 
    double *dev_dxdvdy, double *dev_q_art, double *dev_xm, double *dev_liga_poros, double re, double **p, double **ru, int imax, 
    int j, double b_art, double darcy_number, double cf, double g){
    
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;    

    if(i <= imax){
        //compute x-direction velocity component un
        dev_fn[idx] = 0.5 * (vm[i][j+1]+vm[i-1][j+1]) * dev_areau_n[i] / dev_epsilon1[idx];
        dev_fs[idx] = 0.5 * (vm[i][j]+vm[i-1][j]) * dev_areau_s[i] / dev_epsilon1[idx];
        dev_fe[idx] = 0.5 * (um[i+1][j]+um[i][j]) * dev_areau_e[j] / dev_epsilon1[idx];
        dev_fw[idx] = 0.5 * (um[i][j]+um[i-1][j]) * dev_areau_w[j] / dev_epsilon1[idx];

        dev_df[idx] = dev_fe[idx] - dev_fw[idx] + dev_fn[idx] - dev_fs[idx];

        dev_dn[idx] = (dev_epsilon1[idx]/re) * dev_areau_n[i] / (dev_y[j+1]-dev_y[j]);
        dev_ds[idx] = (dev_epsilon1[idx]/re) * dev_areau_s[i] / (dev_y[j]-dev_y[j-1]);
        dev_de[idx] = (dev_epsilon1[idx]/re) * dev_areau_e[j] / (dev_xm[i+1]-dev_xm[i]);
        dev_dw[idx] = (dev_epsilon1[idx]/re) * dev_areau_w[j] / (dev_xm[i]-dev_xm[i-1]);

        //upwind
        dev_aw[idx] = dev_dw[idx] + max(dev_fw[idx], 0.0);
        dev_as[idx] = dev_ds[idx] + max(dev_fs[idx], 0.0);
        dev_ae[idx] = dev_de[idx] + max(0.0, -dev_fe[idx]);
        dev_an[idx] = dev_dn[idx] + max(0.0, -dev_fn[idx]);

        dev_ap[idx] = dev_aw[idx] + dev_ae[idx] + dev_as[idx] + dev_an[idx] + dev_df[idx];

        dev_u_w[idx] = um[i-1][j];
        dev_u_e[idx] = um[i+1][j];
        dev_u_s[idx] = um[i][j-1];
        dev_u_n[idx] = um[i][j+1];
        dev_u_p[idx] = um[i][j];
        dev_v_p[idx] = vm[i][j];

        dev_dudxdx[idx] = dev_areau_e[j] * (dev_u_e[idx]-dev_u_p[idx]) / (dev_xm[i+1]-dev_xm[i])
                     - dev_areau_w[j] * (dev_u_p[idx]-dev_u_w[idx]) / (dev_xm[i]-dev_xm[i-1]); 
    
        dev_dxdvdy[idx] = dev_areau_e[j] * (vm[i][j+1]-vm[i][j]) / (dev_ym[j+1]-dev_ym[j])
                     - dev_areau_w[j] * (vm[i-1][j+1]-vm[i-1][j]) / (dev_ym[j+1]-dev_ym[j]);

        //bulk artificial viscosity term from Ramshaw(1990)
        dev_q_art[idx] = dev_epsilon1[idx] * (p[i][j]-p[i-1][j]) / (dev_x[i]-dev_x[i-1]) - iterations.b_art * (dev_dudxdx[idx]+dev_dxdvdy[idx]);

        ru[i][j] = 1.0 / (dev_x[i]-dev_x[i-1]) / (dev_y[j]-dev_y[j-1]) * (-dev_ap[idx]*dev_u_p[idx]
                +  dev_aw[idx] * dev_u_w[idx] + dev_ae[idx] * dev_u_e[idx]
                +  dev_as[idx] * dev_u_s[idx] + dev_an[idx] * dev_u_n[idx])
                -  dev_q_art[idx] - dev_epsilon1[idx] * (dev_u_p[idx]/(re*darcy_number) 
                +  cf/pow((dev_epsilon1[idx]*darcy_number),0.5) * dev_u_p[idx]
                *  (pow((pow(dev_u_p[idx],2.0) + pow(dev_v_p[idx],2.0)),0.5))) 
                *  dev_liga_poros[idx] - g * dev_epsilon1[idx];
    }
}

//--- upwind_U ---
void upwind_Ui(double **um, double **vm, double **p, double **ru, int j){
    int threads = 256;
    dim3 blocks = grid_1d(imax, threads);

    calc_upwind_Ui<<<blocks, threads>>>(dev_fn, dev_fs, dev_fe, dev_fw, vm, um, dev_areau_n, dev_areau_s,
         dev_areau_e, dev_areau_w, dev_epsilon1, dev_df, dev_dn, dev_ds, dev_de, dev_dw, dev_ym, dev_x, dev_y, dev_aw,
         dev_as, dev_ae, dev_an, dev_ap, dev_u_w, dev_u_e, dev_u_s, dev_u_n, dev_v_p, dev_u_p, dev_dudxdx, dev_dxdvdy,
         dev_q_art, dev_xm, dev_liga_poros, re, p, ru, imax, j, iterations.b_art, darcy_number, cf, g);

    cudaDeviceSynchronize();
}