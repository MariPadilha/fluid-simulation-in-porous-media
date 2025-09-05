#include "comum.h"
#define idx i*(jmax+1)+j

//--- upwind_V ---
__global__ void calc_upwind_Vi(double *dev_fn, double *dev_fs, double *dev_fe, double *dev_fw, double **vm, double **um, double *dev_areav_n, 
    double *dev_areav_s, double *dev_areav_e, double *dev_areav_w, double *dev_epsilon1, double *dev_df, double *dev_dn, double *dev_ds, 
    double *dev_de, double *dev_dw, double *dev_ym, double *dev_x, double *dev_y, double *dev_aw, double *dev_as, double *dev_ae, double *dev_an, 
    double *dev_ap, double *dev_v_w, double *dev_v_e, double *dev_v_s, double *dev_v_n, double *dev_v_p, double *dev_u_p, double *dev_dvdydy, 
    double *dev_dydudx, double *dev_q_art, double *dev_xm, double *dev_liga_poros, double re, double **p, double **t, double **rv, int imax, 
    int j, double b_art, double invfr2, double darcy_number, double cf){
    
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;    

    if(i <= imax-1){
        dev_fn[idx] = 0.5 * (vm[i][j]+vm[i][j+1]) * dev_areav_n[i] / dev_epsilon1[idx]; 
        dev_fs[idx] = 0.5 * (vm[i][j]+vm[i][j-1]) * dev_areav_s[i] / dev_epsilon1[idx]; 
        dev_fe[idx] = 0.5 * (um[i+1][j]+um[i+1][j-1]) * dev_areav_e[j] / dev_epsilon1[idx];  
        dev_fw[idx] = 0.5 * (um[i][j] + um[i][j-1]) * dev_areav_w[j] / dev_epsilon1[idx];

        dev_df[idx] = dev_fe[idx] - dev_fw[idx] + dev_fn[idx] - dev_fs[idx];

        dev_dn[idx] = (dev_epsilon1[idx]/re) * dev_areav_n[i] / (dev_ym[j+1]-dev_ym[j]);
        dev_ds[idx] = (dev_epsilon1[idx]/re) * dev_areav_s[i] / (dev_ym[j]-dev_ym[j-1]);
        dev_de[idx] = (dev_epsilon1[idx]/re) * dev_areav_e[j] / (dev_x[i+1]-dev_x[i]);
        dev_dw[idx] = (dev_epsilon1[idx]/re) * dev_areav_w[j] / (dev_x[i]-dev_x[i-1]);

        //upwind
        dev_aw[idx] = dev_dw[idx] + max(dev_fw[idx], 0.0);
        dev_as[idx] = dev_ds[idx] + max(dev_fs[idx], 0.0);
        dev_ae[idx] = dev_de[idx] + max(0.0, -dev_fe[idx]);
        dev_an[idx] = dev_dn[idx] + max(0.0, -dev_fn[idx]);

        dev_ap[idx] = dev_aw[idx] + dev_ae[idx] + dev_as[idx] + dev_an[idx] + dev_df[idx];

        dev_v_w[idx] = vm[i-1][j];
        dev_v_e[idx] = vm[i+1][j];
        dev_v_s[idx] = vm[i][j-1];
        dev_v_n[idx] = vm[i][j+1];
        dev_v_p[idx] = vm[i][j];
        dev_u_p[idx] = um[i][j];

        dev_dvdydy[idx] = dev_areav_n[i] * (dev_v_n[idx]-dev_v_p[idx]) / (dev_ym[j+1]-dev_ym[j])
                     - dev_areav_s[i] * (dev_v_p[idx]-dev_v_s[idx]) / (dev_ym[j]-dev_ym[j-1]); 
            
        dev_dydudx[idx] = dev_areav_n[i] * (um[i+1][j]-um[i][j]) / (dev_xm[i+1]-dev_xm[i])
                     - dev_areav_s[i] * (um[i+1][j-1]-um[i][j-1]) / (dev_xm[i+1]-dev_xm[i]);
            
        //bulk artificial viscosity term from Ramshdev_aw(1dx)
        dev_q_art[idx] = dev_epsilon1[idx] * (p[i][j]-p[i][j-1]) / (dev_y[j]-dev_y[j-1]) - b_art * (dev_dydudx[idx]+dev_dvdydy[idx]);
        rv[i][j] = 1.0 / (dev_x[i]-dev_x[i-1]) / (dev_y[j]-dev_y[j-1]) * (-dev_ap[idx] * dev_v_p[idx]
                +  dev_aw[idx] * dev_v_w[idx] + dev_ae[idx] * dev_v_e[idx] 
                +  dev_as[idx] * dev_v_s[idx] + dev_an[idx] * dev_v_n[idx])  
                -  dev_q_art[idx] + invfr2 * (1.0 - 1.0 / ((t[i][j]+t[i][j-1]) * 0.5))
                -  dev_epsilon1[idx] * (dev_v_p[idx]/(re*darcy_number) 
                +  cf/(pow((dev_epsilon1[idx]*darcy_number), 0.5)) * dev_v_p[idx]
                *  (pow((pow(dev_u_p[idx], 2.0) + pow(dev_v_p[idx], 2.0)), 0.5))) * dev_liga_poros[idx]; 
    }
}

void upwind_Vi(double **um, double **vm, double **p, double **rv, double **t, int j){
    int threads = 256;
    dim3 blocks = grid_1d((imax-1), threads);

    calc_upwind_Vi<<<blocks, threads>>>(dev_fn, dev_fs, dev_fe, dev_fw, vm, um, dev_areav_n, dev_areav_s,  dev_areav_e,
        dev_areav_w, dev_epsilon1, dev_df, dev_dn, dev_ds, dev_de, dev_dw, dev_ym, dev_x, dev_y, dev_aw, dev_as, dev_ae, dev_an, 
        dev_ap, dev_v_w, dev_v_e, dev_v_s, dev_v_n, dev_v_p, dev_u_p, dev_dvdydy, dev_dydudx, dev_q_art, dev_xm, dev_liga_poros,
        re, p, t, rv, imax, j, iterations.b_art, invfr2, darcy_number, cf);

    cudaDeviceSynchronize();
}
