#include "comum.h"
#define idx i*(jmax+1)+j

__global__ void calc_upwind_ui(double **um, double **vm, double **p, double **ru, int j){

    for(i = 2; i <= imax; i++){
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

}