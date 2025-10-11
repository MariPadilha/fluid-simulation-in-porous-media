#include "comum.h"
#define idx i*(jmax+1)+j

__global__ void calc_resv(double *dev_fn, double *dev_fs, double *dev_fe, double *dev_fw, double *dev_df, double *dev_dn, double *dev_ds, double *dev_de, 
	double *dev_dw, double *dev_areav_e, double *dev_areav_n, double *dev_areav_s, double *dev_areav_w, double *dev_epsilon1, double *dev_y, double *dev_x,
	double *dev_xm, double *dev_afw, double *dev_afe, double *dev_afn, double *dev_afs, double *dev_aw, double *dev_ae, double *dev_as, double *dev_an,
        double *dev_aww, double *dev_aee, double *dev_ass, double *dev_ann, double *dev_ap, double *dev_v_e, double *dev_v_ee, double *dev_v_n, double *dev_v_nn, 
	double *dev_v_p, double *dev_v_s, double *dev_v_ss, double *dev_v_w, double *dev_v_ww, double *dev_u_p, double *dev_ym, double *dev_dvdydy, double *dev_dydudx,
	double *dev_q_art, double *dev_artdivv, double *dev_liga_poros, double b_art, int imax, int jmax, double re, double darcy_number, double cf, double invfr2,
	double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_rv){

    int i = blockIdx.x * blockDim.x + threadIdx.x + 3;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
    if(i <= imax-2 && j <= jmax-1){
        dev_fn[idx] = 0.5 * (dev_vm[i*(jmax+2)+j]+dev_vm[i*(jmax+2)+(j+1)]) * dev_areav_n[i] / dev_epsilon1[idx];
        dev_fs[idx] = 0.5 * (dev_vm[i*(jmax+2)+j]+dev_vm[i*(jmax+2)+(j-1)]) * dev_areav_s[i] / dev_epsilon1[idx];
        dev_fe[idx] = 0.5 * (dev_um[(i+1)*(jmax+1)+j]+dev_um[(i+1)*(jmax+1)+(j-1)]) * dev_areav_e[j] / dev_epsilon1[idx];
        dev_fw[idx] = 0.5 * (dev_um[idx]+dev_um[i*(jmax+1)+(j-1)]) * dev_areav_w[j] / dev_epsilon1[idx];

        dev_df[idx] = dev_fe[idx] - dev_fw[idx] + dev_fn[idx] - dev_fs[idx];

        dev_dn[idx] = (dev_epsilon1[idx]/re) * dev_areav_n[i] / (dev_ym[j+1]-dev_ym[j]);
        dev_ds[idx] = (dev_epsilon1[idx]/re) * dev_areav_s[i] / (dev_ym[j]-dev_ym[j-1]);
        dev_de[idx] = (dev_epsilon1[idx]/re) * dev_areav_e[j] / (dev_x[i+1]-dev_x[i]);
        dev_dw[idx] = (dev_epsilon1[idx]/re) * dev_areav_w[j] / (dev_x[i]-dev_x[i-1]);

        //quick
        dev_afw[idx] = (double)(dev_fw[idx] > 0.0);
        dev_afe[idx] = (double)(dev_fe[idx] > 0.0);
        dev_afn[idx] = (double)(dev_fn[idx] > 0.0);
        dev_afs[idx] = (double)(dev_fs[idx] > 0.0);
        
        
        dev_aw[idx] = dev_dw[idx] + 0.75 * dev_afw[idx] * dev_fw[idx]
                +  0.125 * dev_afe[idx] * dev_fe[idx]
                +  0.375 * (1.0-dev_afw[idx]) * dev_fw[idx];

        dev_ae[idx] = dev_de[idx] - 0.375 * dev_afe[idx] * dev_fe[idx]
                - 0.75 * (1.0-dev_afe[idx]) * dev_fe[idx]
                - 0.125 * (1.0-dev_afw[idx]) * dev_fw[idx];

        dev_as[idx] = dev_ds[idx] + 0.75  * dev_afs[idx] * dev_fs[idx] 
                +  0.125 * dev_afn[idx] * dev_fn[idx] 
                +  0.375 * (1.0-dev_afs[idx]) * dev_fs[idx];

        dev_an[idx] = dev_dn[idx] - 0.375 * dev_afn[idx] * dev_fn[idx]
                - 0.75 * (1.0-dev_afn[idx]) * dev_fn[idx]
                - 0.125 * (1.0-dev_afs[idx]) * dev_fs[idx];

        dev_aww[idx] = -0.125 * dev_afw[idx] * dev_fw[idx];
        dev_aee[idx] =  0.125 * (1.0-dev_afe[idx]) * dev_fe[idx];
        dev_ass[idx] = -0.125 * dev_afs[idx] * dev_fs[idx];
        dev_ann[idx] =  0.125 * (1.0-dev_afn[idx]) * dev_fn[idx];

        dev_ap[idx] = dev_aw[idx] + dev_ae[idx] + dev_as[idx] + dev_an[idx] 
                +  dev_aww[idx] + dev_aee[idx] + dev_ass[idx] + dev_ann[idx] + dev_df[idx];
        //end Quick///////////////////////////////////////////////

        dev_v_w[idx]  = dev_vm[(i-1)*(jmax+2)+j];
        dev_v_ww[idx] = dev_vm[(i-2)*(jmax+2)+j];
        dev_v_e[idx]  = dev_vm[(i+1)*(jmax+2)+j];
        dev_v_ee[idx] = dev_vm[(i+2)*(jmax+2)+j];
        dev_v_s[idx]  = dev_vm[i*(jmax+2)+(j-1)];
        dev_v_ss[idx] = dev_vm[i*(jmax+2)+(j-2)];
        dev_v_n[idx]  = dev_vm[i*(jmax+2)+(j+1)];
        dev_v_nn[idx] = dev_vm[i*(jmax+2)+(j+2)];         
        dev_v_p[idx]  = dev_vm[i*(jmax+2)+j];
        dev_u_p[idx]  = dev_um[idx];

        dev_dvdydy[idx] = dev_areav_n[i] * (dev_v_n[idx]-dev_v_p[idx]) / (dev_ym[j+1]-dev_ym[j])
                    -  dev_areav_s[i] * (dev_v_p[idx]-dev_v_s[idx]) / (dev_ym[j]-dev_ym[j-1]);

        dev_dydudx[idx] = dev_areav_n[i] * (dev_um[(i+1)*(jmax+1)+j]-dev_um[idx]) / (dev_xm[i+1]-dev_xm[i])
                    -  dev_areav_s[i] * (dev_um[(i+1)*(jmax+1)+(j-1)]-dev_um[i*(jmax+1)+(j-1)]) / (dev_xm[i+1]-dev_xm[i]);

        dev_artdivv[idx] = -(b_art) * (dev_dydudx[idx]+dev_dvdydy[idx]);            

        //bulk artificial viscosity term from Ramshaw(1990)
        dev_q_art[idx] = dev_epsilon1[idx] * (dev_p[idx]-dev_p[i*(jmax+1)+(j-1)]) / (dev_y[j]-dev_y[j-1]) + dev_artdivv[idx];

        dev_rv[i*(jmax+2)+j] = 1.0 / (dev_x[i]-dev_x[i-1]) / (dev_y[j]-dev_y[j-1]) 
                * (-dev_ap[idx] * dev_v_p[idx] + dev_aww[idx] * dev_v_ww[idx] + dev_aw[idx] 
                * dev_v_w[idx] + dev_aee[idx] * dev_v_ee[idx] + dev_ae[idx] * dev_v_e[idx] 
                + dev_ass[idx] * dev_v_ss[idx] + dev_as[idx] * dev_v_s[idx] 
                + dev_ann[idx] * dev_v_nn[idx] + dev_an[idx] * dev_v_n[idx]) 
                - dev_q_art[idx] + invfr2 * (1.0 - 1.0 / ((dev_t[idx]+dev_t[i*(jmax+1)+(j-1)]) * 0.5))
                - dev_epsilon1[idx]*(dev_v_p[idx]/(re*darcy_number) 
                + cf/(pow((dev_epsilon1[idx]*darcy_number), 0.5)) * dev_v_p[idx] 
                * (pow((pow(dev_u_p[idx], 2.0) + pow(dev_v_p[idx], 2.0)), 0.5))) * dev_liga_poros[idx]; 
    }
}

//--- ResV ---
void RESV(double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_rv){
    dim3 blockDim(16,16);
    dim3 gridDim((imax-5 + blockDim.x - 1)/blockDim.x, (jmax-4 + blockDim.y - 1)/blockDim.y);

    calc_resv<<<gridDim, blockDim>>>(dev_fn, dev_fs, dev_fe, dev_fw, dev_df, dev_dn, dev_ds,
    dev_de, dev_dw, dev_areav_e, dev_areav_n, dev_areav_s, dev_areav_w, dev_epsilon1, dev_y, 
    dev_x, dev_xm, dev_afw, dev_afe, dev_afn, dev_afs, dev_aw, dev_ae, dev_as, 
    dev_an, dev_aww, dev_aee, dev_ass, dev_ann, dev_ap, dev_v_e, dev_v_ee, dev_v_n, 
    dev_v_nn, dev_v_p, dev_v_s, dev_v_ss, dev_v_w, dev_v_ww, dev_u_p, dev_ym, dev_dvdydy,
    dev_dydudx, dev_q_art, dev_artdivv, dev_liga_poros, iterations.b_art, imax, jmax, re, darcy_number, 
    cf, invfr2, dev_um, dev_vm, dev_p, dev_t, dev_rv);

    for(int i = 1; i <= imax; i++){
        for(int j = 1; j <= jmax; j++){
            printf("[%i][%i]epilson = %lf\n", i, j, dev_epsilon1[i*(jmax+1)+j]);
        }
    }

    upwind_Vi(dev_um,dev_vm, dev_p,dev_rv, dev_t, 2);
    upwind_Vi(dev_um,dev_vm, dev_p,dev_rv, dev_t, jmax);    
    upwind_Vj(dev_um,dev_vm, dev_p,dev_rv, dev_t, 2);
    upwind_Vj(dev_um,dev_vm, dev_p,dev_rv, dev_t, imax-1);
}