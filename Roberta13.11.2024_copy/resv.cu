#include "comum.h"
#define idx 1

__global__ void calc_resv_parte1(){
    for(j = 3; j <= jmax-1; j++){
        for(i = 3; i <= imax-2; i++){
            fn[i][j] = 0.5 * (vm[i][j]+vm[i][j+1]) * areav_n[i] / epsilon1[i][j];
            fs[i][j] = 0.5 * (vm[i][j]+vm[i][j-1]) * areav_s[i] / epsilon1[i][j];
            fe[i][j] = 0.5 * (um[i+1][j]+um[i+1][j-1]) * areav_e[j] / epsilon1[i][j];
            fw[i][j] = 0.5 * (um[i][j]+um[i][j-1]) * areav_w[j] / epsilon1[i][j];

            df[i][j] = fe[i][j] - fw[i][j] + fn[i][j] - fs[i][j];

            dn[i][j] = (epsilon1[i][j]/re) * areav_n[i] / (ym[j+1]-ym[j]);
            ds[i][j] = (epsilon1[i][j]/re) * areav_s[i] / (ym[j]-ym[j-1]);
            de[i][j] = (epsilon1[i][j]/re) * areav_e[j] / (x[i+1]-x[i]);
            dw[i][j] = (epsilon1[i][j]/re) * areav_w[j] / (x[i]-x[i-1]);
        }
    }
}

__global__ void calc_resv_parte2(){
    for(j = 3; j <= jmax-1; j++){
        for(i = 3; i <= imax-2; i++){
            aw[i][j] = dw[i][j] + 0.75 * afw[i][j] * fw[i][j]
                    +  0.125 * afe[i][j] * fe[i][j]
                    +  0.375 * (1.0-afw[i][j]) * fw[i][j];

            ae[i][j] = de[i][j] - 0.375 * afe[i][j] * fe[i][j]
                    - 0.75 * (1.0-afe[i][j]) * fe[i][j]
                    - 0.125 * (1.0-afw[i][j]) * fw[i][j];

            as[i][j] = ds[i][j] + 0.75  * afs[i][j] * fs[i][j] 
                    +  0.125 * afn[i][j] * fn[i][j] 
                    +  0.375 * (1.0-afs[i][j]) * fs[i][j];

            an[i][j] = dn[i][j] - 0.375 * afn[i][j] * fn[i][j]
                    - 0.75 * (1.0-afn[i][j]) * fn[i][j]
                    - 0.125 * (1.0-afs[i][j]) * fs[i][j];

            aww[i][j] = -0.125 * afw[i][j] * fw[i][j];
            aee[i][j] =  0.125 * (1.0-afe[i][j]) * fe[i][j];
            ass[i][j] = -0.125 * afs[i][j] * fs[i][j];
            ann[i][j] =  0.125 * (1.0-afn[i][j]) * fn[i][j];

            ap[i][j] = aw[i][j] + ae[i][j] + as[i][j] + an[i][j] 
                    +  aww[i][j] + aee[i][j] + ass[i][j] + ann[i][j] + df[i][j];
            //end Quick///////////////////////////////////////////////

            v_w[i][j]  = vm[i-1][j];
            v_ww[i][j] = vm[i-2][j];
            v_e[i][j]  = vm[i+1][j];
            v_ee[i][j] = vm[i+2][j];
            v_s[i][j]  = vm[i][j-1];
            v_ss[i][j] = vm[i][j-2];
            v_n[i][j]  = vm[i][j+1];
            v_nn[i][j] = vm[i][j+2];         
            v_p[i][j]  = vm[i][j];
            u_p[i][j]  = um[i][j];

            dvdydy[i][j] = areav_n[i] * (v_n[i][j]-v_p[i][j]) / (ym[j+1]-ym[j])
                        -  areav_s[i] * (v_p[i][j]-v_s[i][j]) / (ym[j]-ym[j-1]);

            dydudx[i][j] = areav_n[i] * (um[i+1][j]-um[i][j]) / (xm[i+1]-xm[i])
                        -  areav_s[i] * (um[i+1][j-1]-um[i][j-1]) / (xm[i+1]-xm[i]);

            artdivv[i][j] = -(iterations.b_art) * (dydudx[i][j]+dvdydy[i][j]);            

            //bulk artificial viscosity term from Ramshaw(1990)
            q_art[i][j] = epsilon1[i][j] * (p[i][j]-p[i][j-1]) / (y[j]-y[j-1]) + artdivv[i][j];

            rv[i][j] = 1.0 / (x[i]-x[i-1]) / (y[j]-y[j-1]) 
                    * (-ap[i][j] * v_p[i][j] + aww[i][j] * v_ww[i][j] + aw[i][j] 
                    * v_w[i][j] + aee[i][j] * v_ee[i][j] + ae[i][j] * v_e[i][j] 
                    + ass[i][j] * v_ss[i][j] + as[i][j] * v_s[i][j] 
                    + ann[i][j] * v_nn[i][j] + an[i][j] * v_n[i][j]) 
                    - q_art[i][j] + invfr2 * (1.0 - 1.0 / ((t[i][j]+t[i][j-1]) * 0.5))
                    - epsilon1[i][j]*(v_p[i][j]/(re*darcy_number) 
                    + cf/(pow((epsilon1[i][j]*darcy_number), 0.5)) * v_p[i][j] 
                    * (pow((pow(u_p[i][j], 2.0) + pow(v_p[i][j], 2.0)), 0.5))) * liga_poros[i][j]; 
        }
    }
}


//--- ResV ---
void RESV(double **um, double **vm, double **p, double **t, double **rv){
    int i, j;

    calc_resv_parte1<<<thread, blocks>>>

    for(j = 3; j <= jmax-1; j++){
        for(i = 3; i <= imax-2; i++){
            //quick
            if(fw[i][j] > 0.0){
                afw[i][j] = 1.0;
            }else if(fw[i][j] < 0.0){
                afw[i][j] = 0.0;
            }

            if(fe[i][j] > 0.0){
                afe[i][j] = 1.0;
            }else if(fe[i][j] < 0.0){
                afe[i][j] = 0.0;
            }

            if(fn[i][j] > 0.0){
                afn[i][j] = 1.0;
            }else if(fn[i][j] < 0.0){
                afn[i][j] = 0.0;
            }

            if(fs[i][j] > 0.0){
                afs[i][j] = 1.0;
            }else if(fs[i][j] < 0.0){
                afs[i][j] = 0.0;
            }
    }

    calc_resv_parte2<<<thread, blocks>>>

    upwind_Vi(um,vm,p,rv,t,2);
    upwind_Vi(um,vm,p,rv,t,jmax);    
    upwind_Vj(um,vm,p,rv,t,2);
    upwind_Vj(um,vm,p,rv,t,imax-1);
}