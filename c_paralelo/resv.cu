#include "comum.h"

__global__ void calc_resv(
	double *dev_areav_e, double *dev_areav_n, double *dev_areav_s, double *dev_areav_w, double *dev_epsilon1, double *dev_y, double *dev_x,
	double *dev_xm, double *dev_ym, double *dev_liga_poros, double b_art, int imax, int jmax, double re, double darcy_number, double cf, double invfr2,
	double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_rv){

    int i = blockIdx.x * blockDim.x + threadIdx.x + 3;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
    int idx = i*(jmax+1)+j;
    double df, dn, ds, de, dw;
    double fn, fs, fe, fw;
    double afw, afe, afn, afs;
    double aw, ae, as, an, ap;
    double aww, aee, ass, ann;
    double v_w, v_e, v_n, v_s, v_p, u_p;
    double v_ww, v_ee, v_nn, v_ss;
    double dydudx, dvdydy;
    double q_art, artdivv;
    double aux = dev_epsilon1[idx]/re;

    if(i <= imax-2 && j <= jmax-1){
        fn = 0.5 * (dev_vm[i*(jmax+2)+j]+dev_vm[i*(jmax+2)+(j+1)]) * dev_areav_n[i] / dev_epsilon1[idx];
        fs = 0.5 * (dev_vm[i*(jmax+2)+j]+dev_vm[i*(jmax+2)+(j-1)]) * dev_areav_s[i] / dev_epsilon1[idx];
        fe = 0.5 * (dev_um[(i+1)*(jmax+1)+j]+dev_um[(i+1)*(jmax+1)+(j-1)]) * dev_areav_e[j] / dev_epsilon1[idx];
        fw = 0.5 * (dev_um[idx]+dev_um[i*(jmax+1)+(j-1)]) * dev_areav_w[j] / dev_epsilon1[idx];

        df = fe - fw + fn - fs;
        dn = aux * dev_areav_n[i] / (dev_ym[j+1]-dev_ym[j]);
        ds = aux * dev_areav_s[i] / (dev_ym[j]-dev_ym[j-1]);
        de = aux * dev_areav_e[j] / (dev_x[i+1]-dev_x[i]);
        dw = aux * dev_areav_w[j] / (dev_x[i]-dev_x[i-1]);

        //quick
        afw = (double)(fw > 0.0);
        afe = (double)(fe > 0.0);
        afn = (double)(fn > 0.0);
        afs = (double)(fs > 0.0);


        aw = dw + 0.75 * afw * fw
                +  0.125 * afe * fe
                +  0.375 * (1.0-afw) * fw;

        ae = de - 0.375 * afe * fe
                - 0.75 * (1.0-afe) * fe
                - 0.125 * (1.0-afw) * fw;

        as = ds + 0.75  * afs * fs 
                +  0.125 * afn * fn 
                +  0.375 * (1.0-afs) * fs;

        an = dn - 0.375 * afn * fn
                - 0.75 * (1.0-afn) * fn
                - 0.125 * (1.0-afs) * fs;

        aww = -0.125 * afw * fw;
        aee =  0.125 * (1.0-afe) * fe;
        ass = -0.125 * afs * fs;
        ann =  0.125 * (1.0-afn) * fn;

        ap = aw + ae + as + an + aww + aee + ass + ann + df;
        //end Quick///////////////////////////////////////////////

        v_w  = dev_vm[(i-1)*(jmax+2)+j];
        v_ww = dev_vm[(i-2)*(jmax+2)+j];
        v_e  = dev_vm[(i+1)*(jmax+2)+j];
        v_ee = dev_vm[(i+2)*(jmax+2)+j];
        v_s  = dev_vm[i*(jmax+2)+(j-1)];
        v_ss = dev_vm[i*(jmax+2)+(j-2)];
        v_n  = dev_vm[i*(jmax+2)+(j+1)];
        v_nn = dev_vm[i*(jmax+2)+(j+2)];
        v_p  = dev_vm[i*(jmax+2)+j];
        u_p  = dev_um[idx];

        dvdydy = dev_areav_n[i] * (v_n-v_p) / (dev_ym[j+1]-dev_ym[j])
                    -  dev_areav_s[i] * (v_p-v_s) / (dev_ym[j]-dev_ym[j-1]);

        dydudx = dev_areav_n[i] * (dev_um[(i+1)*(jmax+1)+j]-dev_um[idx]) / (dev_xm[i+1]-dev_xm[i])
                    -  dev_areav_s[i] * (dev_um[(i+1)*(jmax+1)+(j-1)]-dev_um[i*(jmax+1)+(j-1)]) / (dev_xm[i+1]-dev_xm[i]);

        artdivv = -(b_art) * (dydudx+dvdydy);

        //bulk artificial viscosity term from Ramshaw(1990)
        q_art = dev_epsilon1[idx] * (dev_p[idx]-dev_p[i*(jmax+1)+(j-1)]) / (dev_y[j]-dev_y[j-1]) + artdivv;

        dev_rv[i*(jmax+2)+j] = 1.0 / (dev_x[i]-dev_x[i-1]) / (dev_y[j]-dev_y[j-1]) 
                * (-ap * v_p + aww * v_ww + aw * v_w + aee * v_ee + ae * v_e
                + ass * v_ss + as * v_s
                + ann * v_nn + an * v_n) 
                - q_art + invfr2 * (1.0 - 1.0 / ((dev_t[idx]+dev_t[i*(jmax+1)+(j-1)]) * 0.5))
                - dev_epsilon1[idx]*(v_p/(re*darcy_number) 
                + cf/(pow((dev_epsilon1[idx]*darcy_number), 0.5)) * v_p 
                * (pow((pow(u_p, 2.0) + pow(v_p, 2.0)), 0.5))) * dev_liga_poros[idx]; 
    }
}

//--- ResV ---
void RESV(double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_rv){
    dim3 blockDim(16,16);
    dim3 gridDim((imax-5 + blockDim.x - 1)/blockDim.x, (jmax-4 + blockDim.y - 1)/blockDim.y);

    calc_resv<<<gridDim, blockDim>>>(
    dev_areav_e, dev_areav_n, dev_areav_s, dev_areav_w, dev_epsilon1, dev_y, 
    dev_x, dev_xm, dev_ym, dev_liga_poros, iterations.b_art, imax, jmax, re, darcy_number, 
    cf, invfr2, dev_um, dev_vm, dev_p, dev_t, dev_rv);


    upwind_Vi(dev_um,dev_vm, dev_p,dev_rv, dev_t, 2);
    upwind_Vi(dev_um,dev_vm, dev_p,dev_rv, dev_t, jmax);    
    upwind_Vj(dev_um,dev_vm, dev_p,dev_rv, dev_t, 2);
    upwind_Vj(dev_um,dev_vm, dev_p,dev_rv, dev_t, imax-1);
}