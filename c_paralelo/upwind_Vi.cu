#include "comum.h"
#define idx i*(jmax+1)+j

///fazer funcao device max
//--- upwind_V ---
__global__ void calc_upwind_Vi(double *dev_vm, double *dev_um, double *dev_areav_n, 
    double *dev_areav_s, double *dev_areav_e, double *dev_areav_w, double *dev_epsilon1,
    double *dev_ym, double *dev_x, double *dev_y, double *dev_xm, double *dev_liga_poros, double re,
    double *dev_p, double *dev_t, double *dev_rv, int imax, int jmax, 
    int j, double b_art, double invfr2, double darcy_number, double cf){

    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    double df, dn, ds, de, dw;
    double fn, fs, fe, fw;
    double aw, ae, as, an, ap;
    double v_w, v_e, v_s, v_n, v_p, u_p;
    double dvdydy, dydudx;
    double q_art;  

    if(i <= imax-1){
        fn = 0.5 * (dev_vm[i*(jmax+2)+j]+dev_vm[i*(jmax+2)+(j+1)]) * dev_areav_n[i] / dev_epsilon1[idx]; 
        fs = 0.5 * (dev_vm[i*(jmax+2)+j]+dev_vm[i*(jmax+2)+(j-1)]) * dev_areav_s[i] / dev_epsilon1[idx]; 
        fe = 0.5 * (dev_um[(i+1)*(jmax+1)+j]+dev_um[(i+1)*(jmax+1)+(j-1)]) * dev_areav_e[j] / dev_epsilon1[idx];  
        fw = 0.5 * (dev_um[idx] + dev_um[i*(jmax+1)+(j-1)]) * dev_areav_w[j] / dev_epsilon1[idx];
        df = fe - fw + fn - fs;
        dn = (dev_epsilon1[idx]/re) * dev_areav_n[i] / (dev_ym[j+1]-dev_ym[j]);
        ds = (dev_epsilon1[idx]/re) * dev_areav_s[i] / (dev_ym[j]-dev_ym[j-1]);
        de = (dev_epsilon1[idx]/re) * dev_areav_e[j] / (dev_x[i+1]-dev_x[i]);
        dw = (dev_epsilon1[idx]/re) * dev_areav_w[j] / (dev_x[i]-dev_x[i-1]);

        //upwind
        aw = dw + fmax(fw, 0.0);
        as = ds + fmax(fs, 0.0);
        ae = de + fmax(0.0, -fe);
        an = dn + fmax(0.0, -fn);
        ap = aw + ae + as + an + df;

        v_w = dev_vm[(i-1)*(jmax+2)+j];
        v_e = dev_vm[(i+1)*(jmax+2)+j];
        v_s = dev_vm[i*(jmax+2)+(j-1)];
        v_n = dev_vm[i*(jmax+2)+(j+1)];
        v_p = dev_vm[i*(jmax+2)+j];
        u_p = dev_um[idx];

        dvdydy = dev_areav_n[i] * (v_n-v_p) / (dev_ym[j+1]-dev_ym[j])
                     - dev_areav_s[i] * (v_p-v_s) / (dev_ym[j]-dev_ym[j-1]); 

        dydudx = dev_areav_n[i] * (dev_um[(i+1)*(jmax+1)+j]-dev_um[idx]) / (dev_xm[i+1]-dev_xm[i])
                     - dev_areav_s[i] * (dev_um[(i+1)*(jmax+1)+(j-1)]-dev_um[i*(jmax+1)+(j-1)]) / (dev_xm[i+1]-dev_xm[i]);

        //bulk artificial viscosity term from Ramshaw(1990)
        q_art = dev_epsilon1[idx] * (dev_p[idx]-dev_p[i*(jmax+1)+(j-1)]) / (dev_y[j]-dev_y[j-1]) - b_art * (dydudx + dvdydy);
        dev_rv[i*(jmax+2)+j] = 1.0 / (dev_x[i]-dev_x[i-1]) / (dev_y[j]-dev_y[j-1]) * (-ap * v_p
                +  aw * v_w + ae * v_e
                +  as * v_s + an * v_n)
                -  q_art + invfr2 * (1.0 - 1.0 / ((dev_t[idx]+dev_t[i*(jmax+1)+(j-1)]) * 0.5))
                -  dev_epsilon1[idx] * (v_p/(re*darcy_number)
                +  cf/(pow((dev_epsilon1[idx]*darcy_number), 0.5)) * v_p
                *  (pow((pow(u_p, 2.0) + pow(v_p, 2.0)), 0.5))) * dev_liga_poros[idx];
    }
}

void upwind_Vi(double *dev_um, double *dev_vm, double *dev_p, double *dev_rv, double *dev_t, int j){
    int threads = 256;
    dim3 blocks = grid_1d((imax-1), threads);

    calc_upwind_Vi<<<blocks, threads>>>(dev_vm, dev_um, dev_areav_n, dev_areav_s,  dev_areav_e,
        dev_areav_w, dev_epsilon1, dev_ym, dev_x, dev_y,
        dev_xm, dev_liga_poros,
        re, dev_p, dev_t, dev_rv, imax, jmax, j, iterations.b_art, invfr2, darcy_number, cf);

}
