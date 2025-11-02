#include "comum.h"

//--- upwind_U ---
__global__ void calc_upwind_Uj(double *dev_vm, double *dev_um, double *dev_areau_n, 
    double *dev_areau_s, double *dev_areau_e, double *dev_areau_w, double *dev_epsilon1,
    double *dev_ym, double *dev_x, double *dev_y,
    double *dev_xm, double *dev_liga_poros, double re, double *dev_p, double *dev_ru, int jmax, int imax,
    int i, double b_art, double darcy_number, double cf, double g){

    int j = blockIdx.x * blockDim.x + threadIdx.x + 2;  
    int idx = i*(jmax+1)+j;
    double df, dn, ds, de, dw;
    double fn, fs, fe, fw;
    double aw, ae, as, an, ap;
    double u_w, u_e, u_s, u_n, u_p, v_p;
    double dudxdx, dxdvdy;
    double q_art;
    double aux = dev_epsilon1[idx]/re;
    
    if(j <= jmax-1){
        //compute x-direction velocity component un
        fn = 0.5 * (dev_vm[i*(jmax+2)+(j+1)]+dev_vm[(i-1)*(jmax+2)+(j+1)]) * dev_areau_n[i] / dev_epsilon1[idx];
        fs = 0.5 * (dev_vm[i*(jmax+2)+j]+dev_vm[(i-1)*(jmax+2)+j]) * dev_areau_s[i] / dev_epsilon1[idx];
        fe = 0.5 * (dev_um[(i+1)*(jmax+1)+j]+dev_um[idx]) * dev_areau_e[j] / dev_epsilon1[idx];
        fw = 0.5 * (dev_um[idx]+dev_um[(i-1)*(jmax+1)+j]) * dev_areau_w[j] / dev_epsilon1[idx];

        df = fe - fw + fn - fs;
        dn = aux * dev_areau_n[i] / (dev_y[j+1]-dev_y[j]);
        ds = aux * dev_areau_s[i] / (dev_y[j]-dev_y[j-1]);
        de = aux * dev_areau_e[j] / (dev_xm[i+1]-dev_xm[i]);
        dw = aux * dev_areau_w[j] / (dev_xm[i]-dev_xm[i-1]);

        //upwind
        aw = dw + fmax(fw, 0.0);
        as = ds + fmax(fs, 0.0);
        ae = de + fmax(0.0, -fe);
        an = dn + fmax(0.0, -fn);
        ap = aw + ae + as + an + df;

        u_w = dev_um[(i-1)*(jmax+1)+j];
        u_e = dev_um[(i+1)*(jmax+1)+j];
        u_s = dev_um[i*(jmax+1)+(j-1)];
        u_n = dev_um[i*(jmax+1)+(j+1)];
        u_p = dev_um[idx];
        v_p = dev_vm[i*(jmax+2)+j];

        dudxdx = dev_areau_e[j] * (u_e-u_p) / (dev_xm[i+1]-dev_xm[i])
                     - dev_areau_w[j] * (u_p-u_w) / (dev_xm[i]-dev_xm[i-1]); 

        dxdvdy = dev_areau_e[j] * (dev_vm[i*(jmax+2)+(j+1)]-dev_vm[i*(jmax+2)+j]) / (dev_ym[j+1]-dev_ym[j])
                     - dev_areau_w[j] * (dev_vm[(i-1)*(jmax+2)+(j+1)]-dev_vm[(i-1)*(jmax+2)+j]) / (dev_ym[j+1]-dev_ym[j]);

        //bulk artificial viscosity term from Ramshaw(1990)
        q_art = dev_epsilon1[idx] * (dev_p[idx]-dev_p[(i-1)*(jmax+1)+j]) / (dev_x[i]-dev_x[i-1])
                    - (b_art) * (dudxdx+dxdvdy);

        dev_ru[idx] = 1.0 / (dev_x[i]-dev_x[i-1]) / (dev_y[j]-dev_y[j-1]) * (-ap * u_p
                +  aw * u_w + ae * u_e
                +  as * u_s + an * u_n)
                -  q_art - dev_epsilon1[idx] * (u_p/(re*darcy_number)
                +  cf/(pow((dev_epsilon1[idx]*darcy_number),0.5)) * u_p
                *  (pow((pow(u_p,2.0) + pow(v_p,2.0)),0.5)))
                *  dev_liga_poros[idx] - g * dev_epsilon1[idx];
    }
}

void upwind_Uj(double *dev_um, double *dev_vm, double *dev_p, double *dev_ru, int i){
    int threads = 256;
    dim3 blocks = grid_1d((jmax-1), threads);

    calc_upwind_Uj<<<blocks, threads>>>(dev_vm, dev_um, dev_areau_n, dev_areau_s,
         dev_areau_e, dev_areau_w, dev_epsilon1, dev_ym, dev_x, dev_y,
         dev_xm, dev_liga_poros, re, dev_p, dev_ru, jmax, imax, i, iterations.b_art, darcy_number, cf, g);

}
