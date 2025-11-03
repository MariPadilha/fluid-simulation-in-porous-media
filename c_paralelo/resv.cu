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
    
    double epsilon_idx = dev_epsilon1[idx];
    double aux = epsilon_idx/re;
    double areav_e_j = dev_areav_e[j];
    double areav_w_j = dev_areav_w[j]; 
    double areav_n_i = dev_areav_n[i];
    double areav_s_i = dev_areav_s[i];
    double xm_ip1 = dev_xm[i+1];
    double xm_i = dev_xm[i];
    double xm_im1 = dev_xm[i-1];
    double ym_jp1 = dev_ym[j+1];
    double ym_j = dev_ym[j];
    double ym_jm1 = dev_ym[j-1];
    double y_j = dev_y[j];
    double y_jm1 = dev_y[j-1];
    double x_i = dev_x[i];
    double x_im1 = dev_x[i-1];

    if(i <= imax-2 && j <= jmax-1){
        double inv_epsilon = 1.0 / epsilon_idx;
        fn = 0.5 * (dev_vm[i*(jmax+2)+j]+dev_vm[i*(jmax+2)+(j+1)]) * areav_n_i * inv_epsilon;
        fs = 0.5 * (dev_vm[i*(jmax+2)+j]+dev_vm[i*(jmax+2)+(j-1)]) * areav_s_i * inv_epsilon;
        fe = 0.5 * (dev_um[(i+1)*(jmax+1)+j]+dev_um[(i+1)*(jmax+1)+(j-1)]) * areav_e_j * inv_epsilon;
        fw = 0.5 * (dev_um[idx]+dev_um[i*(jmax+1)+(j-1)]) * areav_w_j * inv_epsilon;

        df = fe - fw + fn - fs;
        
        double dy_n = ym_jp1 - ym_j;
        double dy_s = ym_j - ym_jm1;
        double dx_e = dev_x[i+1] - x_i;
        double dx_w = x_i - x_im1;
        
        dn = aux * areav_n_i / dy_n;
        ds = aux * areav_s_i / dy_s;
        de = aux * areav_e_j / dx_e;
        dw = aux * areav_w_j / dx_w;

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

        dvdydy = areav_n_i * (v_n-v_p) / dy_n - areav_s_i * (v_p-v_s) / dy_s;

        double dx_xm = xm_ip1 - xm_i;
        double um_e = dev_um[(i+1)*(jmax+1)+j];
        double um_es = dev_um[(i+1)*(jmax+1)+(j-1)];
        double um_ws = dev_um[i*(jmax+1)+(j-1)];
        
        dydudx = areav_n_i * (um_e - u_p) / dx_xm - areav_s_i * (um_es - um_ws) / dx_xm;

        artdivv = -(b_art) * (dydudx+dvdydy);

        //bulk artificial viscosity term from Ramshaw(1990)
        double dy_main = y_j - y_jm1;
        double dx_main = x_i - x_im1;
        q_art = epsilon_idx * (dev_p[idx]-dev_p[i*(jmax+1)+(j-1)]) / dy_main + artdivv;

        double inv_vol = 1.0 / (dx_main * dy_main);
        double velocidade_mag = sqrt(u_p*u_p + v_p*v_p);
        double darcy_term = v_p / (re * darcy_number);
        double forchheimer_coef = cf / sqrt(epsilon_idx * darcy_number);
        double porous_drag = epsilon_idx * (darcy_term + forchheimer_coef * v_p * velocidade_mag) * dev_liga_poros[idx];
        double temp_avg = (dev_t[idx] + dev_t[i*(jmax+1)+(j-1)]) * 0.5;
        double buoyancy_term = invfr2 * (1.0 - 1.0 / temp_avg);

        dev_rv[i*(jmax+2)+j] = inv_vol * (-ap * v_p + aww * v_ww + aw * v_w + aee * v_ee + ae * v_e
                + ass * v_ss + as * v_s + ann * v_nn + an * v_n) 
                - q_art + buoyancy_term - porous_drag; 
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


    upwind_V_pair(dev_um, dev_vm, dev_p, dev_t, dev_rv);
}