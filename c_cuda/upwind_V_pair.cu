#include "comum.h"

// Fusão SEGURA de Vi + Vj em um kernel
__global__ void calc_upwind_V_pair(
    double *dev_vm, double *dev_um, double *dev_areav_n, double *dev_areav_s, 
    double *dev_areav_e, double *dev_areav_w, double *dev_epsilon1, 
    double *dev_ym, double *dev_x, double *dev_y, double *dev_xm, double *dev_liga_poros, 
    double re, double *dev_p, double *dev_t, double *dev_rv, int imax, int jmax,
    double b_art, double invfr2, double darcy_number, double cf, 
    int mode, int fixed_index) {
    
    int thread_id = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int i, j;
    
    if (mode == 0) {         
        i = thread_id;
        j = fixed_index;
        if (i > imax-1) return;
    } else {            
        i = fixed_index;
        j = thread_id;
        if (j > jmax-1) return;
    }
    
    int idx = i*(jmax+1)+j;
    
    double epsilon_idx = dev_epsilon1[idx];
    double aux = epsilon_idx/re;
    double inv_epsilon = 1.0 / epsilon_idx;
    
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
    double x_ip1 = dev_x[i+1];
    double x_i = dev_x[i];
    double x_im1 = dev_x[i-1];
    double y_j = dev_y[j];
    double y_jm1 = dev_y[j-1];
    
    double fn = 0.5 * (dev_vm[i*(jmax+2)+j] + dev_vm[i*(jmax+2)+(j+1)]) * areav_n_i * inv_epsilon;
    double fs = 0.5 * (dev_vm[i*(jmax+2)+j] + dev_vm[i*(jmax+2)+(j-1)]) * areav_s_i * inv_epsilon;
    double fe = 0.5 * (dev_um[(i+1)*(jmax+1)+j] + dev_um[(i+1)*(jmax+1)+(j-1)]) * areav_e_j * inv_epsilon;
    double fw = 0.5 * (dev_um[idx] + dev_um[i*(jmax+1)+(j-1)]) * areav_w_j * inv_epsilon;
    
    double df = fe - fw + fn - fs;
    
    double dy_n = ym_jp1 - ym_j;
    double dy_s = ym_j - ym_jm1;
    double dx_e = x_ip1 - x_i;
    double dx_w = x_i - x_im1;
    
    double dn = aux * areav_n_i / dy_n;
    double ds = aux * areav_s_i / dy_s;
    double de = aux * areav_e_j / dx_e;
    double dw = aux * areav_w_j / dx_w;
    
    double aw = dw + fmax(fw, 0.0);
    double as = ds + fmax(fs, 0.0);
    double ae = de + fmax(0.0, -fe);
    double an = dn + fmax(0.0, -fn);
    double ap = aw + ae + as + an + df;
    
    double v_w = dev_vm[(i-1)*(jmax+2)+j];
    double v_e = dev_vm[(i+1)*(jmax+2)+j];
    double v_s = dev_vm[i*(jmax+2)+(j-1)];
    double v_n = dev_vm[i*(jmax+2)+(j+1)];
    double v_p = dev_vm[i*(jmax+2)+j];
    double u_p = dev_um[idx];
    
    double dvdydy = areav_n_i * (v_n - v_p) / dy_n - areav_s_i * (v_p - v_s) / dy_s;
    
    double dx_xm = xm_ip1 - xm_i;
    double um_e = dev_um[(i+1)*(jmax+1)+j];
    double um_es = dev_um[(i+1)*(jmax+1)+(j-1)];
    double um_ws = dev_um[i*(jmax+1)+(j-1)];
    
    double dydudx = areav_n_i * (um_e - u_p) / dx_xm - areav_s_i * (um_es - um_ws) / dx_xm;
    
    double dy_main = y_j - y_jm1;
    double dx_main = x_i - x_im1;
    double q_art = epsilon_idx * (dev_p[idx] - dev_p[i*(jmax+1)+(j-1)]) / dy_main - b_art * (dydudx + dvdydy);
    
    double inv_vol = 1.0 / (dx_main * dy_main);
    double velocidade_mag = sqrt(u_p*u_p + v_p*v_p);
    double darcy_term = v_p / (re * darcy_number);
    double forchheimer_coef = cf / sqrt(epsilon_idx * darcy_number);
    double porous_drag = epsilon_idx * (darcy_term + forchheimer_coef * v_p * velocidade_mag) * dev_liga_poros[idx];
    double temp_avg = (dev_t[idx] + dev_t[i*(jmax+1)+(j-1)]) * 0.5;
    double buoyancy_term = invfr2 * (1.0 - 1.0 / temp_avg);
    
    dev_rv[i*(jmax+2)+j] = inv_vol * (-ap * v_p + aw * v_w + ae * v_e + as * v_s + an * v_n)
            - q_art + buoyancy_term - porous_drag;
}

void upwind_V_pair(double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_rv) {
    int threads = 256;
    
    // Vi para j=2 (varia i de 2 até imax-1)
    dim3 blocks_i = dim3((imax-3 + threads - 1) / threads);
    calc_upwind_V_pair<<<blocks_i, threads>>>(
        dev_vm, dev_um, dev_areav_n, dev_areav_s, dev_areav_e, dev_areav_w, 
        dev_epsilon1, dev_ym, dev_x, dev_y, dev_xm, dev_liga_poros, re, 
        dev_p, dev_t, dev_rv, imax, jmax, iterations.b_art, invfr2, darcy_number, cf, 0, 2);
    
    // Vi para j=jmax-1 (varia i de 2 até imax-1)
    calc_upwind_V_pair<<<blocks_i, threads>>>(
        dev_vm, dev_um, dev_areav_n, dev_areav_s, dev_areav_e, dev_areav_w,
        dev_epsilon1, dev_ym, dev_x, dev_y, dev_xm, dev_liga_poros, re,
        dev_p, dev_t, dev_rv, imax, jmax, iterations.b_art, invfr2, darcy_number, cf, 0, jmax-1);
    
    // Vj para i=2 (varia j de 2 até jmax-1)
    dim3 blocks_j = dim3((jmax-3 + threads - 1) / threads);
    calc_upwind_V_pair<<<blocks_j, threads>>>(
        dev_vm, dev_um, dev_areav_n, dev_areav_s, dev_areav_e, dev_areav_w,
        dev_epsilon1, dev_ym, dev_x, dev_y, dev_xm, dev_liga_poros, re,
        dev_p, dev_t, dev_rv, imax, jmax, iterations.b_art, invfr2, darcy_number, cf, 1, 2);
    
    // Vj para i=imax-1 (varia j de 2 até jmax-1)  
    calc_upwind_V_pair<<<blocks_j, threads>>>(
        dev_vm, dev_um, dev_areav_n, dev_areav_s, dev_areav_e, dev_areav_w,
        dev_epsilon1, dev_ym, dev_x, dev_y, dev_xm, dev_liga_poros, re,
        dev_p, dev_t, dev_rv, imax, jmax, iterations.b_art, invfr2, darcy_number, cf, 1, imax-1);
}