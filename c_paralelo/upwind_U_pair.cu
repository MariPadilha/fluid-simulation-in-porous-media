#include "comum.h"

// Fusão SEGURA de Ui + Uj em um kernel
__global__ void calc_upwind_U_pair(
    double *dev_vm, double *dev_um, double *dev_areau_n, double *dev_areau_s, 
    double *dev_areau_e, double *dev_areau_w, double *dev_epsilon1, 
    double *dev_ym, double *dev_x, double *dev_y, double *dev_xm, double *dev_liga_poros, 
    double re, double *dev_p, double *dev_ru, int imax, int jmax,
    double b_art, double darcy_number, double cf, double g, 
    int mode, int fixed_index) {
    
    int thread_id = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int i, j;
    
    if (mode == 0) {         
        i = thread_id;
        j = fixed_index;
        if (i > imax) return;
    } else {                
        i = fixed_index;
        j = thread_id;
        if (j > jmax-1) return;
    }
    
    int idx = i*(jmax+1)+j;
    
    double epsilon_idx = dev_epsilon1[idx];
    double aux = epsilon_idx/re;
    double inv_epsilon = 1.0 / epsilon_idx;
    
    double areau_e_j = dev_areau_e[j];
    double areau_w_j = dev_areau_w[j];
    double areau_n_i = dev_areau_n[i];
    double areau_s_i = dev_areau_s[i];
    
    double xm_ip1 = dev_xm[i+1];
    double xm_i = dev_xm[i];
    double xm_im1 = dev_xm[i-1];
    double y_jp1 = dev_y[j+1];
    double y_j = dev_y[j];
    double y_jm1 = dev_y[j-1];
    double x_i = dev_x[i];
    double x_im1 = dev_x[i-1];
    double ym_jp1 = dev_ym[j+1];
    double ym_j = dev_ym[j];
    
    double fn = 0.5 * (dev_vm[i*(jmax+2)+(j+1)] + dev_vm[(i-1)*(jmax+2)+(j+1)]) * areau_n_i * inv_epsilon;
    double fs = 0.5 * (dev_vm[i*(jmax+2)+j] + dev_vm[(i-1)*(jmax+2)+j]) * areau_s_i * inv_epsilon;
    double fe = 0.5 * (dev_um[(i+1)*(jmax+1)+j] + dev_um[idx]) * areau_e_j * inv_epsilon;
    double fw = 0.5 * (dev_um[idx] + dev_um[(i-1)*(jmax+1)+j]) * areau_w_j * inv_epsilon;
    
    double df = fe - fw + fn - fs;
    
    double dy_n = y_jp1 - y_j;
    double dy_s = y_j - y_jm1;
    double dx_e = xm_ip1 - xm_i;
    double dx_w = xm_i - xm_im1;
    
    double dn = aux * areau_n_i / dy_n;
    double ds = aux * areau_s_i / dy_s;
    double de = aux * areau_e_j / dx_e;
    double dw = aux * areau_w_j / dx_w;
    
    double aw = dw + fmax(fw, 0.0);
    double as = ds + fmax(fs, 0.0);
    double ae = de + fmax(0.0, -fe);
    double an = dn + fmax(0.0, -fn);
    double ap = aw + ae + as + an + df;
    
    double u_w = dev_um[(i-1)*(jmax+1)+j];
    double u_e = dev_um[(i+1)*(jmax+1)+j];
    double u_s = dev_um[i*(jmax+1)+(j-1)];
    double u_n = dev_um[i*(jmax+1)+(j+1)];
    double u_p = dev_um[idx];
    double v_p = dev_vm[i*(jmax+2)+j];
    
    double dudxdx = areau_e_j * (u_e - u_p) / dx_e - areau_w_j * (u_p - u_w) / dx_w;
    
    double dy_ym = ym_jp1 - ym_j;
    double vm_n = dev_vm[i*(jmax+2)+(j+1)];
    double vm_s = dev_vm[i*(jmax+2)+j];
    double vm_nw = dev_vm[(i-1)*(jmax+2)+(j+1)];
    double vm_sw = dev_vm[(i-1)*(jmax+2)+j];
    
    double dxdvdy = areau_e_j * (vm_n - vm_s) / dy_ym - areau_w_j * (vm_nw - vm_sw) / dy_ym;
    
    double dx_main = x_i - x_im1;
    double dy_main = y_j - y_jm1;
    double q_art = epsilon_idx * (dev_p[idx] - dev_p[(i-1)*(jmax+1)+j]) / dx_main - b_art * (dudxdx + dxdvdy);
    
    double inv_vol = 1.0 / (dx_main * dy_main);
    double velocidade_mag = sqrt(u_p*u_p + v_p*v_p);
    double darcy_term = u_p / (re * darcy_number);
    double forchheimer_coef = cf / sqrt(epsilon_idx * darcy_number);
    double porous_drag = epsilon_idx * (darcy_term + forchheimer_coef * u_p * velocidade_mag) * dev_liga_poros[idx];
    
    dev_ru[idx] = inv_vol * (-ap * u_p + aw * u_w + ae * u_e + as * u_s + an * u_n)
            - q_art - porous_drag - g * epsilon_idx;
}

void upwind_U_pair(double *dev_um, double *dev_vm, double *dev_p, double *dev_ru) {
    int threads = 256;
    
    // Ui para j=2 (varia i de 2 até imax)
    dim3 blocks_i = dim3((imax-2 + threads - 1) / threads);
    calc_upwind_U_pair<<<blocks_i, threads>>>(
        dev_vm, dev_um, dev_areau_n, dev_areau_s, dev_areau_e, dev_areau_w, 
        dev_epsilon1, dev_ym, dev_x, dev_y, dev_xm, dev_liga_poros, re, 
        dev_p, dev_ru, imax, jmax, iterations.b_art, darcy_number, cf, g, 0, 2);
    
    // Ui para j=jmax-1 (varia i de 2 até imax)
    calc_upwind_U_pair<<<blocks_i, threads>>>(
        dev_vm, dev_um, dev_areau_n, dev_areau_s, dev_areau_e, dev_areau_w,
        dev_epsilon1, dev_ym, dev_x, dev_y, dev_xm, dev_liga_poros, re,
        dev_p, dev_ru, imax, jmax, iterations.b_art, darcy_number, cf, g, 0, jmax-1);
    
    // Uj para i=2 (varia j de 2 até jmax-1)
    dim3 blocks_j = dim3((jmax-3 + threads - 1) / threads);
    calc_upwind_U_pair<<<blocks_j, threads>>>(
        dev_vm, dev_um, dev_areau_n, dev_areau_s, dev_areau_e, dev_areau_w,
        dev_epsilon1, dev_ym, dev_x, dev_y, dev_xm, dev_liga_poros, re,
        dev_p, dev_ru, imax, jmax, iterations.b_art, darcy_number, cf, g, 1, 2);
    
    // Uj para i=imax (varia j de 2 até jmax-1)  
    calc_upwind_U_pair<<<blocks_j, threads>>>(
        dev_vm, dev_um, dev_areau_n, dev_areau_s, dev_areau_e, dev_areau_w,
        dev_epsilon1, dev_ym, dev_x, dev_y, dev_xm, dev_liga_poros, re,
        dev_p, dev_ru, imax, jmax, iterations.b_art, darcy_number, cf, g, 1, imax);
}