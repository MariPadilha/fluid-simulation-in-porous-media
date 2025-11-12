#ifndef COMUM_H
#define COMUM_H

#include <math.h>
#include "functions.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <device_launch_parameters.h>

struct Iterations{
    int itc_max, nc, n_tr, n_out, n_vort;
    double beta, b_art, dtau_f, final_time, eps, eps_mass;
    int start_mode;
};
extern struct Iterations iterations;

struct Ref{
    double tnu, yf_b, yo_oo, ts, tn_too;
};
extern struct Ref ref;

extern int restart_mode;
extern double dtau, dt, tempo;
extern double porosidade, darcy_number, cf, temp_cylinder, concentracao_inicial;
extern double px_grid, py_grid, q_grid, lhori, y_up, y_down, hvert;
extern int imax, jmax;    
extern double dx_c;
extern double *dev_x, *dev_y, *dev_xm, *dev_ym;
extern double *dev_vol_u, *dev_vol_v, *dev_vol_p;
extern double *dev_areau_n, *dev_areau_s,*dev_areau_e, *dev_areau_w, *dev_areav_n, *dev_areav_s, *dev_areav_e, *dev_areav_w;       
extern double *liga_poros, *epsilon1, *dev_epsilon1, *dev_liga_poros;
extern double rad1;
extern int c_i, c_b, c_f, c_bs;     
extern int *dev_flag;
extern double g, ao, l_c, v_i, v_c, fr, invfr2, s, lf, lo;
extern double too, tsup, tinf, q_dim, q;
extern double cp_tot, rho_tot, k_tot, nu_tot, alpha_tot, re, pr, pe, sc;
extern double *dev_dudx, *dev_dvdy, *dev_rp, *dev_pi;
extern double *dev_res_p, *dev_res_z, *dev_res_c;

extern double **dcdx2, **dcdy2;
extern double *dev_ci, *dev_zi, *dev_dx, *dev_dy; 

dim3 grid_2d(int imax, int jmax);
void calcular(int n_imax, int n_itc);
void alocar_globais();
void desalocar_globais();
#endif