#ifndef COMUM_H
#define COMUM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

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
extern double dtau, dt, time;
extern double porosidade, darcy_number, cf, temp_cylinder, concentracao_inicial;
extern double px_grid, py_grid, q_grid, lhori, y_up, y_down, hvert;
extern int imax, jmax;    
extern double dx_c;
extern double *dev_x, *dev_y, *dev_xm, *dev_ym;
extern double *dev_vol_u, *dev_vol_v, *dev_vol_p;
extern double *dev_areau_n, *dev_areau_s,*dev_areau_e, *dev_areau_w, *dev_areav_n, *dev_areav_s, *dev_areav_e, *dev_areav_w;       
extern double *dev_dx, *dev_dy;
extern double *liga_poros, *epsilon1, *dev_epsilon1, *dev_liga_poros;
extern double rad1;
extern int c_i, c_b, c_f, c_bs;     
extern int **flag;
extern double g, ao, l_c, v_i, v_c, fr, invfr2, s, lf, lo;
extern double too, tsup, tinf, q_dim, q;
extern double cp_tot, rho_tot, k_tot, nu_tot, alpha_tot, re, pr, pe, sc;
extern double *dev_fw, *dev_fe, *dev_fs, *dev_fn;
extern double *dev_df, *dev_dn, *dev_ds, *dev_de, *dev_dw;
extern double *dev_aw, *dev_as, *dev_ae, *dev_an, *dev_ap;
extern double *dev_u_w, *dev_u_e, *dev_u_s, *dev_u_n, *dev_u_p, *dev_v_p;
extern double *dev_dudxdx, *dev_dxdvdy;
extern double *dev_q_art;
extern double *dev_dvdydy, *dev_dydudx;
extern double *dev_v_w, *dev_v_e, *dev_v_s, *dev_v_n;
extern double *dev_aww, *dev_aee, *dev_ass;
extern double *dev_ann, *dev_u_ww, *dev_artdivu;
extern double *dev_u_ee, *dev_u_ss, *dev_u_nn;
extern double *dev_afw, *dev_afe, *dev_afn, *dev_afs;
extern double *dev_v_ww, *dev_artdivv, *dev_dcudx;
extern double *dev_v_ee, *dev_v_ss, *dev_dcvdy, *dev_v_nn;
extern double *dev_dzudx, *dev_dzvdy, *dev_dp;
extern double *dev_dudx, *dev_dvdy, *dev_rp, *dev_pi;
extern double *dev_res_p, *dev_res_z, *dev_res_c;

extern double **dcdx2, **dcdy2;
extern double *dev_ci, *dev_zi; 

void calcular(int n_imax, int n_itc);
void alocar_globais();
void desalocar_globais();
#endif