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

extern double **aww, **aee, **ass;
extern double **ann, **u_ww;
extern double **u_ee, **u_ss, **u_nn, **v_w, **v_ww, **v_e;
extern double **v_ee, **v_s, **v_ss, **v_n, **v_nn;
extern double **dvdydy, **dydudx;
extern double **afw, **afe, **afn, **afs, **dudx, **dvdy, **dzudx, **dzvdy, **dcudx;
extern double **dcvdy, **dcdx2, **dcdy2, **dp, **rp, **pi, **res_p, **rz, **zi;
extern double **rc, **ci, **artdivu, **artdivv, **res_z, **res_c; 

void calcular();
void alocar_globais();
void desalocar_globais();
#endif