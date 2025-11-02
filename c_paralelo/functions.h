#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include "comum.h"

void calcula_vol();
void upwind_U_fusion(double *dev_um, double *dev_vm, double *dev_p, double *dev_ru, int i, int j);
void upwind_U_pair(double *dev_um, double *dev_vm, double *dev_p, double *dev_ru);
void upwind_V_pair(double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_rv);
void upwind_Uj(double *dev_um, double *dev_vm, double *dev_p, double *dev_ru, int i);
void upwind_Ui(double *dev_um, double *dev_vm, double *dev_p, double *dev_ru, int j);
void upwind_Vi(double *dev_um, double *dev_vm, double *dev_p, double *dev_rv, double *dev_t, int j);
void upwind_Vj(double *dev_um, double *dev_vm, double *dev_p, double *dev_rv, double *dev_t, int i);
void RESU(double *dev_um, double *dev_vm, double *dev_p, double *dev_ru);
void RESV(double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_rv);
void RESZ(double *dev_um_n, double *dev_vm_n, double *dev_z, double *dev_rz);
void RESC(double *dev_um_n, double *dev_vm_n, double *dev_c, double *dev_rc);
void solve_Z(double *dev_um_n, double *dev_vm_n, double *dev_z, double *dev_z_n_tau, double *dev_z_tau);
void solve_C(double *dev_um_n, double *dev_vm_n, double *dev_c, double *dev_c_n_tau, double *dev_c_tau);
double solve_U(double *dev_um, double *dev_vm, double *dev_um_n, double *dev_um_tau, double *dev_vm_tau, double *dev_um_n_tau, double *dev_p);
double solve_V(double *dev_um, double *dev_vm, double *dev_vm_n, double *dev_um_tau, double *dev_vm_tau, double *dev_vm_n_tau, double *dev_p, double *dev_t);
double solve_P(double *dev_p, double *dev_um_n, double *dev_vm_n, double *dev_pn);
void bcUV(double *dev_um, double *dev_vm);
void bcP(double *dev_pn);
void bcZ(double *dev_zt);
void bcC(double *dev_cn);
void convergence(int itc, double error, double residual_p, double residual_u, double residual_v);
void transient(double *dev_u, double *dev_v, double *dev_p, double *dev_t, double *dev_c, int tr);
void output(double *dev_um, double *dev_vm, double *dev_u, double *dev_v, double *dev_p, double *dev_t, double *dev_c, int k);
void comp_mean(double *dev_u, double *dev_v, double *dev_um, double *dev_vm);
void init();
void IC(double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_c, double *dev_pn);
void restart(double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_c);
void restart_dom(double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_z, double *dev_h);
void mesh();
void calcula_area_das_fases();
void calcula_x_y();
void calcula_xm_ym();
void calcula_dx_dy();
int grid_1d(int tamanho, int threads);
double max_reduce(double* dev_matriz_linearizada, int imax, int jmax);
#endif
