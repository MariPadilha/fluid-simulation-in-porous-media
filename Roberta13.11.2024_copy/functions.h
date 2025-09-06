#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void upwind_Ui(double **um, double **vm, double **p, double **ru, int j);
void upwind_Uj(double **um, double **vm, double **p, double **ru, int i);
void upwind_Vi(double **um, double **vm, double **p, double **rv, double **t, int j);
void upwind_Vj(double **um, double **vm, double **p, double **rv, double **t, int i);
void RESU(double **um, double **vm, double **p, double **ru);
void solve_U(double **um, double **vm, double **um_n, double **um_tau, double **vm_tau, double **um_n_tau, double **p, double *residual_u);
void RESV(double **um, double **vm, double **p, double **t, double **rv);
void solve_V(double **um, double **vm, double **vm_n, double **um_tau, double **vm_tau, double **vm_n_tau, double **p, double **T, double *residual_v);
void solve_P(double **p, double **um_n, double **vm_n, double **pn, double *residual_p);
void RESZ(double **um_n, double **vm_n, double **z, double **rz);
void solve_Z(double **um_n, double **vm_n, double **z, double **z_n_tau, double **z_tau);
void RESC(double **um_n, double **vm_n, double **c, double **rc);
void solve_C(double **um_n, double **vm_n, double **c, double **c_n_tau, double **c_tau);
void bcUV(double **um, double **vm);
void bcP(double **pn);
void bcZ(double **Zt);
void bcC(double **Cn);
void convergence(int itc, double error, double residual_p, double residual_u, double residual_v);
void transient(double **u, double **v, double **p, double **t, double **c, int tr);
void output(double **um, double **vm, double **u, double **v, double **p, double **t, double **c, int k);
void comp_mean(double **u, double **v, double **um, double **vm);
void init();
void IC(double **um, double **vm, double **p, double **t, double **c, double **pn);
void restart(double **um, double **vm, double **p, double **t, double **c);
void restart_dom(double **um, double **vm, double **p, double **t, double **z, double **h);
void mesh();
double max(double a, double b);
void calcula_area_das_fases();
void calcula_x_y();
void calcula_xm_ym();
void calcula_dx_dy();
dim3 grid_1d(int n, int threads);
dim3 grid_2d(int imax, int jmax, int threads_x, int threads_y);

#endif
