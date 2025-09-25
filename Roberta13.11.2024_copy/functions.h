#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void upwind_Ui(double *dev_um, double *dev_vm, double **p, double *dev_ru, int j);
void upwind_Uj(double *dev_um, double *dev_vm, double **p, double *dev_ru, int i);
void upwind_Vi(double *dev_um, double *dev_vm, double **p, double *dev_rv, double **t, int j);
void upwind_Vj(double *dev_um, double *dev_vm, double **p, double *dev_rv, double **t, int i);
void RESU(double *dev_um, double *dev_vm, double **p, double *dev_ru);
void solve_U(double *dev_um, double *dev_vm, double *dev_um_n, double *dev_um_tau, double *dev_vm_tau, double *dev_um_n_tau, double **p, double *residual_u);
void RESV(double *dev_um, double *dev_vm, double **p, double **t, double *dev_rv);
void solve_V(double *dev_um, double *dev_vm, double *dev_vm_n, double *dev_um_tau, double *dev_vm_tau, double *dev_vm_n_tau, double **p, double **T, double *residual_v);
void solve_P(double **p, double *dev_um_n, double *dev_vm_n, double **pn, double *residual_p);
void RESZ(double *dev_um_n, double *dev_vm_n, double **z, double *dev_rz);
void solve_Z(double *dev_um_n, double *dev_vm_n, double **z, double **z_n_tau, double **z_tau);
void RESC(double *dev_um_n, double *dev_vm_n, double **c, double *dev_rc);
void solve_C(double *dev_um_n, double *dev_vm_n, double **c, double **c_n_tau, double **c_tau);
void bcUV(double *dev_um, double *dev_vm);
void bcP(double **pn);
void bcZ(double **Zt);
void bcC(double **Cn);
void convergence(int itc, double error, double residual_p, double residual_u, double residual_v);
void transient(double **u, double **v, double **p, double **t, double **c, int tr);
void output(double *dev_um, double *dev_vm, double **u, double **v, double **p, double **t, double **c, int k);
void comp_mean(double **u, double **v, double *dev_um, double *dev_vm);
void init();
void IC(double *dev_um, double *dev_vm, double **p, double **t, double **c, double **pn);
void restart(double *dev_um, double *dev_vm, double **p, double **t, double **c);
void restart_dom(double *dev_um, double *dev_vm, double **p, double **t, double **z, double **h);
void mesh();
void calcula_area_das_fases();
void calcula_x_y();
void calcula_xm_ym();
void calcula_dx_dy();
int grid_1d(int tamanho, int threads);
dim3 grid_2d(int imax, int jmax, dim3 threads);

#endif
