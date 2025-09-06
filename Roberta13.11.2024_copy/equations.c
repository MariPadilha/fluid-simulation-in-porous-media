#include "comum.h"

double maior_valor(double **array, int linha, int coluna, int i, int j){
    double max = fabs(array[i][j]);
    for(int i = 3; i < linha; i++){
        for(int j = 2; j < coluna; j++){
            if(array[i][j] > max)
                max = fabs(array[i][j]);
        }
    }
    return max;
}

//--- solve_U ---
void solve_U(double **um, double **vm, double **um_n, double **um_tau, double **vm_tau, double **um_n_tau, double **p, double *residual_u){
    int i, j;
    double **ui, **ru, **res_u;

    ui = (double**)malloc(sizeof(double*)*(imax+2));
    ru = (double**)malloc(sizeof(double*)*(imax+2));
    for(i = 0; i < (imax+2); i++){
        ui[i] = (double*)malloc(sizeof(double)*(jmax+1));
        ru[i] = (double*)malloc(sizeof(double)*(jmax+1));
    }

    res_u = (double**)malloc(sizeof(double*)*(imax+2));
    for(i = 0; i < (imax+2); i++){
        res_u[i] = (double*)malloc(sizeof(double)*(jmax+1));
    }

    RESU(um_tau,vm_tau,p,ru);

    for(j = 2; j <= jmax-1; j++){
        for(i = 3; i <= imax-1; i++){
            res_u[i][j] = ((um[i][j]-um_tau[i][j]) + ru[i][j]*dt) * dtau;
            ui[i][j] = um_tau[i][j] + res_u[i][j];
        }
    }

    bcUV(ui,vm_tau);
    RESU(ui,vm_tau,p,ru);

    for(j = 2; j <= jmax-1; j++){
        for(i = 3; i <= imax-1; i++){
            res_u[i][j] = ((um[i][j]-um_tau[i][j]) + ru[i][j]*dt) * dtau;
            ui[i][j] = 0.75 * um_tau[i][j] + 0.25 * (ui[i][j]+res_u[i][j]);            
        }
    }

    bcUV(ui,vm_tau);
    RESU(ui,vm_tau,p,ru);

    for(j = 2; j <= jmax-1; j++){
        for(i = 3; i <= imax-1; i++){
            res_u[i][j] = ((um[i][j]-um_tau[i][j]) + ru[i][j]*dt) * dtau;
            um_n_tau[i][j] = 1.0 / 3.0 * um_tau[i][j] + 2.0 / 3.0 * (ui[i][j]+res_u[i][j]); 
        }
    }

    bcUV(um_n_tau,vm_tau);

    (*residual_u) =  maior_valor(res_u, imax+2, jmax+1, 3, 2); 

    //desalocando
    for(i = 0; i < (imax+2); i++){
        free(ui[i]);
        free(ru[i]);
    }
    free(ui);
    free(ru);

    for(i = 0; i < (imax+2); i++){
        free(res_u[i]);
    }
    free(res_u);
}

//--- solve_V ---
void solve_V(double **um, double **vm, double **vm_n, double **um_tau, double **vm_tau, double **vm_n_tau, double **p, double **T, double *residual_v){
    int i, j;

    double **vi, **RV, **res_v;
    vi = (double**)malloc(sizeof(double*)*(imax+1));
    RV = (double**)malloc(sizeof(double*)*(imax+1));
    res_v = (double**)malloc(sizeof(double*)*(imax+1));
    for(i = 0; i < imax+1; i++){
        vi[i] = (double*)malloc(sizeof(double)*(jmax+2));
        RV[i] = (double*)malloc(sizeof(double)*(jmax+2));
        res_v[i] = (double*)malloc(sizeof(double)*(jmax+2));
    }

    /*for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){                
        printf("vm_n_tau [%d][%d] = %.14lf\n", i, j, vm_n_tau[i][j]);
        }
    }*/

    RESV(um_tau,vm_tau,p,T,RV);

    for(j = 3; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_v[i][j] = ((vm[i][j]-vm_tau[i][j]) + RV[i][j]*dt) * dtau;
            vi[i][j] = (vm_tau[i][j] + res_v[i][j]);
        }
    }

    bcUV(um_tau,vi);
    RESV(um_tau,vi,p,T,RV);

    for(j = 3; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_v[i][j] = ((vm[i][j]-vm_tau[i][j]) + RV[i][j]*dt) * dtau;
            vi[i][j] = ((double)(3.0/4.0) * vm_tau[i][j] + (double)(1.0/4.0) * (vi[i][j] + res_v[i][j]));
        }   
    }

    bcUV(um_tau,vi);
    RESV(um_tau,vi,p,T,RV);

    for(j = 3; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_v[i][j] =((vm[i][j]-vm_tau[i][j]) +  RV[i][j]*dt) * dtau;
            vm_n_tau[i][j] = (double)(1.0 / 3.0) * vm_tau[i][j] + (double)(2.0 / 3.0) * (vi[i][j] + res_v[i][j]); 
        }
    }

    bcUV(um_tau,vm_n_tau);
    (*residual_v) = maior_valor(res_v, imax+1, jmax+2, 2, 3); 

    //desalocando
    for(i = 0; i < imax+1; i++){
        free(vi[i]);
        free(RV[i]);
        free(res_v[i]);
    }
    free(vi);
    free(RV);
    free(res_v);
}

//--- Solve Continuity Equation ---
void solve_P(double **p, double **um_n, double **vm_n, double **pn, double *residual_p){
    int i, j;

    // RALSTON'S METHOD (Second Order Runge-Kutta) mass conservation
    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            dudx[i][j] = um_n[i+1][j] * areau_e[j] - um_n[i][j] * areau_w[j];
            dvdy[i][j] = vm_n[i][j+1] * areav_n[i] - vm_n[i][j] * areav_s[i];
            rp[i][j] = - (dudx[i][j]+dvdy[i][j]); 
            pi[i][j] = p[i][j] + dtau * rp[i][j] * iterations.beta;
        }
    }

    bcP(pi);

    for(j = 3; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){ 
            pi[i][j] = (double)(3.0/4.0) * p[i][j] + (double)(1.0/4.0) * (pi[i][j] + dtau * rp[i][j] * iterations.beta);
        }
    }

    bcP(pi);

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_p[i][j] = dtau * rp[i][j] * iterations.beta;
            pn[i][j] = (double)(1.0 / 3.0) * p[i][j] + (double)(2.0 / 3.0) * (pi[i][j] + res_p[i][j]);
        }
    }


    bcP(pn);
            //for(i = 1; i <= imax; i++){
            //    for(j = 1; j <= jmax; j++){
            //        printf("pn [%d][%d] = %lf\n", i, j, pn[i][j]);
            //    }
            //}
    (*residual_p) = maior_valor(res_p, imax+1, jmax+1, 1, 1);
}

//--- Solve Mixture Fraction ---
void solve_Z(double **um_n, double **vm_n, double **z, double **z_n_tau, double **z_tau){
    int i, j;

    // RALSTON'S METHOD (Second Order Runge-Kutta)
    RESZ(um_n,vm_n, z_tau, rz);

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){    
            res_z[i][j] = ((z[i][j]-z_tau[i][j]) + rz[i][j]*dt) * dtau; 
            zi[i][j] = z_tau[i][j] + res_z[i][j];             
        }
    }

    bcZ(zi);
    RESZ(um_n,vm_n,zi,rz);

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){    
            res_z[i][j] = ((z[i][j]-z_tau[i][j]) + rz[i][j]*dt) * dtau;
            zi[i][j] = 0.75 * z_tau[i][j] + 0.25 * (zi[i][j]+res_z[i][j]);
        }
    }
    
    bcZ(zi);
    RESZ(um_n,vm_n,zi,rz);

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_z[i][j] = ((z[i][j]-z_tau[i][j]) + rz[i][j]*dt) * dtau; 
            z_n_tau[i][j] = 1.0 / 3.0 * z_tau[i][j] + 2.0 / 3.0 * (zi[i][j]+res_z[i][j]);
        }
    }

    bcZ(z_n_tau);
}

//--- solve_C - concentration ---
void solve_C(double **um_n, double **vm_n, double **c, double **c_n_tau, double **c_tau){
    int i, j;

    // RALSTON'S METHOD (Second Order Runge-Kutta)
    RESC(um_n, vm_n, c_tau, rc);

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_c[i][j] = ((c[i][j]-c_tau[i][j]) + rc[i][j]*dt) * dtau; 
            ci[i][j] = c_tau[i][j] + res_c[i][j];
        }
    }

    bcC(ci);
    RESC(um_n,vm_n,ci,rc);
    
    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_c[i][j] = ((c[i][j]-c_tau[i][j]) + rc[i][j]*dt) * dtau; 
            ci[i][j] = 0.75 * c_tau[i][j] + 0.25 * (ci[i][j] + res_c[i][j]);   
        }
    }

    bcC(ci);
    RESC(um_n, vm_n, ci, rc);

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_c[i][j] = ((c[i][j]-c_tau[i][j]) + rc[i][j]*dt) * dtau; 
            c_n_tau[i][j] = 1.0 / 3.0 * c_tau[i][j] + 2.0 / 3.0 * (ci[i][j] + res_c[i][j]);            
        }
    }

    bcC(c_n_tau);
}