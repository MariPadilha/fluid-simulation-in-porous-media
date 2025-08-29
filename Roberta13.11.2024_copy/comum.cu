#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "comum.h"
#include "functions.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

struct Iterations iterations;
struct Ref ref;
//Parâmetros de iteração e controle
//int itc_max;        //numero de iterações

/*// Frequência dos outputs:
int nc;                //erros
int n_tr;              //plota parte transiente
int n_out;             //salva resultados preliminares
int n_vort;            //salva dados do vortice
double beta;           //parâmetro de compressibilidade        
double b_art;          //coeficiente da dissipação artificial
double dtau_f;         //fator de correcao para calc de dt        
double final_time;     //tempo máximo de duração do tempo da simulação
double eps, eps_mass;  //criterio de convergencia  
*/
int restart_mode;      //tipo de start, se eh CI ou solucao anterior

// Passo de tempo
double dtau, dt;
double time;

// Parâmetros físicos e geométricos
double porosidade = 0.5;                               //Lido em main, mesh, nonsymetric_mesh
double darcy_number = 1.0e-2;                          //Lido em equations
double cf;                                  //1.75 / pow((150.0 * pow(porosidade,3.0)), 0.5);
double temp_cylinder = 1.0, concentracao_inicial = 1.0;

    
// Parâmetros de refinamento P e Q
// colocar P=1 para a malha uniforme
double px_grid = 1.0, py_grid = 1.60, q_grid = 1.0; 


// Tamanho do domínio e da malha
double lhori = 12.0;            //largura total
double y_up = 15.0;             //altura em y+
double y_down = 5.0;            //altura em y-
double hvert;            //y_up + y_down;   //altura total
int imax;                //numero de pontos da malha em x
double dx_c;         //Lhori / (imax-1); //dita o tamanho de dy e dx
int jmax;  //(Hvert / dx_c) + 1;  //numero de pontos da malha em y


// Main data structures for control of the mesh
double *dev_x, *dev_y;       //malha principal
double *dev_xm, *dev_ym; //malha deslocada
double *dev_vol_u;    //volume de controle de u
double *dev_vol_v;     //volume de controle de v
double *dev_vol_p;       //volume de controle de p

//Áreas para u e v
double *dev_areau_n;  //area n de u
double *dev_areau_s;  //area s de u
double *dev_areau_e;  //area e de u
double *dev_areau_w;  //area w de u
double *dev_areav_n;  //area n de v
double *dev_areav_s;  //area s de v
double *dev_areav_e;  //area e de v
double *dev_areav_w;  //area w de v

// Variáveis auxiliares        
double *dev_dx, *dev_dy; 

double **epsilon1, **liga_poros; 

double rad1 = 1.0;     //raio do cilindro

// flags for obstacle interior, boundary, fluid cells, and close to the boundary
int c_i = 2, c_b = 1, c_f = 0, c_bs = 3;     
int **flag;

// Constantes físicas e parâmetros de fluidos
double g = 9.80665;               //gravitational constant [m/s^2]
double ao  = 1.0e-3;              //initial radius [m]
double v_i = 0.5;
double l_c;
double v_c;    
double fr = 1.0;
double invfr2;
double s;
double lf = 1.0;
double lo = 1.0;

/*//NAMELIST ref
double Tnu;    // for Methane , 3.51d0 for n-Heptane
double YF_b; 
double YO_oo; 
double Ts;     // Tb  = boiling temperature [k]
double TnToo;
*/

//Compute in main, after init. Depends of Ts
double too;                //dimen ambient temp [k]
double tsup;
double tinf;
double q_dim = 5.015e7;    //used only to q calculation  ! q = combustion heat [J/kg] for Methane CH4
double q;                  //Compute in initial, equations, boundary, Depends of properties


//Parameters computed in subroutine properties 
double cp_tot;          // [J/kgK]        // Affects q
double rho_tot;         // [kg/m^3]       // No affect
double k_tot;           // [W/mK]         // No affect
double nu_tot;          // [m^2/s]        // Affects Pr
double alpha_tot;       // [m^2/s]        // Affects Pr
double re = 1.0;
double pr;              // Affects Pe
double pe = 1.0;
double sc = 1.0;

//Temporary variables
double **fw, **fe, **fs, **fn, **df, **aw, **aww, **ae, **aee, **as, **ass, **an;
double **ann, **ap, **dn, **ds, **de, **dw, **u_w, **u_ww, **u_e, **u_ee, **u_s;
double **u_ss, **u_n, **u_nn, **u_p, **v_w, **v_ww, **v_e, **v_ee, **v_s, **v_ss;
double **v_n, **v_nn, **v_p, **q_art, **dudxdx, **dvdydy, **dxdvdy, **dydudx, **afw;
double **afe, **afn, **afs, **dudx, **dvdy, **dzudx, **dzvdy, **dcudx, **dcvdy;
double **dcdx2, **dcdy2, **dp, **rp, **pi, **res_p, **rz, **zi, **rc, **ci;
double **artdivu, **artdivv, **res_z, **res_c; 

void calcular(){
    l_c = ao;
    v_c = v_i;
    cf = 1.75 / (150.0 * pow(pow(porosidade,3.0), 0.5));                                  //1.75 / pow((150.0 * pow(porosidade,3.0)), 0.5);
    hvert = y_up + y_down;
    imax = 10;
    dx_c = lhori / (imax-1);
    jmax = (int)((hvert / dx_c) + 1);
    invfr2 = 1.0 / (fr * fr);
}

void alocar_globais(){
    cudaMallocManaged((void**)&dev_x, sizeof(double)*(imax+1));
    cudaMallocManaged((void**)&dev_y, sizeof(double)*(jmax+1));
    cudaMalloc((void**)&dev_xm, sizeof(double)*(imax+2));
    cudaMalloc((void**)&dev_ym, sizeof(double)*(jmax+2));
    cudaMalloc((void**)&dev_vol_u, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_vol_v, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_vol_p, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_areau_n, sizeof(double)*(imax+2));
    cudaMalloc((void**)&dev_areau_s, sizeof(double)*(imax+2));
    cudaMalloc((void**)&dev_areau_e, sizeof(double)*(jmax+1));
    cudaMalloc((void**)&dev_areau_w, sizeof(double)*(jmax+1));
    cudaMalloc((void**)&dev_areav_n, sizeof(double)*(imax+1));
    cudaMalloc((void**)&dev_areav_s, sizeof(double)*(imax+1)); 
    cudaMalloc((void**)&dev_areav_e, sizeof(double)*(jmax+2));
    cudaMalloc((void**)&dev_areav_w, sizeof(double)*(jmax+2));
    cudaMalloc((void**)&dev_dx, sizeof(double)*(imax+2));
    cudaMalloc((void**)&dev_dy, sizeof(double)*(jmax+2));

    epsilon1 = (double**)malloc(sizeof(double*)*(imax+1));
    liga_poros = (double**)malloc(sizeof(double*)*(imax+1));
    for(int i = 0; i < (imax+1); i++){
        epsilon1[i] = (double*)malloc(sizeof(double)*(jmax+1));
        liga_poros[i] = (double*)malloc(sizeof(double)*(jmax+1));
    }

    flag = (int**)malloc(sizeof(int*)*(imax+1));
    for(int i = 0; i < (imax+1); i++){
        flag[i] = (int*)malloc(sizeof(int)*(jmax+1));
    }

    fw = (double**)malloc(sizeof(double*)*(imax+2));
    fe = (double**)malloc(sizeof(double*)*(imax+2));
    fs = (double**)malloc(sizeof(double*)*(imax+2));
    fn = (double**)malloc(sizeof(double*)*(imax+2));
    df = (double**)malloc(sizeof(double*)*(imax+2));
    aw = (double**)malloc(sizeof(double*)*(imax+2));
    aww = (double**)malloc(sizeof(double*)*(imax+2));
    ae = (double**)malloc(sizeof(double*)*(imax+2));
    aee = (double**)malloc(sizeof(double*)*(imax+2));
    as = (double**)malloc(sizeof(double*)*(imax+2));
    ass = (double**)malloc(sizeof(double*)*(imax+2));
    an = (double**)malloc(sizeof(double*)*(imax+2));
    ann = (double**)malloc(sizeof(double*)*(imax+2));
    ap = (double**)malloc(sizeof(double*)*(imax+2));
    dn = (double**)malloc(sizeof(double*)*(imax+2));
    ds = (double**)malloc(sizeof(double*)*(imax+2));
    de = (double**)malloc(sizeof(double*)*(imax+2));
    dw = (double**)malloc(sizeof(double*)*(imax+2));
    u_w = (double**)malloc(sizeof(double*)*(imax+2));
    u_ww = (double**)malloc(sizeof(double*)*(imax+2));
    u_e = (double**)malloc(sizeof(double*)*(imax+2));
    u_ee = (double**)malloc(sizeof(double*)*(imax+2));
    u_s = (double**)malloc(sizeof(double*)*(imax+2));
    u_ss = (double**)malloc(sizeof(double*)*(imax+2));
    u_n = (double**)malloc(sizeof(double*)*(imax+2));
    u_nn = (double**)malloc(sizeof(double*)*(imax+2));
    u_p = (double**)malloc(sizeof(double*)*(imax+2));
    v_w = (double**)malloc(sizeof(double*)*(imax+2));
    v_ww = (double**)malloc(sizeof(double*)*(imax+2));
    v_e = (double**)malloc(sizeof(double*)*(imax+2));
    v_ee = (double**)malloc(sizeof(double*)*(imax+2));
    v_s = (double**)malloc(sizeof(double*)*(imax+2));
    v_ss = (double**)malloc(sizeof(double*)*(imax+2));
    v_n = (double**)malloc(sizeof(double*)*(imax+2));
    v_nn = (double**)malloc(sizeof(double*)*(imax+2));
    v_p = (double**)malloc(sizeof(double*)*(imax+2));
    q_art = (double**)malloc(sizeof(double*)*(imax+2));
    dudxdx = (double**)malloc(sizeof(double*)*(imax+2));
    dvdydy = (double**)malloc(sizeof(double*)*(imax+2));
    dxdvdy = (double**)malloc(sizeof(double*)*(imax+2));
    dydudx = (double**)malloc(sizeof(double*)*(imax+2));
    for(int i = 0; i < (imax+2); i++){
        fw[i] = (double*)malloc(sizeof(double)*(jmax+1));
        fe[i] = (double*)malloc(sizeof(double)*(jmax+1));
        fs[i] = (double*)malloc(sizeof(double)*(jmax+1));
        fn[i] = (double*)malloc(sizeof(double)*(jmax+1));
        df[i] = (double*)malloc(sizeof(double)*(jmax+1));
        aw[i] = (double*)malloc(sizeof(double)*(jmax+1));
        aww[i] = (double*)malloc(sizeof(double)*(jmax+1));
        ae[i] = (double*)malloc(sizeof(double)*(jmax+1));
        aee[i] = (double*)malloc(sizeof(double)*(jmax+1));
        as[i] = (double*)malloc(sizeof(double)*(jmax+1));
        ass[i] = (double*)malloc(sizeof(double)*(jmax+1));
        an[i] = (double*)malloc(sizeof(double)*(jmax+1));
        ann[i] = (double*)malloc(sizeof(double)*(jmax+1));
        ap[i] = (double*)malloc(sizeof(double)*(jmax+1));
        dn[i] = (double*)malloc(sizeof(double)*(jmax+1));
        ds[i] = (double*)malloc(sizeof(double)*(jmax+1));
        de[i] = (double*)malloc(sizeof(double)*(jmax+1));
        dw[i] = (double*)malloc(sizeof(double)*(jmax+1));
        u_w[i] = (double*)malloc(sizeof(double)*(jmax+1));
        u_ww[i] = (double*)malloc(sizeof(double)*(jmax+1));
        u_e[i] = (double*)malloc(sizeof(double)*(jmax+1));
        u_ee[i] = (double*)malloc(sizeof(double)*(jmax+1));
        u_s[i] = (double*)malloc(sizeof(double)*(jmax+1));
        u_ss[i] = (double*)malloc(sizeof(double)*(jmax+1));
        u_n[i] = (double*)malloc(sizeof(double)*(jmax+1));
        u_nn[i] = (double*)malloc(sizeof(double)*(jmax+1));
        u_p[i] = (double*)malloc(sizeof(double)*(jmax+1));
        v_w[i] = (double*)malloc(sizeof(double)*(jmax+1));
        v_ww[i] = (double*)malloc(sizeof(double)*(jmax+1));
        v_e[i] = (double*)malloc(sizeof(double)*(jmax+1));
        v_ee[i] = (double*)malloc(sizeof(double)*(jmax+1));
        v_s[i] = (double*)malloc(sizeof(double)*(jmax+1));
        v_ss[i] = (double*)malloc(sizeof(double)*(jmax+1));
        v_n[i] = (double*)malloc(sizeof(double)*(jmax+1));
        v_nn[i] = (double*)malloc(sizeof(double)*(jmax+1));
        v_p[i] = (double*)malloc(sizeof(double)*(jmax+1));
        q_art[i] = (double*)malloc(sizeof(double)*(jmax+1));
        dudxdx[i] = (double*)malloc(sizeof(double)*(jmax+1));
        dvdydy[i] = (double*)malloc(sizeof(double)*(jmax+1));
        dxdvdy[i] = (double*)malloc(sizeof(double)*(jmax+1));
        dydudx[i] = (double*)malloc(sizeof(double)*(jmax+1));
    }


    afw = (double**)malloc(sizeof(double*)*(imax+1));                                  
    afe = (double**)malloc(sizeof(double*)*(imax+1));                                  
    afn = (double**)malloc(sizeof(double*)*(imax+1));                                  
    afs = (double**)malloc(sizeof(double*)*(imax+1));
    dudx = (double**)malloc(sizeof(double*)*(imax+1));
    dvdy = (double**)malloc(sizeof(double*)*(imax+1));
    dzudx = (double**)malloc(sizeof(double*)*(imax+1));
    dzvdy = (double**)malloc(sizeof(double*)*(imax+1));
    dcudx = (double**)malloc(sizeof(double*)*(imax+1));
    dcvdy = (double**)malloc(sizeof(double*)*(imax+1));
    dcdx2 = (double**)malloc(sizeof(double*)*(imax+1));
    dcdy2 = (double**)malloc(sizeof(double*)*(imax+1));
    dp = (double**)malloc(sizeof(double*)*(imax+1));
    rp = (double**)malloc(sizeof(double*)*(imax+1));
    pi = (double**)malloc(sizeof(double*)*(imax+1));
    res_p = (double**)malloc(sizeof(double*)*(imax+1));
    rz = (double**)malloc(sizeof(double*)*(imax+1));
    zi = (double**)malloc(sizeof(double*)*(imax+1));
    rc = (double**)malloc(sizeof(double*)*(imax+1));
    ci = (double**)malloc(sizeof(double*)*(imax+1));
    artdivu = (double**)malloc(sizeof(double*)*(imax+1));
    artdivv = (double**)malloc(sizeof(double*)*(imax+1));
    res_z = (double**)malloc(sizeof(double*)*(imax+1));
    res_c = (double**)malloc(sizeof(double*)*(imax+1));
    for(int i = 0; i < (imax+1); i++){
        afw[i] = (double*)malloc(sizeof(double)*(jmax+1));                                 
        afe[i] = (double*)malloc(sizeof(double)*(jmax+1));                                  
        afn[i] = (double*)malloc(sizeof(double)*(jmax+1));                                  
        afs[i] = (double*)malloc(sizeof(double)*(jmax+1)); 
        dudx[i] = (double*)malloc(sizeof(double)*(jmax+1));
        dvdy[i] = (double*)malloc(sizeof(double)*(jmax+1));
        dzudx[i] = (double*)malloc(sizeof(double)*(jmax+1));
        dzvdy[i] = (double*)malloc(sizeof(double)*(jmax+1));
        dcudx[i] = (double*)malloc(sizeof(double)*(jmax+1));  
        dcvdy[i] = (double*)malloc(sizeof(double)*(jmax+1));
        dcdx2[i] = (double*)malloc(sizeof(double)*(jmax+1));
        dcdy2[i] = (double*)malloc(sizeof(double)*(jmax+1));
        dp[i] = (double*)malloc(sizeof(double)*(jmax+1));
        rp[i] = (double*)malloc(sizeof(double)*(jmax+1));
        pi[i] = (double*)malloc(sizeof(double)*(jmax+1));
        res_p[i] = (double*)malloc(sizeof(double)*(jmax+1));
        rz[i] = (double*)malloc(sizeof(double)*(jmax+1));
        zi[i] = (double*)malloc(sizeof(double)*(jmax+1));
        rc[i] = (double*)malloc(sizeof(double)*(jmax+1));
        ci[i] = (double*)malloc(sizeof(double)*(jmax+1));
        artdivu[i] = (double*)malloc(sizeof(double)*(jmax+1));
        artdivv[i] = (double*)malloc(sizeof(double)*(jmax+1));
        res_z[i] = (double*)malloc(sizeof(double)*(jmax+1));
        res_c[i] = (double*)malloc(sizeof(double)*(jmax+1));
    }
}

void desalocar_globais(){
    cudaFree(dev_x);
    cudaFree(dev_y);
    cudaFree(dev_dx);
    cudaFree(dev_dy);
    cudaFree(dev_vol_u);
    cudaFree(dev_vol_v);
    cudaFree(dev_vol_p);
    cudaFree(dev_xm);
    cudaFree(dev_ym);
    cudaFree(dev_areau_n);
    cudaFree(dev_areau_s);
    cudaFree(dev_areau_e);
    cudaFree(dev_areau_w);
    cudaFree(dev_areav_n);
    cudaFree(dev_areav_s);
    cudaFree(dev_areav_e);
    cudaFree(dev_areav_w);

    for(int i = 0; i < (imax+2); i++){
        free(fw[i]);
        free(fe[i]);
        free(fs[i]);
        free(fn[i]);
        free(df[i]);
        free(aw[i]);
        free(aww[i]);
        free(ae[i]);
        free(aee[i]);
        free(as[i]);
        free(ass[i]);
        free(an[i]);
        free(ann[i]);
        free(ap[i]);
        free(dn[i]);
        free(ds[i]);
        free(de[i]);
        free(dw[i]);
        free(u_w[i]);
        free(u_ww[i]);
        free(u_e[i]);
        free(u_ee[i]);
        free(u_s[i]);
        free(u_ss[i]);
        free(u_n[i]);
        free(u_nn[i]);
        free(u_p[i]);
        free(v_w[i]);
        free(v_ww[i]);
        free(v_e[i]);
        free(v_ee[i]);
        free(v_s[i]);
        free(v_ss[i]);
        free(v_n[i]);
        free(v_nn[i]);
        free(v_p[i]);
        free(q_art[i]);
        free(dudxdx[i]);
        free(dvdydy[i]);
        free(dxdvdy[i]);
        free(dydudx[i]);
    }
    free(fw);
    free(fe);
    free(fs);
    free(fn);
    free(df);
    free(aw);
    free(aww);
    free(ae);
    free(aee);
    free(as);
    free(ass);
    free(an);
    free(ann);
    free(ap);
    free(dn);
    free(ds);
    free(de);
    free(dw);
    free(u_w);
    free(u_ww);
    free(u_e);
    free(u_ee);
    free(u_s);
    free(u_ss);
    free(u_n);
    free(u_nn);
    free(u_p);
    free(v_w);
    free(v_ww);
    free(v_e);
    free(v_ee);
    free(v_s);
    free(v_ss);
    free(v_n);
    free(v_nn);
    free(v_p);
    free(q_art);
    free(dudxdx);
    free(dvdydy);
    free(dxdvdy);
    free(dydudx);

    for(int i = 0; i < (imax+1); i++){
        free(epsilon1[i]);
        free(liga_poros[i]);
        free(flag[i]);
        free(afw[i]);                         
        free(afe[i]);                         
        free(afn[i]); 
        free(afs[i]); 
        free( dudx[i]);
        free(dvdy[i]);
        free(dzudx[i]);
        free(dzvdy[i]);
        free(dcudx[i]);  
        free(dcvdy[i]);
        free(dcdx2[i]);
        free(dcdy2[i]);
        free(dp[i]);
        free(rp[i]);
        free(pi[i]);
        free(res_p[i]);
        free(rz[i]);
        free(zi[i]);
        free(rc[i]);
        free(ci[i]);
        free(artdivu[i]);
        free(artdivv[i]);
        free(res_z[i]);
        free(res_c[i]);
    }
    free(epsilon1);
    free(liga_poros);
    free(flag);
    free(afw);                             
    free(afe);                             
    free(afn);                             
    free(afs);
    free(dudx);
    free(dvdy);
    free(dzudx);
    free(dzvdy);
    free(dcudx);
    free(dcvdy);
    free(dcdx2);
    free(dcdy2);
    free(dp);
    free(rp);
    free(pi);
    free(res_p);
    free(rz);
    free(zi);
    free(rc);
    free(ci);
    free(artdivu);
    free(artdivv);
    free(res_z);
    free(res_c);
}