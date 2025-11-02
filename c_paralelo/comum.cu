#include "comum.h"
#include "functions.h"

//#define DEBUG 1

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
double dtau, dt, tempo;

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
double *dev_epsilon1, *dev_liga_poros;
double rad1 = 1.0;     //raio do cilindro

// flags for obstacle interior, boundary, fluid cells, and close to the boundary
int c_i = 2, c_b = 1, c_f = 0, c_bs = 3;     
int *dev_flag;

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
double *dev_res_p;
double *dev_dudx, *dev_dvdy, *dev_rp, *dev_pi;
double *dev_res_z, *dev_res_c; 
double **dcdx2, **dcdy2, *dev_zi, *dev_ci;

void calcular(int n_imax, int n_itc){
    l_c = ao;
    v_c = v_i;
    cf = 1.75 / (150.0 * pow(pow(porosidade,3.0), 0.5));
    hvert = y_up + y_down;
    imax = n_imax;
    dx_c = lhori / (imax-1);
    jmax = (int)((hvert / dx_c) + 1);
    invfr2 = 1.0 / (fr * fr);
    iterations.itc_max = n_itc;
}

dim3 grid_2d(int imax, int jmax){
    return dim3(
        ((imax) + 16 - 1) / 16,
        ((jmax) + 16 - 1) / 16
    );
}

void alocar_globais(){
    //vetores
    cudaMallocManaged((void**)&dev_x, sizeof(double)*(imax+1));
    cudaMallocManaged((void**)&dev_y, sizeof(double)*(jmax+1));
    cudaMallocManaged((void**)&dev_xm, sizeof(double)*(imax+2));
    cudaMallocManaged((void**)&dev_ym, sizeof(double)*(jmax+2));
    cudaMallocManaged((void**)&dev_vol_u, sizeof(double)*(imax+2)*(jmax+1));
    cudaMallocManaged((void**)&dev_vol_v, sizeof(double)*(imax+1)*(jmax+2));
    cudaMallocManaged((void**)&dev_vol_p, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_areau_n, sizeof(double)*(imax+2));
    cudaMallocManaged((void**)&dev_areau_s, sizeof(double)*(imax+2));
    cudaMallocManaged((void**)&dev_areau_e, sizeof(double)*(jmax+1));
    cudaMallocManaged((void**)&dev_areau_w, sizeof(double)*(jmax+1));
    cudaMallocManaged((void**)&dev_areav_n, sizeof(double)*(imax+1));
    cudaMallocManaged((void**)&dev_areav_s, sizeof(double)*(imax+1)); 
    cudaMallocManaged((void**)&dev_areav_e, sizeof(double)*(jmax+2));
    cudaMallocManaged((void**)&dev_areav_w, sizeof(double)*(jmax+2));
    cudaMallocManaged((void**)&dev_dx, sizeof(double)*(imax+2));
    cudaMallocManaged((void**)&dev_dy, sizeof(double)*(jmax+2));

////////////////////////////matrizes linearizadas////////////////////////////////////

    cudaMallocManaged((void**)&dev_epsilon1, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_liga_poros, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_dudx, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_dvdy, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_rp, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_pi, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_res_p, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_res_z, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_res_c, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_zi, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_flag, sizeof(int)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_ci, sizeof(double*)*(imax+1)*(jmax+1));

////////////////////////////////////////////////////////

    dcdx2 = (double**)malloc(sizeof(double*)*(imax+1));
    dcdy2 = (double**)malloc(sizeof(double*)*(imax+1));
    for(int i = 0; i < (imax+1); i++){
        dcdx2[i] = (double*)malloc(sizeof(double)*(jmax+1));
        dcdy2[i] = (double*)malloc(sizeof(double)*(jmax+1));
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
    cudaFree(dev_epsilon1);
    cudaFree(dev_liga_poros);
    cudaFree(dev_dudx);
    cudaFree(dev_dvdy);
    cudaFree(dev_rp);
    cudaFree(dev_pi);
    cudaFree(dev_res_p);
    cudaFree(dev_res_z);
    cudaFree(dev_res_c);
    cudaFree(dev_zi);
    cudaFree(dev_flag);
    cudaFree(dev_ci);

    for(int i = 0; i < (imax+1); i++){
        free(dcdx2[i]);
        free(dcdy2[i]);
    }
    free(dcdx2);
    free(dcdy2);
}