#include "comum.h"

void init(){
    FILE *arquivo;
    int resultado_arquivo;

    //--- iterations values ---
    arquivo = fopen("input/iterations.dat", "r+");
    resultado_arquivo = fscanf(arquivo, "%d %d %d %d %lf %lf %lf %lf %lf %lf %d",
           &iterations.nc, &iterations.n_tr, &iterations.n_out, &iterations.n_vort,
           &iterations.beta, &iterations.b_art, &iterations.dtau_f, &iterations.final_time,
           &iterations.eps, &iterations.eps_mass, &iterations.start_mode);
    fclose(arquivo);
    printf("itc_max: %d, nc: %d, n_tr: %d, n_out: %d, n_vort: %d\nbeta: %lf, b_art: %lf, dtau_f: %lf, final_time: %lf,\neps: %lf, eps_mass: %lf, start_mode: %d\n",
        iterations.itc_max, iterations.nc, iterations.n_tr, iterations.n_out, iterations.n_vort,
        iterations.beta, iterations.b_art, iterations.dtau_f, iterations.final_time,
        iterations.eps, iterations.eps_mass, iterations.start_mode);

    //--- reference values ---
    arquivo = fopen("input/reference.dat", "r+");
    resultado_arquivo = fscanf(arquivo, "%lf %lf %lf %lf %lf", &ref.tnu, &ref.yf_b, &ref.yo_oo, &ref.ts, &ref.tn_too);
    fclose(arquivo);
    printf("tnu: %lf, yf_b: %lf, yo_oo: %lf, ts: %lf, tn_too: %lf\n", ref.tnu, ref.yf_b, ref.yo_oo, ref.ts, ref.tn_too);

}

void IC(double **um, double **vm, double **p, double **t, double **c, double **pn){ //condicoes iniciais
    int i, j;
    double eta_e;
    
    for(i = 1; i <= imax+1; i++){
        for(j = 1; j <= jmax; j++){
            um[i][j] = 0.0;
        }
    }

    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax+1; j++){
            vm[i][j] = (1.0e-3)*(v_i);
        }
    }

    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            p[i][j] = 1.0;
            //pn[i][j] = 0.0;
            t[i][j] = tinf;
            c[i][j] = 0.0;
        }
    }
    
    for(j = 1; j <= jmax; j++){
        for(i = 1; i <= imax; i++){
            if(flag[i][j] != c_f){
                t[i][j] = temp_cylinder;
                c[i][j] = concentracao_inicial;
            }
        }
    } 

    remove("data/flametip.dat");
    remove("data/error.dat");
}

void restart(double **um, double **vm, double **p, double **t, double **c){
    int i, j;
    FILE *arquivo;
    int resultado_arquivo;

    printf("RESTARTING PROGRAM\n");
    
    arquivo = fopen("data/restart/restartU.dat", "r");
    if (!arquivo){
        perror("Erro ao abrir data/restart/restartU.dat");
        exit(1);
    }
    for(j = 1; j <= jmax; j++){
        for(i = 1; i <= imax+1; i++){
            resultado_arquivo = fscanf(arquivo, "%lf", &um[i][j]);
        }
    }
    fclose(arquivo);

    arquivo = fopen("data/restart/restartV.dat", "r");
    if (!arquivo){
        perror("Erro ao abrir data/restart/restartV.dat");
        exit(1);
    }
    for(j = 1; j <= jmax+1; j++){
        for(i = 1; i <= imax; i++){
            resultado_arquivo = fscanf(arquivo, "%lf", &vm[i][j]);
        }
    }
    fclose(arquivo);

    arquivo = fopen("data/restart/restartPTC.dat", "r");
    if (!arquivo){
        perror("Erro ao abrir data/restart/restartPTC.dat");
        exit(1);
    }
    for(j = 1; j <= jmax; j++){
        for(i = 1; i <= imax; i++){
            resultado_arquivo = fscanf(arquivo, "%lf %lf %lf", &p[i][j], &t[i][j], &c[i][j]);
        }
    }
    fclose(arquivo);
}

void restart_dom(double **um, double **vm, double **p, double **t, double **z, double **h){
    int i, j;
    int rimax = 41;  //em x
    int rjmax = 321; //em y
    double **umr, **vmr, **pres, **zr, **tr, **hr, **h_res;
    FILE *arquivo;
    int resultado_arquivo;

    umr = (double**)malloc(sizeof(double*)*(rimax+2));
    for(i = 0; i < (rimax+2); i++){
        umr[i] = (double*)malloc(sizeof(double)*(rjmax+1));
    }

    vmr = (double**)malloc(sizeof(double*)*(rimax+1));
    for(i = 0; i < (rimax+1); i++){
        vmr[i] = (double*)malloc(sizeof(double)*(rjmax+2));
    }
    
    pres = (double**)malloc(sizeof(double*)*(rimax+1));
    zr = (double**)malloc(sizeof(double*)*(rimax+1));
    tr = (double**)malloc(sizeof(double*)*(rimax+1));
    hr = (double**)malloc(sizeof(double*)*(rimax+1));
    h_res = (double**)malloc(sizeof(double*)*(rimax+1));
    for(i = 0; i < rimax+1; i++){
        pres[i] = (double*)malloc(sizeof(double)*(rjmax+1));
        zr[i] = (double*)malloc(sizeof(double)*(rjmax+1));
        tr[i] = (double*)malloc(sizeof(double)*(rjmax+1));
        hr[i] = (double*)malloc(sizeof(double)*(rjmax+1));
        h_res[i] = (double*)malloc(sizeof(double)*(rjmax+1));
    }        

    printf("RESTARTING PROGRAM\n");

    arquivo = fopen("data/restart/restartU.dat", "r");
    for(j = 1; j <= rjmax; j++){
        for(i = 1; i <= rimax+1; i++){
            resultado_arquivo = fscanf(arquivo, "%lf", &umr[i][j]);
        }
    }
    fclose(arquivo);

    arquivo = fopen("data/restart/restartV.dat", "r");
    for(j = 1; j <= rjmax+1; j++){
        for(i = 1; j <= rimax; j++){
            resultado_arquivo = fscanf(arquivo, "%lf", &vmr[i][j]);
        }
    }
    fclose(arquivo);

    arquivo = fopen("data/restart/restartPTZH.dat", "r");
    for(j = 1; j <= rjmax; j++){
        for(i = 1; i <= rimax; i++){
            resultado_arquivo = fscanf(arquivo, "%lf %lf %lf %lf", &pres[i][j], &tr[i][j], &zr[i][j], &h_res[i][j]);
            hr[i][j] = h_res[i][j] + (((s + 1.0) * lf * tinf / q + 1.0) - h_res[i][j]);
        }
    }
    fclose(arquivo);

    for(j = 1; j <= jmax; j++){
        for(i = 1; i <= imax; i++){
            if(j <= rjmax){
                p[i][j] = pres[i][j];
                t[i][j] = tr[i][j];
                z[i][j] = zr[i][j];
                h[i][j] = hr[i][j];
            }else{
                p[i][j] = p[i][j-1];
                t[i][j] = t[i][j-1];
                z[i][j] = z[i][j-1];
                h[i][j] = h[i][j-1];
            }
        }
    }

    for(j = 1; j <= jmax; j++){
        for(i = 1; i <= imax+1; i++){
            if(j <= rjmax){
                um[i][j] = umr[i][j];
            }else{
                um[i][j] = um[i][j-1];
            }
        }
    }

    for(j = 1; j <= jmax+1; j++){
        for(i = 1; i <= imax; i++){
            if(j <= rjmax){
               vm[i][j] = vmr[i][j];
            }else{
               vm[i][j] = vm[i][j-1];
            }
        }
    }

//desalocando
    for(i = 0; i < (rimax+2); i++){
        free(umr[i]);
    }
    free(umr);

    for(i = 0; i < (rimax+1); i++){
        free(vmr[i]);
    }
    free(vmr);

    for(i = 0; i < rimax+1; i++){
        free(pres[i]);
        free(zr[i]);
        free(tr[i]);
        free(hr[i]);
        free(h_res[i]);
    }  
    free(pres);
    free(zr);
    free(tr);
    free(hr);
    free(h_res);
}