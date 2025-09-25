#include "comum.h"

void output(double **um, double **vm, double **u, double **v, double **p, double **t, double **c, int k){
    FILE *arquivo;
    
    arquivo = fopen("data/output_variables.dat", "w");
    for(int i = 1; i <= imax; i++){
        for(int j = 1; j <= jmax; j++)
            fprintf(arquivo, "%lf %lf %lf %lf %lf %lf %lf\n", dev_x[i], dev_y[j], u[i][j], v[i][j], p[i][j], t[i][j], c[i][j]);
    }
    fclose(arquivo);

    //--- RESTART/RESTART.dat ---
    arquivo = fopen("data/restart/restartU.dat", "w");
    for(int i = 1; i <= (imax+1); i++){
        for(int j = 1; j <= jmax; j++)
            fprintf(arquivo, "%lf\n", um[i][j]);
    }
    fclose(arquivo);

    arquivo = fopen("data/restart/restartV.dat","w");
    for(int i = 1; i <= imax; i++){
        for(int j = 1; j <= (jmax+1); j++)
            fprintf(arquivo, "%lf\n", vm[i][j]);
    }
    fclose(arquivo);

    arquivo = fopen("data/restart/restartPTC.dat", "w");
    for(int i = 1; i <= imax; i++){
        for(int j = 1; j <= jmax; j++)
            fprintf(arquivo, "%lf %lf %lf\n", p[i][j], t[i][j], c[i][j]);       
    }
    fclose(arquivo);
}