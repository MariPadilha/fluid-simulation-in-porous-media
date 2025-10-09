#include "comum.h"

void comp_mean(double **u, double **v, double **um, double **vm){
    int i, j;

    //---CALCULA OS PONTOS MEDIOS DAS VELOCIDADES ---
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            u[i][j] = (um[i+1][j]+um[i][j])*0.50;
        }
    }

    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            v[i][j] = (vm[i][j+1]+vm[i][j])*0.50;
        }
    }
}