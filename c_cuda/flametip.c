#include "comum.h"

void flametip(double **Z, int itc){
    int i = 1;
    double yf;

    for(int j = (int)(y_down/dx_c); j <= jmax; j++){
        if(Z[i][j] > 1.0){
            yf = ((1.0-Z[i][j-1])*(dev_y[j]- dev_y[j-1] )) / (Z[i][j]-Z[i][j-1]) + dev_y[j-1];
            break;
        }
    }

    FILE *arquivo = fopen("data/flametip.dat", "a");
        
    if(yf < y_up)
        fprintf(arquivo, "%d %lf\n", itc, yf);
    else
        fprintf(arquivo, "%d %lf\n", itc, y_up);        
    fclose(arquivo);
}