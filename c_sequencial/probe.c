#include "comum.h"

void probe(double **v, int itc){
    double st_pos = 0.0;

    for(int j=1; j <= jmax; j++){
        if(v[1][j] < 0.0 && v[1][j+1] > 0.0)
            st_pos = 0.5 * (y[j] + y[j+1]);
    }

    FILE *arquivo = fopen("data/probe.dat", "a");
    fprintf(arquivo, "%d %lf\n", itc, st_pos-0.5);
    fclose(arquivo);
}
