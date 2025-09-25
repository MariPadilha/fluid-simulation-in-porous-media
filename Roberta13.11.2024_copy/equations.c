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

//disparar ao mesmo tempo solve u e solve v
