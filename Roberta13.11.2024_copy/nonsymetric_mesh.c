#include "comum.h"
//#define OUTPUT 1 //Para debug

double maxval(double *array, int n){
    double max = array[0];
    for(int i = 1; i < n; i++){
        if (array[i] > max)
            max = array[i];
    }
    return max;
}

void mesh(){
    int i, j, ii=0, jj=0;

    calcula_x_y();
    calcula_xm_ym();
    calcula_area_das_fases();
    calcula_dy_dx();

//------------bloco L=1 ----------------------------------------
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            if(dev_x[i] <= -4.0 && dev_y[j] >= 0.5){
                flag[i][j] = c_i;
            }else{
                flag[i][j] = c_f;
            }
        } 
    }
//--------- sem objeto no dominio -------------
    for(i = 2; i <= imax-1; i++){
        for(j = 2; j <= jmax-1; j++){
            if((flag[i][j] == c_i && flag[i-1][j] == c_f)
            || (flag[i][j] == c_i && flag[i][j-1] == c_f)
            || (flag[i][j] == c_i && flag[i+1][j] == c_f)
            || (flag[i][j] == c_i && flag[i][j+1] == c_f)){
                flag[i][j] = c_b;
            }
        }
    }

    i = 1;
    for(j = 2; j <= jmax-1; j++){
        if(flag[i][j] == c_i && flag[i][j-1] == c_f
        || flag[i][j] == c_i && flag[i][j+1] == c_f){
            flag[i][j] = c_b;
        }
    }

    for(i = 1; i <= imax-1; i++){
        for(j = 1; j <= jmax-1; j++){
            if((flag[i][j] == c_f && flag[i-1][j] == c_b)
            || (flag[i][j] == c_f && flag[i][j-1] == c_b)
            || (flag[i][j] == c_f && flag[i+1][j] == c_b)
            || (flag[i][j] == c_f && flag[i][j+1] == c_b)){
                flag[i][j] = c_bs;
            }
        }
    }

//on/off for porosity
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            if(flag[i][j] != c_f){
                epsilon1[i][j] = porosidade;       
                liga_poros[i][j] = 1.0;
            }else{ 
                epsilon1[i][j] = 1.0;
                liga_poros[i][j] = 0.0;
            }
        } 
    }

//////////////////////////////////////////////////////////////////////////
//debug
    #ifdef OUTPUT
        FILE *arquivo;
        arquivo = fopen("data/grid_droplet.dat", "w");
        for(i = 1; i <= imax; i++){
            for(int j = 1; j <= jmax; j++){
                if(flag[i][j] == c_i){ 
                    fprintf(arquivo, "%lf %lf\n", dev_x[i], dev_y[j]);
                }
            }
        }
        fclose(arquivo);

        arquivo = fopen("data/grid_boundary.dat", "w");
        for(int i = 1; i <= imax; i++){
            for(int j = 1; j <= jmax; j++){
                if(flag[i][j] == c_b){ 
                    fprintf(arquivo, "%lf %lf\n", dev_x[i], dev_y[j]);
                }
            }
        }
        fclose(arquivo);

        arquivo = fopen("data/grid_boundary_side.dat", "w");
        for(int i = 1; i <= imax; i++){
            for(int j = 1; j <= jmax; j++){
                if(flag[i][j] == c_bs){ 
                    fprintf(arquivo, "%lf %lf\n", dev_x[i], dev_y[j]);
                }
            }
        }
        fclose(arquivo);
    #endif
}