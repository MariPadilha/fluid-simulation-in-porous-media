#include "comum.h"
#define idx i*(jmax+1)+j

void mesh(){
    int i, j;

    calcula_x_y();
    calcula_xm_ym();
    calcula_vol();
    calcula_area_das_fases();
    calcula_dx_dy();

//------------bloco L=1 ----------------------------------------
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            if(dev_x[i] <= -4.0 && dev_y[j] >= 0.5){
                dev_flag[idx] = c_i;
            }else{
                dev_flag[idx] = c_f;
            }
        } 
    }
//--------- sem objeto no dominio -------------
    for(i = 2; i <= imax-1; i++){
        for(j = 2; j <= jmax-1; j++){
            if((dev_flag[idx] == c_i && dev_flag[(i-1)*(jmax+1)+j] == c_f)
            || (dev_flag[idx] == c_i && dev_flag[i*(jmax+1)+(j-1)] == c_f)
            || (dev_flag[idx] == c_i && dev_flag[(i+1)*(jmax+1)+j] == c_f)
            || (dev_flag[idx] == c_i && dev_flag[i*(jmax+1)+(j+1)] == c_f)){
                dev_flag[idx] = c_b;
            }
        }
    }

    i = 1;
    for(j = 2; j <= jmax-1; j++){
        if(dev_flag[idx] == c_i && dev_flag[i*(jmax+1)+(j-1)] == c_f
        || dev_flag[idx] == c_i && dev_flag[i*(jmax+1)+(j+1)] == c_f){
            dev_flag[idx] = c_b;
        }
    }

    for(i = 1; i <= imax-1; i++){
        for(j = 1; j <= jmax-1; j++){
            if((dev_flag[idx] == c_f && dev_flag[(i-1)*(jmax+1)+j] == c_b)
            || (dev_flag[idx] == c_f && dev_flag[i*(jmax+1)+(j-1)] == c_b)
            || (dev_flag[idx] == c_f && dev_flag[(i+1)*(jmax+1)+j] == c_b)
            || (dev_flag[idx] == c_f && dev_flag[i*(jmax+1)+(j+1)] == c_b)){
                dev_flag[idx] = c_bs;
            }
        }
    }

//on/off for porosity
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            if(dev_flag[idx] != c_f){
                dev_epsilon1[i*(jmax+1) + j] = porosidade;       
                dev_liga_poros[i*(jmax+1) + j] = 1.0;
            }else{ 
                dev_epsilon1[i*(jmax+1) + j] = 1.0;
                dev_liga_poros[i*(jmax+1) + j] = 0.0;
            }
        } 
    }

//////////////////////////////////////////////////////////////////////////
    #ifdef DEBUG
        FILE *arquivo;
        arquivo = fopen("data/grid_droplet.dat", "w");
        for(i = 1; i <= imax; i++){
            for(j = 1; j <= jmax; j++){
                if(dev_flag[idx] == c_i){ 
                    fprintf(arquivo, "%lf %lf\n", dev_x[i], dev_y[j]);
                }
            }
        }
        fclose(arquivo);

        arquivo = fopen("data/grid_boundary.dat", "w");
        for(i = 1; i <= imax; i++){
            for(j = 1; j <= jmax; j++){
                if(dev_flag[idx] == c_b){ 
                    fprintf(arquivo, "%lf %lf\n", dev_x[i], dev_y[j]);
                }
            }
        }
        fclose(arquivo);

        arquivo = fopen("data/grid_boundary_side.dat", "w");
        for(i = 1; i <= imax; i++){
            for(j = 1; j <= jmax; j++){
                if(dev_flag[idx] == c_bs){ 
                    fprintf(arquivo, "%lf %lf\n", dev_x[i], dev_y[j]);
                }
            }
        }
        fclose(arquivo);
    #endif
}