#include "comum.h"

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
    FILE *arquivo;

    calcula_x_y();
    calcula_xm_ym();
    calcula_vol_u_v_p();
    calcula_area_das_fases();

    //##########################################################

    for(j = 2; j <= jmax; j++){
        dy[j] = y[j]-y[j-1];
    }

    for(i = 2; i <= imax; i++){
        dx[i] = x[i]-x[i-1];
    }

    printf("dx = %f [mm]\n", maxval(dx, imax+2));
    printf("dy = %f [mm]\n", maxval(dy, jmax+2));

    arquivo = fopen("data/mesh.dat", "w");
    fprintf(arquivo, "%i %i %lf\n", imax, jmax, hvert);
    fclose(arquivo);

    arquivo = fopen("data/grid.dat", "w");
    for(j = 1; j <= jmax; j+=2){ //plotar na direção de i
        for(i = 1; i < imax; i++){
            fprintf(arquivo, "%lf %lf\n", x[i],y[j]);
        } 
        if(j<jmax){
            jj=j+1;
            for(i = imax; i >= 1; i--){
                fprintf(arquivo, "%lf %lf\n", x[i], y[jj]);
            } 
        }
    }

    for(i = 1; i <= imax; i += 2){ //plotar na direção de j
        for(j = jmax; j >= 1; j--){
            fprintf(arquivo, "%lf %lf\n", x[i], y[j]);
        }
        if(i<imax){
            ii=i+1;
            for(j = 1; j <= jmax; j++){
                fprintf(arquivo, "%lf %lf\n", x[ii], y[j]);
            }
        }
    } 
    fclose(arquivo);
 
    arquivo = fopen("data/grid_u.dat", "w");
    for(i = 1; i <= imax+1; i++){
        for(j = 1; j <= jmax; j++){
            fprintf(arquivo, "%lf %lf\n", xm[i], y[j]);
        }
    }
    fclose(arquivo);

    arquivo = fopen("data/grid_v.dat", "w");
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax+1; j++){
            fprintf(arquivo, "%lf %lf\n", x[i], ym[j]);
        }
    }
    fclose(arquivo);

//------------bloco L=1 ----------------------------------------
    for(j = 1; j <= jmax; j++){
        for(i = 1; i <= imax; i++){
            if(x[i] <= -4.0 && y[j] >= 0.5){
                flag[i][j] = c_i;
            }else{
                flag[i][j] = c_f;
            }
        } 
    }
//--------- sem objeto no dominio -------------
    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
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

    for(j = 1; j <= jmax-1; j++){
        for(i = 1; i <= imax-1; i++){
            if((flag[i][j] == c_f && flag[i-1][j] == c_b)
            || (flag[i][j] == c_f && flag[i][j-1] == c_b)
            || (flag[i][j] == c_f && flag[i+1][j] == c_b)
            || (flag[i][j] == c_f && flag[i][j+1] == c_b)){
                flag[i][j] = c_bs;
            }
        }
    }

//on/off for porosity
    for(j = 1; j <= jmax; j++){
        for(i = 1; i <= imax; i++){
            if(flag[i][j] != c_f){
                epsilon1[i][j] = porosidade;       
                liga_poros[i][j] = 1.0;
            }else{ 
                epsilon1[i][j] = 1.0;
                liga_poros[i][j] = 0.0;
            }
        } 
    }

    arquivo = fopen("data/grid_droplet.dat", "w");
    for(i = 1; i <= imax; i++){
        for(int j = 1; j <= jmax; j++){
            if(flag[i][j] == c_i){ 
                fprintf(arquivo, "%lf %lf\n", x[i], y[j]);
            }
        }
    }
    fclose(arquivo);

    arquivo = fopen("data/grid_boundary.dat", "w");
    for(int i = 1; i <= imax; i++){
        for(int j = 1; j <= jmax; j++){
            if(flag[i][j] == c_b){ 
                fprintf(arquivo, "%lf %lf\n", x[i], y[j]);
            }
        }
    }
    fclose(arquivo);

    arquivo = fopen("data/grid_boundary_side.dat", "w");
    for(int i = 1; i <= imax; i++){
        for(int j = 1; j <= jmax; j++){
            if(flag[i][j] == c_bs){ 
                fprintf(arquivo, "%lf %lf\n", x[i], y[j]);
            }
        }
    }
    fclose(arquivo);
}