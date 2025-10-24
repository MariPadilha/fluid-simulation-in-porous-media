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
    double etax, etay, sx, sy;
    FILE *arquivo;

    for(i = 1; i <= imax; i++){ //o termo abaixo tira a simetria da malha
        x[i] = (((double)(i) - 1.0) * dx_c) - lhori/(double)(2.0);
    }

    for(j = 1; j <= jmax; j++){
        y[j] = (((double)(j) - 1.0) * dx_c) - y_down;
    }

    //--- APLICA O REFINAMENTO NA MALHA CONSTANTE ---
    for(i = (int)(imax/2+1); i <= imax; i++){
        etax = (x[i] - x[imax]) / (x[(int)(imax/2+1)] - x[imax]);
        sx = px_grid * etax + (1.0 - px_grid)
           * (1.0 - ((tanh(q_grid * (1.0 - etax))) / tanh(q_grid)));
        x[i] = x[imax] - sx * (x[imax] - x[(int)(imax/2+1)]);
    }

    for(i = (int)(imax/2+1); i >= 1; i--){
        etax = (x[i] + x[imax]) / (x[imax/2+1] + x[imax]);
        sx = px_grid * etax + (1.0 - px_grid)
           * (1.0 - (tanh(q_grid * (1.0 - etax)) / tanh(q_grid)));
        x[i] = -x[imax] - sx * (-x[imax] - x[imax/2+1]);
    }

    //aplica a função de stretching para y acima do cilindro
    for(j = ((int)(y_down/dx_c)+1); j <= jmax; j++){
        etay = (y[j] - y[jmax]) / (y[((int)(y_down/dx_c)+1)] - y[jmax]);
        sy = py_grid * etay + (1.0 - py_grid)
           * (1.0 - (tanh((q_grid * (1.0 - etay))) / tanh(q_grid)));
        y[j] = y[jmax] - sy * (y[jmax] - y[(int)(y_down/dx_c)+1]);
    }
    
    for(j = (int)(y_down/dx_c)+1; j >= 1; j--){
        etay = (y[j] + y[jmax]) / (y[(int)(y_down/dx_c)+1] + y[jmax]);
        sy = py_grid * etay + (1.0 - py_grid)
           * (1.0 - (tanh(q_grid * (1.0 - etay)) / tanh(q_grid)));
        y[j] = -y[jmax] - sy * (-y[jmax] - y[(int)((y_down/dx_c)+1)]);
    }

    //--- CALCULA OS PONTOS MEDIOS DA MALHA XM E YM ---
    xm[1] = -(x[2] + x[1]) * (double)(1.0/2.0);
    xm[imax+1] = (x[imax] + x[imax-1]) * (double)(1.0/2.0) + (x[imax] - x[imax-1]);

    for(i = 2; i <= imax; i++){
        xm[i] = (x[i] + x[i-1]) * (double)(1.0/2.0);
    }

    for(i = 2; i <= imax; i++){
        printf("xm[%i]: %.30lf\n", i, xm[i]);
    }
    ym[1] = y[1] - (y[2] - y[1]) * 0.5;
    ym[jmax+1] = (y[jmax] + y[jmax-1]) * 0.5 + (y[jmax] - y[jmax-1]);

    for(j = 2; j <= jmax; j++){
        ym[j] = (y[j] + y[j-1]) * 0.5; 
    }

    //--- CALCULA OS VOLUMES DE U, V E P ---
    for(i = 2; i <= imax; i++){
    for(j = 1 ; j <= jmax; j++){
            vol_u[i][j] = (x[i]-x[i-1]) * (ym[j+1]-ym[j]); 
        }
    }

    for(i = 1; i <= imax; i++){
    for(j = 2; j <= jmax; j++){
            vol_v[i][j] = (y[j]-y[j-1]) * (xm[i+1]-xm[i]); 
        }
    }

    for(i = 1; i <= imax; i++){
    for(j = 2; j <= jmax; j++){
            vol_p[i][j] = (ym[j+1]-ym[j]) * (xm[i+1]-xm[i]); 
        }
    }


    // CALCULA AS AREAS DAS FACES W,E,N,S DE U E V ###################
    for(i = 2; i <= imax; i++){
        areau_n[i] = x[i] - x[i-1];
        areau_s[i] = x[i] - x[i-1];
    }

    for(j = 1; j <= jmax; j++){
        areau_e[j] = ym[j+1] - ym[j];
        areau_w[j] = ym[j+1] - ym[j];
    }

    for(i = 1; i <= imax; i++){
        areav_n[i] = xm[i+1] - xm[i];
        areav_s[i] = xm[i+1] - xm[i];
    }

    for(j = 2; j <= jmax; j++){
        areav_e[j] = y[j] - y[j-1];
        areav_w[j] = y[j] - y[j-1];
    }

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
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            if(x[i] <= -4.0 && y[j] >= 0.5){
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