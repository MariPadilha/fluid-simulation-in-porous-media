#include "comum.h"

void bcUV(double **um, double **vm){
    int i, j;

    //------ contorno inferior -----------
    for(i = 1; i <= imax; i++){
        vm[i][1] = v_i;
        vm[i][2] = vm[i][1];     
    }

    for(i = 2; i <= (imax); i++){
        um[i][1] = 0.0;    //DARCY
    }

    /*//------ contorno superior ----------
    for(i = 2; i <= (imax); i++){
        vm[i][jmax] = 2.0*vm[i][jmax-2] - vm[i][jmax-1];
        vm[i][(jmax+1)] = 2.0*vm[i][jmax-1] - vm[i][jmax];
    }

    for(i = 2; i <= (imax); i++){
        um[i][jmax] = um[i][jmax-1]; 
    }
    
    //------ contorno esquerdo --------
    for(j = 1; j <= jmax; j++){
        um[2][j] = 0.0; //symmetry
        um[1][j] = 0.0; //symmetry
    }

    for(j = 1; j <= (jmax+1); j++){
        vm[1][j] = vm[2][j]; //darcy
    }

    //------ contorno direito --------
    for(j = 1; j <= jmax; j++){
        um[imax][j] = 0.0;
        um[(imax+1)][j] = 0.0;        
    }

    for(j = 2; j <= (jmax); j++){
        vm[imax][j] =  vm[imax-1][j]; //Darcy
    }*/
}

void bcP(double **pn){
    int i, j;
    
    //--------contorno inferior e superior --------
    for(i = 1; i <= imax; i++){
        pn[i][1] = pn[i][2];
        pn[i][jmax] = pn[i][jmax-1] + 1.0*(pn[i][jmax-1]-pn[i][jmax-2]);
    }
    
    //-------contorno esquerdo e direito -------    
    for(j = 1; j <= jmax; j++){
        pn[1][j] =  pn[2][j];
        pn[imax][j] = pn[imax-1][j];
    }
}

void bcZ(double **zt){
    int i, j;

    double **gradz_x, **gradz_y;
    gradz_x = (double**)malloc(sizeof(double*)*(imax+1));
    gradz_y = (double**)malloc(sizeof(double*)*(imax+1));
    for(i = 0; i < (imax+1); i++){
        gradz_x[i] = (double*)malloc(sizeof(double)*(jmax+1));
        gradz_y[i] = (double*)malloc(sizeof(double)*(jmax+1));  
    }

    //--------contorno inferior e superior --------
    for(i = 1; i <= imax; i++){
        zt[i][1] = tinf; //Zn(i,2)
        zt[i][jmax] = zt[i][jmax-1] + 1.0*(zt[i][jmax-1]-zt[i][jmax-2]);
    }
    
    //-------contorno esquerdo e direito -------    
    for(j = 1; j <= jmax; j++){
        zt[1][j] = zt[2][j];
        zt[imax][j] = zt[imax-1][j];
    }

    i=1;
    for(j = 2; j <= jmax-1; j++){
        gradz_x[i][j] = x[i] * (zt[i+1][j]-zt[i][j]) * dx[i];
        gradz_y[i][j] = y[j] * (zt[i][j+1]-zt[i][j]) * dy[j];
    }

    j=1;
    for(i = 1; i <= imax-1; i++){
        gradz_x[i][j] = x[i] * (zt[i+1][j]-zt[i][j]) * dx[i+1];
        gradz_y[i][j] = y[i] * (zt[i][j+1]-zt[i][j]) * dy[j+1];
    }
        
    for(j = 2; j <= (jmax); j++){
        for(i = 2; i <= (imax); i++){
            gradz_x[i][j] = x[i] * (zt[i][j]-zt[i-1][j]) * dx[i];
            gradz_y[i][j] = y[j] * (zt[i][j]-zt[i][j-1]) * dy[j];
        }
    }

    //desalocando
    for(i = 0; i < (imax+1); i++){
        free(gradz_x[i]);
        free(gradz_y[i]);  
    }
    free(gradz_x);
    free(gradz_y);
}

void bcC(double **cn){
    int i, j;
    double **gradc_x, **gradc_y;
    gradc_x = (double**)malloc(sizeof(double*)*(imax+1));
    gradc_y = (double**)malloc(sizeof(double*)*(imax+1));
    for(i = 0; i < (imax+1); i++){
        gradc_x[i] = (double*)malloc(sizeof(double)*(jmax+1));
        gradc_y[i] = (double*)malloc(sizeof(double)*(jmax+1));  
    }

    //-----contorno inferior e superior --------
    for(i = 1; i <= imax; i++){
        cn[i][1] = 0.0;  //Zn(i,2)
        cn[i][jmax] = cn[i][jmax-1] + 1.0*(cn[i][jmax-1]-cn[i][jmax-2]);
    }
    
    //-------contorno esquerdo e direito -------    
    for(j = 1; j <= jmax; j++){
        cn[1][j] = cn[2][j];
        cn[imax][j] = cn[imax-1][j];
    }

    i=1;
    for(j = 2; j <= jmax-1; j++){
        gradc_x[i][j] = x[i] * (cn[i+1][j]-cn[i][j]) * dx[i];
        gradc_y[i][j] = y[j] * (cn[i][j+1]-cn[i][j]) * dy[j];
    }

    j=1;
    for(i = 1; i <= imax-1; i++){
        gradc_x[i][j] = x[i] * (cn[i+1][j]-cn[i][j]) * dx[i+1];
        gradc_y[i][j] = y[i] * (cn[i][j+1]-cn[i][j]) * dy[j+1];       
    }
        
    for(j = 2; j <= (jmax); j++){
        for(i = 2; i <= (imax); i++){
            gradc_x[i][j] = x[i] * (cn[i][j]-cn[i-1][j]) * dx[i];
            gradc_y[i][j] = y[j] * (cn[i][j]-cn[i][j-1]) * dy[j];
        }
    }

    //desalocando
    for(i = 0; i < (imax+1); i++){
        free(gradc_x[i]);
        free(gradc_y[i]);  
    }
    free(gradc_x);
    free(gradc_y);
}