#include "comum.h"
#define idx i*(jmax+1)+j 

void transient(double *dev_u, double *dev_v, double *dev_p, double *dev_t, double *dev_c, int tr){    
    char filename2[100];
    char filepath[256];
    FILE *arquivo;

    int num = (tr / iterations.n_tr) + 100, i, j;
    sprintf(filename2, "%d", num);
    
    snprintf(filepath, sizeof(filepath), "transient/data/u%s.dat", filename2);
    arquivo = fopen(filepath, "w");
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            fprintf(arquivo, "%lf ", dev_u[idx]);
        }
        fprintf(arquivo, "\n");
    }
    fclose(arquivo);

    snprintf(filepath, sizeof(filepath), "transient/data/v%s.dat", filename2);
    arquivo = fopen(filepath, "w");
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            fprintf(arquivo, "%lf ", dev_v[idx]);
        }
        fprintf(arquivo, "\n");
    }
    fclose(arquivo);
    
    snprintf(filepath, sizeof(filepath), "transient/data/P%s.dat", filename2);
    arquivo = fopen(filepath, "w");
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            fprintf(arquivo, "%lf ", dev_p[idx]); // espaço entre números     
        }
        fprintf(arquivo, "\n");
    }
    fclose(arquivo);
    
    snprintf(filepath, sizeof(filepath), "transient/data/T%s.dat", filename2);
    arquivo = fopen(filepath, "w");
    for(i = 1; i <=imax; i++){
        for(j = 1; j <= jmax; j++){
            fprintf(arquivo, "%lf ", dev_t[idx]);
        }
        fprintf(arquivo, "\n");
    }
    fclose(arquivo);
    
    snprintf(filepath, sizeof(filepath), "transient/data/C%s.dat", filename2);
    arquivo = fopen(filepath, "w");
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            fprintf(arquivo, "%lf ",dev_c[idx]);
        }
        fprintf(arquivo, "\n");
    }
    fclose(arquivo);

    snprintf(filepath, sizeof(filepath), "transient/data/time%s.dat", filename2);
    arquivo = fopen(filepath, "w");
    fprintf(arquivo, "%lf\n", tempo);
    fclose(arquivo);
    
    arquivo = fopen("transient/data/grid.dat", "w");
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            fprintf(arquivo, "%lf %lf\n", dev_x[i], dev_y[j]);
        }
    }
    fclose(arquivo);

//########## OUTPUT TO USE IN PARAVIEW APP ################################
    snprintf(filepath, sizeof(filepath), "data/paraview_output%s.vtk", filename2);
    arquivo = fopen(filepath, "w");

    fprintf(arquivo, "# vtk DataFile Version 3.0\n");
    fprintf(arquivo, "Droplet Combustion\n");
    fprintf(arquivo, "ASCII\n");
    fprintf(arquivo, "DATASET STRUCTURED_GRID\n");
    fprintf(arquivo, "DIMENSIONS %d %d 1\n", jmax, imax);
    fprintf(arquivo, "POINTS %d float\n", imax * jmax);

    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            fprintf(arquivo, "%14.4f %14.4f %14.4f\n", dev_x[i], dev_y[j], 0.0);
        }
    }

    // Velocidade magnitude
    fprintf(arquivo, "POINT_DATA %d\n", imax * jmax);
    fprintf(arquivo, "SCALARS VEL_MAGNITUDE float\n");
    fprintf(arquivo, "LOOKUP_TABLE default\n");
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            double vel = sqrt(dev_u[idx]*dev_u[idx] + dev_v[idx]*dev_v[idx]);
            fprintf(arquivo, "%14.4f\n", vel);
        }
    }

    // Temperatura
    fprintf(arquivo, "SCALARS T float\n");
    fprintf(arquivo, "LOOKUP_TABLE default\n");
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            fprintf(arquivo, "%14.4f ",dev_t[idx]);
        }
        fprintf(arquivo, "\n");
    }

    // Pressão
    fprintf(arquivo, "SCALARS P float\n");
    fprintf(arquivo, "LOOKUP_TABLE default\n");
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            fprintf(arquivo, "%14.4f ", dev_p[idx]);
        }
        fprintf(arquivo, "\n");
    }


    fprintf(arquivo, "\n\n");
    fprintf(arquivo, "VECTORS Vectors float\n");
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            fprintf(arquivo, "%14.4f %14.4f %14.4f ", dev_u[idx], dev_v[idx], 0.0);
        }
        fprintf(arquivo, "\n");
    }
    fclose(arquivo);
}