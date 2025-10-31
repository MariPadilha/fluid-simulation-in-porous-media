#include "comum.h"

__global__ void x_parte1(double *dev_x, double dx_c, double lhori, int imax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    if(i <= imax){
        dev_x[i] = (((double)(i) - 1.0) * dx_c) - lhori/(double)(2.0);
    }
}

__global__ void x_parte2(double *dev_x, double px_grid, double q_grid, int imax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + (int)(imax/2+1);
    if(i <= imax){
        double etax = (dev_x[i] - dev_x[imax]) / (dev_x[(int)(imax/2+1)] - dev_x[imax]);
        double sx = px_grid * etax + (1.0 - px_grid)
           * (1.0 - ((tanh(q_grid * (1.0 - etax))) / tanh(q_grid)));
        dev_x[i] = dev_x[imax] - sx * (dev_x[imax] - dev_x[(int)(imax/2+1)]);
    }
}

__global__ void x_parte3(double *dev_x, double px_grid, double q_grid, int imax){
    int idx = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int i_start = (int)(imax/2+1);

    if(idx <= i_start){
        int i = i_start - idx;
        double etax = (dev_x[i] + dev_x[imax]) / (dev_x[i_start] + dev_x[imax]);
        double sx = px_grid * etax + (1.0 - px_grid)
            *  (1.0 - (tanh(q_grid * (1.0 - etax)) / tanh(q_grid)));
        dev_x[i] = -dev_x[imax] - sx * (-dev_x[imax] - dev_x[i_start]);
    }
}

__global__ void y_parte1(double *dev_y, double dx_c, double y_down, int jmax){
    int j = blockIdx.x * blockDim.x + threadIdx.x + 1;
    if(j <= jmax)
        dev_y[j] = (((double)(j) - 1.0) * dx_c) - y_down;
}

__global__ void y_parte2(double *dev_y, double y_down, double dx_c, double py_grid, double q_grid, int jmax){
    int j = blockIdx.x * blockDim.x + threadIdx.x + ((int)(y_down/dx_c)+1);
    if(j > jmax)return;

    double etay = (dev_y[j] - dev_y[jmax]) / (dev_y[((int)(y_down/dx_c)+1)] - dev_y[jmax]);
    double sy = py_grid * etay + (1.0 - py_grid)
           * (1.0 - (tanh((q_grid * (1.0 - etay))) / tanh(q_grid)));
    dev_y[j] = dev_y[jmax] - sy * (dev_y[jmax] - dev_y[(int)(y_down/dx_c)+1]);
}

__global__ void y_parte3(double *dev_y, double y_down, double dx_c, double py_grid, double q_grid, int jmax){
    int idx = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j_start = (int)(y_down/dx_c) + 1;

    if(idx > j_start) return;

    int j = j_start - idx;

    double etay = (dev_y[j] + dev_y[jmax]) / (dev_y[j_start] + dev_y[jmax]);
    double sy = py_grid * etay + (1.0 - py_grid) * (1.0 - (tanh(q_grid * (1.0 - etay)) / tanh(q_grid)));
    dev_y[j] = -dev_y[jmax] - sy * (-dev_y[jmax] - dev_y[j_start]);
}

void calcula_x_y(){
    int threads = 256;
    int blocks = grid_1d(imax, threads);
    x_parte1<<<blocks, threads>>>(dev_x, dx_c, lhori, imax);
    
    blocks = grid_1d(jmax, threads);
    y_parte1<<<blocks, threads>>>(dev_y, dx_c, y_down, jmax);
    
    blocks = grid_1d((imax+1)-(int)(imax/2+1), threads);
    x_parte2<<<blocks, threads>>>(dev_x, px_grid, q_grid, imax);

    blocks = grid_1d(imax/2, threads);
    x_parte3<<<blocks, threads>>>(dev_x, px_grid, q_grid, imax);
    
    blocks = grid_1d((jmax+1)-(int)(y_down/dx_c+1), threads);
    y_parte2<<<blocks, threads>>>(dev_y, y_down, dx_c, py_grid, q_grid, jmax);

    blocks = grid_1d(((int)(y_down/dx_c)+1), threads);
    y_parte3<<<blocks, threads>>>(dev_y, y_down, dx_c, py_grid, q_grid, jmax);
}