#include "comum.h"

//--- ResC com Shared Memory ---
__global__ void calc_resc_shared(double *dev_areau_e, double *dev_areau_w, double *dev_areav_n, double *dev_areav_s, double *dev_xm, 
    double *dev_ym, double *dev_x, double *dev_y, double *dev_liga_poros, double *dev_epsilon1, double re, double sc, int imax,  
    int jmax, double *dev_um_n, double *dev_vm_n, double *dev_c, double *dev_rc){
        
    int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
    
    // Shared memory para cache dos dados mais acessados
    __shared__ double shared_c[18][18];  // 16+2 para halo
    __shared__ double shared_um[19][18]; // 16+3 para halo U
    __shared__ double shared_vm[18][19]; // 16+3 para halo V
    
    int tx = threadIdx.x + 1;
    int ty = threadIdx.y + 1;
    
    // Carregar dados para shared memory
    if(i <= imax-1 && j <= jmax-1) {
        int idx = i*(jmax+1)+j;
        shared_c[tx][ty] = dev_c[idx];
        shared_um[tx][ty] = dev_um_n[idx];
        shared_vm[tx][ty] = dev_vm_n[i*(jmax+2)+j];
        
        // Carregar halo (bordas)
        if(threadIdx.x == 0 && i > 2) {
            shared_c[0][ty] = dev_c[(i-1)*(jmax+1)+j];
            shared_um[0][ty] = dev_um_n[(i-1)*(jmax+1)+j];
        }
        if(threadIdx.x == blockDim.x-1 && i < imax-1) {
            shared_c[tx+1][ty] = dev_c[(i+1)*(jmax+1)+j];
            shared_um[tx+2][ty] = dev_um_n[(i+1)*(jmax+1)+j];
        }
        if(threadIdx.y == 0 && j > 2) {
            shared_c[tx][0] = dev_c[i*(jmax+1)+(j-1)];
            shared_vm[tx][0] = dev_vm_n[i*(jmax+2)+(j-1)];
        }
        if(threadIdx.y == blockDim.y-1 && j < jmax-1) {
            shared_c[tx][ty+1] = dev_c[i*(jmax+1)+(j+1)];
            shared_vm[tx][ty+2] = dev_vm_n[i*(jmax+2)+(j+1)];
        }
    }
    
    __syncthreads();
    
    if(i <= imax-1 && j <= jmax-1){
        int idx = i*(jmax+1)+j;
        
        // Variáveis locais otimizadas
        double xm_ip1 = dev_xm[i+1];
        double xm_i = dev_xm[i];
        double ym_jp1 = dev_ym[j+1];
        double ym_j = dev_ym[j];
        double x_i = dev_x[i];
        double x_ip1 = dev_x[i+1];
        double x_im1 = dev_x[i-1];
        double y_j = dev_y[j];
        double y_jp1 = dev_y[j+1];
        double y_jm1 = dev_y[j-1];
        
        // Usar dados da shared memory
        double c_center = shared_c[tx][ty];
        double c_east = shared_c[tx+1][ty];
        double c_west = shared_c[tx-1][ty];
        double c_north = shared_c[tx][ty+1];
        double c_south = shared_c[tx][ty-1];
        
        double um_center = shared_um[tx][ty];
        double um_east = shared_um[tx+1][ty];
        double vm_north = shared_vm[tx][ty+1];
        double vm_center = shared_vm[tx][ty];
        
        // Cálculos otimizados
        double dcudx = 0.5 * (c_east + c_center) * um_east * dev_areau_e[j]
                     - 0.5 * (c_west + c_center) * um_center * dev_areau_w[j];
                     
        double dcvdy = 0.5 * (c_north + c_center) * vm_north * dev_areav_n[i]
                     - 0.5 * (c_south + c_center) * vm_center * dev_areav_s[i];

        double aux = 1.0/(re*sc);
        double de = (ym_jp1 - ym_j) * aux / (x_ip1 - x_i);
        double dw = (ym_jp1 - ym_j) * aux / (x_i - x_im1);
        double dn = (xm_ip1 - xm_i) * aux / (y_jp1 - y_j);
        double ds = (xm_ip1 - xm_i) * aux / (y_j - y_jm1);
        double dp = de + dw + dn + ds;

        double vol_inv = 1.0 / ((xm_ip1 - xm_i) * (ym_jp1 - ym_j));
        double poros_factor = 1.0 / (dev_liga_poros[idx] * (dev_epsilon1[idx] - 1.0) + 1.0);

        dev_rc[idx] = vol_inv * (-dp * c_center + de * c_east + dw * c_west + 
                                 dn * c_north + ds * c_south - (dcudx + dcvdy) * poros_factor);
    }
}

void RESC_SHARED(double *dev_um_n, double *dev_vm_n, double *dev_c, double *dev_rc){
    dim3 blockDim(16,16);
    dim3 gridDim((imax-3 + blockDim.x - 1)/blockDim.x, (jmax-3 + blockDim.y - 1)/blockDim.y);
    
    calc_resc_shared<<<gridDim, blockDim>>>(dev_areau_e, dev_areau_w, dev_areav_n, dev_areav_s, dev_xm, dev_ym, dev_x, dev_y, 
    dev_liga_poros, dev_epsilon1, re, sc, imax, jmax, dev_um_n, dev_vm_n, dev_c, dev_rc);
}