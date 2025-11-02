#include "comum.h"

// Template para especialização de kernels por tamanho
template<int BLOCK_SIZE>
__global__ void calc_resu_unrolled(double *dev_areau_e, double *dev_areau_w, double *dev_areav_n, double *dev_areav_s, 
    double *dev_epsilon1, double *dev_y, double *dev_x, double *dev_xm, double *dev_ym, 
    double *dev_liga_poros, double b_art, int imax, int jmax, double re, double darcy_number, 
    double cf, double invfr2, double *dev_um, double *dev_vm, double *dev_p, double *dev_ru) {
    
    int i = blockIdx.x * BLOCK_SIZE + threadIdx.x + 2;
    int j = blockIdx.y * BLOCK_SIZE + threadIdx.y + 2;
    
    // Desenrolar loops internos para eliminar overhead
    if(i <= imax-1 && j <= jmax-1) {
        int idx = i*(jmax+1)+j;
        
        // Variáveis locais pré-carregadas (elimina redundância)
        const double re_inv = 1.0/re;
        const double cf_re = cf * re_inv;
        
        // Pre-calcular todas as coordenadas (elimina re-cálculos)
        const double x_coords[3] = {dev_x[i-1], dev_x[i], dev_x[i+1]};
        const double y_coords[3] = {dev_y[j-1], dev_y[j], dev_y[j+1]};
        const double xm_coords[2] = {dev_xm[i], dev_xm[i+1]};
        const double ym_coords[2] = {dev_ym[j], dev_ym[j+1]};
        
        // Desenrolar cálculos de diferenças
        const double dx_e = x_coords[2] - x_coords[1];
        const double dx_w = x_coords[1] - x_coords[0]; 
        const double dy_n = y_coords[2] - y_coords[1];
        const double dy_s = y_coords[1] - y_coords[0];
        const double dx_vol = xm_coords[1] - xm_coords[0];
        const double dy_vol = ym_coords[1] - ym_coords[0];
        
        // Cálculos de coeficientes com constantes pré-calculadas
        const double vol_factor = dy_vol * re_inv;
        const double de = vol_factor / dx_e;
        const double dw = vol_factor / dx_w;
        const double vol_factor_x = dx_vol * re_inv;
        const double dn = vol_factor_x / dy_n;
        const double ds = vol_factor_x / dy_s;
        const double dp = de + dw + dn + ds;
        
        // Desenrolar acessos de velocidade (elimina cálculos de índice redundantes)
        const double u_center = dev_um[idx];
        const double u_east = dev_um[(i+1)*(jmax+1)+j];
        const double u_west = dev_um[(i-1)*(jmax+1)+j]; 
        const double u_north = dev_um[i*(jmax+1)+(j+1)];
        const double u_south = dev_um[i*(jmax+1)+(j-1)];
        
        // Cálculo final otimizado
        const double vol_inv = 1.0 / (dx_vol * dy_vol);
        const double epsilon_factor = dev_epsilon1[idx];
        const double poros_factor = dev_liga_poros[idx] * (epsilon_factor - 1.0) + 1.0;
        
        dev_ru[idx] = vol_inv * (
            -dp * u_center + de * u_east + dw * u_west + 
            dn * u_north + ds * u_south
        ) / poros_factor;
    }
}

// Lançamento otimizado com template
void RESU_UNROLLED(double *dev_um, double *dev_vm, double *dev_p, double *dev_ru) {
    constexpr int BLOCK_SIZE = 16;
    dim3 blockDim(BLOCK_SIZE, BLOCK_SIZE);
    dim3 gridDim((imax-3 + BLOCK_SIZE - 1)/BLOCK_SIZE, (jmax-3 + BLOCK_SIZE - 1)/BLOCK_SIZE);
    
    calc_resu_unrolled<BLOCK_SIZE><<<gridDim, blockDim>>>(
        dev_areau_e, dev_areau_w, dev_areav_n, dev_areav_s, dev_epsilon1, dev_y, 
        dev_x, dev_xm, dev_ym, dev_liga_poros, iterations.b_art, imax, jmax, re, darcy_number, 
        cf, invfr2, dev_um, dev_vm, dev_p, dev_ru);
}