#include "comum.h"

// Versão SIMPLES e SEGURA da fusão upwind
void upwind_U_parallel(double *dev_um, double *dev_vm, double *dev_p, double *dev_ru){
    // Lançar kernels originais em paralelo usando streams
    cudaStream_t stream1, stream2, stream3, stream4;
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);
    cudaStreamCreate(&stream3);
    cudaStreamCreate(&stream4);
    
    // Lançar todos em paralelo (sem fusão, apenas paralelismo)
    upwind_Ui_stream(dev_um, dev_vm, dev_p, dev_ru, 2, stream1);
    upwind_Ui_stream(dev_um, dev_vm, dev_p, dev_ru, jmax-1, stream2);
    upwind_Uj_stream(dev_um, dev_vm, dev_p, dev_ru, 2, stream3);
    upwind_Uj_stream(dev_um, dev_vm, dev_p, dev_ru, imax, stream4);
    
    // Sincronizar
    cudaStreamSynchronize(stream1);
    cudaStreamSynchronize(stream2);
    cudaStreamSynchronize(stream3);
    cudaStreamSynchronize(stream4);
    
    // Cleanup
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);
    cudaStreamDestroy(stream4);
}

// Versões com stream dos kernels originais
void upwind_Ui_stream(double *dev_um, double *dev_vm, double *dev_p, double *dev_ru, int j, cudaStream_t stream){
    int threads = 256;
    dim3 blocks = grid_1d(imax-2, threads);
    
    calc_upwind_Ui<<<blocks, threads, 0, stream>>>(dev_vm, dev_um, dev_areau_n, dev_areau_s,
         dev_areau_e, dev_areau_w, dev_epsilon1, dev_ym, dev_x, dev_y,
         dev_xm, dev_liga_poros, re, dev_p, dev_ru, imax, jmax, j, 
         iterations.b_art, darcy_number, cf, g);
}

void upwind_Uj_stream(double *dev_um, double *dev_vm, double *dev_p, double *dev_ru, int i, cudaStream_t stream){
    int threads = 256;
    dim3 blocks = grid_1d(jmax-2, threads);
    
    calc_upwind_Uj<<<blocks, threads, 0, stream>>>(dev_vm, dev_um, dev_areav_n, dev_areav_s,
         dev_areav_e, dev_areav_w, dev_epsilon1, dev_xm, dev_x, dev_y,
         dev_ym, dev_liga_poros, re, dev_p, dev_ru, imax, jmax, i,
         iterations.b_art, darcy_number, cf, g);
}