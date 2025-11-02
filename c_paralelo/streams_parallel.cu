#include "comum.h"

// Streams globais para paralelização
cudaStream_t stream1, stream2, stream3, stream4;

void init_streams() {
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);
    cudaStreamCreate(&stream3);
    cudaStreamCreate(&stream4);
}

void destroy_streams() {
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);
    cudaStreamDestroy(stream4);
}

// Kernels assíncronos usando streams
void RESU_ASYNC(double *dev_um, double *dev_vm, double *dev_p, double *dev_ru, cudaStream_t stream){
    dim3 blockDim(16,16);
    dim3 gridDim((imax-3 + blockDim.x - 1)/blockDim.x, (jmax-3 + blockDim.y - 1)/blockDim.y);
    
    // Lançar kernel calc_resu no stream especificado
    calc_resu<<<gridDim, blockDim, 0, stream>>>(
        dev_areau_e, dev_areau_w, dev_areav_n, dev_areav_s, dev_epsilon1, dev_y, 
        dev_x, dev_xm, dev_ym, dev_liga_poros, iterations.b_art, imax, jmax, re, darcy_number, 
        cf, invfr2, dev_um, dev_vm, dev_p, dev_ru);
    
    // Lançar upwind em paralelo no mesmo stream
    upwind_U_pair_async<<<gridDim, blockDim, 0, stream>>>(dev_um, dev_vm, dev_p, dev_ru);
}

void RESV_ASYNC(double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_rv, cudaStream_t stream){
    dim3 blockDim(16,16);
    dim3 gridDim((imax-3 + blockDim.x - 1)/blockDim.x, (jmax-3 + blockDim.y - 1)/blockDim.y);
    
    calc_resv<<<gridDim, blockDim, 0, stream>>>(
        dev_areav_e, dev_areav_n, dev_areav_s, dev_areav_w, dev_epsilon1, dev_y, 
        dev_x, dev_xm, dev_ym, dev_liga_poros, iterations.b_art, imax, jmax, re, darcy_number, 
        cf, invfr2, dev_um, dev_vm, dev_p, dev_t, dev_rv);
        
    upwind_V_pair_async<<<gridDim, blockDim, 0, stream>>>(dev_um, dev_vm, dev_p, dev_t, dev_rv);
}

void RESC_ASYNC(double *dev_um_n, double *dev_vm_n, double *dev_c, double *dev_rc, cudaStream_t stream){
    dim3 blockDim(16,16);
    dim3 gridDim((imax-3 + blockDim.x - 1)/blockDim.x, (jmax-3 + blockDim.y - 1)/blockDim.y);
    
    calc_resc_shared<<<gridDim, blockDim, 0, stream>>>(dev_areau_e, dev_areau_w, dev_areav_n, dev_areav_s, dev_xm, dev_ym, dev_x, dev_y, 
    dev_liga_poros, dev_epsilon1, re, sc, imax, jmax, dev_um_n, dev_vm_n, dev_c, dev_rc);
}

// Execução paralela de múltiplos kernels
void SOLVE_PARALLEL(double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_c, 
                    double *dev_um_n, double *dev_vm_n, double *dev_ru, double *dev_rv, double *dev_rc) {
    
    // Executar kernels independentes em streams paralelos
    RESU_ASYNC(dev_um, dev_vm, dev_p, dev_ru, stream1);
    RESV_ASYNC(dev_um, dev_vm, dev_p, dev_t, dev_rv, stream2);
    RESC_ASYNC(dev_um_n, dev_vm_n, dev_c, dev_rc, stream3);
    
    // Sincronizar todos os streams
    cudaStreamSynchronize(stream1);
    cudaStreamSynchronize(stream2);
    cudaStreamSynchronize(stream3);
}