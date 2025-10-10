#include "comum.h"
#define idx i*(jmax+1)+j

__global__ void atualiza_ic_um(double *dev_um, int imax, int jmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    if(i <= imax+1 && j <= jmax){
        dev_um[idx] = 0.0;
    }
}

__global__ void atualiza_ic_vm(double *dev_vm, double v_i, int imax, int jmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    if(i <= imax && j <= jmax+1){
        dev_vm[i*(jmax+2)+j] = (1.0e-3)*(v_i);
    }
}

__global__ void atualiza_ic_ptc(double *dev_p, double *dev_t, double *dev_c, int *dev_flag, double tinf, double temp_cylinder, double concentracao_inicial, int c_f, int imax, int jmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if(i > imax || j > jmax)return;

    dev_p[idx] = 1.0;
    dev_t[idx] = tinf;
    dev_c[idx] = 0.0;
    
    if(dev_flag[idx] != c_f){
        dev_t[idx] = temp_cylinder;
        dev_c[idx] = concentracao_inicial;
    }
}

__global__ void restart_dom_um(double *dev_um, double *dev_umr, int rjmax, int imax, int jmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    if(i <= imax+1 && j <= jmax){
        dev_um[i*(jmax+1)+j] = (j <= rjmax) ? dev_umr[i*(rjmax+1)+j] : dev_um[i*(jmax+1)+(j-1)];
    }
}

__global__ void restart_dom_vm(double *dev_vm, double *dev_vmr, int imax, int jmax, int rjmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    if(i <= imax && j <= jmax+1){
        dev_vm[i*(jmax+2)+j] = (j <= rjmax) ? dev_vmr[i*(rjmax+2)+j] : dev_vm[i*(jmax+2)+(j-1)];
    }
}

__global__ void restart_dom_ptzh(double *dev_p, double *dev_t, double *dev_z, double *dev_h, double *dev_pres, double *dev_tr, double *dev_zr, double *dev_hr, int imax, int jmax, int rjmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if(i > imax || j > jmax)return;

    dev_p[idx] = (j <= rjmax) ? dev_pres[i*(rjmax+1)+j] : dev_p[i*(jmax+1)+(j-1)];
    dev_t[idx] = (j <= rjmax) ? dev_tr[i*(rjmax+1)+j] : dev_t[i*(jmax+1)+(j-1)];
    dev_z[idx] = (j <= rjmax) ? dev_zr[i*(rjmax+1)+j] : dev_z[i*(jmax+1)+(j-1)];
    dev_h[idx] = (j <= rjmax) ? dev_hr[i*(rjmax+1)+j] : dev_h[i*(jmax+1)+(j-1)];       
}

__global__ void atualiza_hr(double *dev_hr, double *dev_h_res, int rimax, int rjmax, double s, double lf, double tinf, double q){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if(i > rimax || j > rjmax)return;

    dev_hr[i*(rjmax+1)+j] = dev_h_res[i*(rjmax+1)+j] + (((s + 1.0) * lf * tinf / q + 1.0) - dev_h_res[i*(rjmax+1)+j]);
}

void init(){
    FILE *arquivo;

    //--- iterations values ---
    arquivo = fopen("input/iterations.dat", "r+");
    (void)fscanf(arquivo, "%d %d %d %d %lf %lf %lf %lf %lf %lf %d",
           &iterations.nc, &iterations.n_tr, &iterations.n_out, &iterations.n_vort,
           &iterations.beta, &iterations.b_art, &iterations.dtau_f, &iterations.final_time,
           &iterations.eps, &iterations.eps_mass, &iterations.start_mode);
    fclose(arquivo);
    printf("itc_max: %d, nc: %d, n_tr: %d, n_out: %d, n_vort: %d\nbeta: %lf, b_art: %lf, dtau_f: %lf, final_time: %lf,\neps: %lf, eps_mass: %lf, start_mode: %d\n",
        iterations.itc_max, iterations.nc, iterations.n_tr, iterations.n_out, iterations.n_vort,
        iterations.beta, iterations.b_art, iterations.dtau_f, iterations.final_time,
        iterations.eps, iterations.eps_mass, iterations.start_mode);

    //--- reference values ---
    arquivo = fopen("input/reference.dat", "r+");
    (void)fscanf(arquivo, "%lf %lf %lf %lf %lf", &ref.tnu, &ref.yf_b, &ref.yo_oo, &ref.ts, &ref.tn_too);
    fclose(arquivo);
    printf("tnu: %lf, yf_b: %lf, yo_oo: %lf, ts: %lf, tn_too: %lf\n", ref.tnu, ref.yf_b, ref.yo_oo, ref.ts, ref.tn_too);

}

void IC(double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_c, double *dev_pn){
    dim3 blockDim(16,16);
    dim3 gridDimUm((imax+1 + blockDim.x - 1)/blockDim.x, (jmax + blockDim.y - 1)/blockDim.y);
    dim3 gridDimVm((imax + blockDim.x - 1)/blockDim.x, (jmax+1 + blockDim.y - 1)/blockDim.y);
    dim3 gridDimPtc((imax + blockDim.x - 1)/blockDim.x, (jmax + blockDim.y - 1)/blockDim.y);

    cudaStream_t s1, s2, s3;
    cudaStreamCreate(&s1);
    cudaStreamCreate(&s2);
    cudaStreamCreate(&s3);

    atualiza_ic_um<<<gridDimUm, blockDim, 0, s1>>>(dev_um, imax, jmax);
    atualiza_ic_vm<<<gridDimVm, blockDim, 0, s2>>>(dev_vm, v_i, imax, jmax);
    atualiza_ic_ptc<<<gridDimPtc, blockDim, 0, s3>>>(dev_p, dev_t, dev_c, dev_flag, tinf, temp_cylinder, concentracao_inicial, c_f, imax, jmax);

    cudaStreamSynchronize(s1);
    cudaStreamSynchronize(s2);
    cudaStreamSynchronize(s3);
    
    cudaStreamDestroy(s1);
    cudaStreamDestroy(s2);
    cudaStreamDestroy(s3);

    remove("data/flametip.dat");
    remove("data/error.dat");
}

void restart(double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_c){
    int i, j;
    FILE *arquivo;

    printf("RESTARTING PROGRAM\n");
    
    arquivo = fopen("data/restart/restartU.dat", "r");
    if (!arquivo){
        perror("Erro ao abrir data/restart/restartU.dat");
        exit(1);
    }
    for(i = 1; i <= imax+1; i++){
        for(j = 1; j <= jmax; j++){
            (void)fscanf(arquivo, "%lf", &(dev_um[idx]));
        }
    }
    fclose(arquivo);

    arquivo = fopen("data/restart/restartV.dat", "r");
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax+1; j++){
            (void)fscanf(arquivo, "%lf", &(dev_vm[i*(jmax+2)+j]));
        }
    }
    fclose(arquivo);

    arquivo = fopen("data/restart/restartPTC.dat", "r");
    if (!arquivo){
        perror("Erro ao abrir data/restart/restartPTC.dat");
        exit(1);
    }
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            (void)fscanf(arquivo, "%lf %lf %lf", &(dev_p[idx]), &(dev_t[idx]), &(dev_c[idx]));
        }
    }
    fclose(arquivo);
}

void restart_dom(double *dev_um, double *dev_vm, double *dev_p, double *dev_t, double *dev_z, double *dev_h){
    int i, j;
    int rimax = 41;  //em x
    int rjmax = 321; //em y
    double *dev_umr, *dev_vmr, *dev_pres, *dev_zr, *dev_tr, *dev_hr, *dev_h_res;
    FILE *arquivo;
    dim3 blocks, threads(16,16);

    cudaMalloc((void**)&dev_umr, sizeof(double)*(rimax+2)*(rjmax+1));
    cudaMalloc((void**)&dev_vmr, sizeof(double)*(rimax+1)*(rjmax+2));
    cudaMalloc((void**)&dev_pres, sizeof(double*)*(rimax+1)*(rjmax+1));
    cudaMalloc((void**)&dev_zr, sizeof(double)*(rimax+1)*(rjmax+1));
    cudaMalloc((void**)&dev_tr, sizeof(double)*(rimax+1)*(rjmax+1));
    cudaMalloc((void**)&dev_hr, sizeof(double)*(rimax+1)*(rjmax+1));
    cudaMalloc((void**)&dev_h_res, sizeof(double)*(rimax+1)*(rjmax+1));

    printf("RESTARTING PROGRAM\n");

    arquivo = fopen("data/restart/restartU.dat", "r");
    for(i = 1; i <= rimax+1; i++){
        for(j = 1; j <= rjmax; j++){
            (void)fscanf(arquivo, "%lf", &(dev_umr[i*(rjmax+1)+j]));
        }
    }
    fclose(arquivo);

    arquivo = fopen("data/restart/restartV.dat", "r");
    for(i = 1; i <= rimax; i++){
        for(j = 1; j <= rjmax+1; j++){
            (void)fscanf(arquivo, "%lf", &(dev_vmr[i*(rjmax+2)+j]));
        }
    }
    fclose(arquivo);

    arquivo = fopen("data/restart/restartPTZH.dat", "r");
    for(i = 1; i <= rimax; i++){
        for(j = 1; j <= rjmax; j++){
            (void)fscanf(arquivo, "%lf %lf %lf %lf", &(dev_pres[i*(rjmax+1)+j]), &(dev_tr[i*(rjmax+1)+j]), &(dev_zr[i*(rjmax+1)+j]), &(dev_h_res[i*(rjmax+1)+j]));
        }
    }
    fclose(arquivo);

    blocks = grid_2d(rimax-1, rjmax-1);
    atualiza_hr<<<blocks, threads>>>(dev_hr, dev_h_res, rimax, rjmax, s, lf, tinf, q);

    cudaStream_t s1, s2, s3;
    cudaStreamCreate(&s1);
    cudaStreamCreate(&s2);
    cudaStreamCreate(&s3);

    blocks = grid_2d(imax, jmax-1);
    restart_dom_um<<<blocks, threads, 0, s1>>>(dev_um, dev_umr, rjmax, imax, jmax);
    blocks = grid_2d(imax-1, jmax);
    restart_dom_vm<<<blocks, threads, 0, s2>>>(dev_vm, dev_vmr, imax, jmax, rjmax);
    blocks = grid_2d(imax-1, jmax-1);
    restart_dom_ptzh<<<blocks, threads, 0, s3>>>(dev_p, dev_t, dev_z, dev_h, dev_pres, dev_tr, dev_zr, dev_hr, imax, jmax, rjmax);

    cudaStreamSynchronize(s1);
    cudaStreamSynchronize(s2);
    cudaStreamSynchronize(s3);
    
    cudaStreamDestroy(s1);
    cudaStreamDestroy(s2);
    cudaStreamDestroy(s3);

//desalocando
    cudaFree(dev_umr);
    cudaFree(dev_vmr);
    cudaFree(dev_pres);
    cudaFree(dev_zr);
    cudaFree(dev_tr);
    cudaFree(dev_hr);
    cudaFree(dev_h_res);
}