//Artificial Compressibility Methods
//Solve Momentum Equation with QUICK Scheme
//BASED ON:
//Versteeg, H. K., and W. Malalasekera. 
//"An introduction to computational Fluid Dynamics, The finite volume control, ed." (1995).
#include "comum.h"
#include <math.h>
#define N_IMAX 51*4
#define N_ITC 800000
#define idx i*(jmax+1)+j

__global__ void atualizar_matrizes_linearizadas(double *origem, double *destino, int tamanhoLinha, int tamanhoColuna, int inicio, int coluna){
    int i = blockIdx.x * blockDim.x + threadIdx.x + inicio;
    int j = blockIdx.y * blockDim.y + threadIdx.y + inicio;
    
    if(i <= tamanhoLinha && j <= tamanhoColuna){
        destino[i*(coluna)+j] = origem[i*(coluna)+j];
    }
}

__global__ void atualiza_tc(double *dev_t, double *dev_c, int *dev_flag, double temp_cylinder, double concentracao_inicial, int c_f, int imax, int jmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    if(i > imax || j > jmax) return;
    
    if(dev_flag[idx] != c_f){
        dev_t[idx] = temp_cylinder;               
        dev_c[idx] = concentracao_inicial;
    }
}

int main(int argc, char *argv[]){
    int itc, tr;
    int n_imax, n_itc;
    dim3 blockDim(16, 16);  
    
    if(argc < 3){
        n_imax = N_IMAX;
        n_itc = N_ITC;
    }else{
        n_imax = atoi(argv[1]);
        n_itc = atoi(argv[2]);
    }
    
    calcular(n_imax, n_itc);         // define imax, jmax, dx_c, etc.
    alocar_globais(); 
    
    dim3 gridDim((imax + blockDim.x - 1)/blockDim.x, (jmax + blockDim.y - 1)/blockDim.y);
    dim3 gridDimUm((imax+1 + blockDim.x - 1)/blockDim.x, (jmax + blockDim.y - 1)/blockDim.y);
    dim3 gridDimVm((imax + blockDim.x - 1)/blockDim.x, (jmax+1 + blockDim.y - 1)/blockDim.y);
    double *dev_rz;
    double *dev_vi, *dev_rv, *dev_res_v;
    double *dev_ui, *dev_ru, *dev_res_u;
    double *dev_rc;
    double *dev_um, *dev_vm;
    double *dev_um_n, *dev_vm_n;
    double *dev_um_tau, *dev_vm_tau;
    double *dev_um_n_tau, *dev_vm_n_tau;
    cudaMalloc((void**)&dev_rc, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_ru, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_ui, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_res_u, sizeof(double)*(imax+2)*(jmax+1));
    cudaMallocManaged((void**)&dev_rv, sizeof(double)*(imax+1)*(jmax+2));
    cudaMallocManaged((void**)&dev_vi, sizeof(double)*(imax+1)*(jmax+2));
    cudaMallocManaged((void**)&dev_res_v, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_rz, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_um, sizeof(double)*(imax+2)*(jmax+1));
    cudaMallocManaged((void**)&dev_vm, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_um_n, sizeof(double)*(imax+2)*(jmax+1));
    cudaMallocManaged((void**)&dev_um_tau, sizeof(double)*(imax+2)*(jmax+1));
    cudaMallocManaged((void**)&dev_um_n_tau, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_vm_n, sizeof(double)*(imax+1)*(jmax+2));
    cudaMallocManaged((void**)&dev_vm_tau, sizeof(double)*(imax+1)*(jmax+2));
    cudaMallocManaged((void**)&dev_vm_n_tau, sizeof(double)*(imax+1)*(jmax+2));
    
    double *dev_u, *dev_v, *dev_p, *dev_pn, *dev_h, *dev_t, *dev_z;
    double *dev_c, *dev_t_n_tau, *dev_t_tau, *dev_c_n_tau, *dev_c_tau; 
    cudaMallocManaged((void**)&dev_u, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_v, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_p, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_t, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_c, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_pn, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_h, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_z, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_t_n_tau, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_t_tau, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_c_n_tau, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_c_tau, sizeof(double)*(imax+1)*(jmax+1));

    double residual_p, residual_u, residual_v, error;

    //=================================================================
    // CUDA Events — Instrumentação de tempo (NÃO altera lógica numérica)
    //=================================================================
    cudaEvent_t ev_ts_start, ev_ts_stop;          // timestep total
    cudaEvent_t ev_solveU_start, ev_solveU_stop;  // solve_U
    cudaEvent_t ev_solveV_start, ev_solveV_stop;  // solve_V
    cudaEvent_t ev_solveP_start, ev_solveP_stop;  // solve_P
    cudaEvent_t ev_solveZ_start, ev_solveZ_stop;  // solve_Z
    cudaEvent_t ev_solveC_start, ev_solveC_stop;  // solve_C
    cudaEvent_t ev_h2d_start, ev_h2d_stop;        // H2D (prefetch)
    cudaEvent_t ev_d2h_start, ev_d2h_stop;        // D2H (sync back)

    cudaEventCreate(&ev_ts_start);    cudaEventCreate(&ev_ts_stop);
    cudaEventCreate(&ev_solveU_start); cudaEventCreate(&ev_solveU_stop);
    cudaEventCreate(&ev_solveV_start); cudaEventCreate(&ev_solveV_stop);
    cudaEventCreate(&ev_solveP_start); cudaEventCreate(&ev_solveP_stop);
    cudaEventCreate(&ev_solveZ_start); cudaEventCreate(&ev_solveZ_stop);
    cudaEventCreate(&ev_solveC_start); cudaEventCreate(&ev_solveC_stop);
    cudaEventCreate(&ev_h2d_start);   cudaEventCreate(&ev_h2d_stop);
    cudaEventCreate(&ev_d2h_start);   cudaEventCreate(&ev_d2h_stop);

    float time_solveU, time_solveV, time_solveP, time_solveZ, time_solveC;
    float time_h2d, time_d2h, time_timestep;
    // Acumuladores por timestep (soma das sub-iterações pseudo-tempo)
    float acc_solveU, acc_solveV, acc_solveP, acc_solveZ, acc_solveC;
    float acc_h2d, acc_d2h;

    //duration = omp_get_wtime()
    /*  character(len=128) :: pwd
        REAL ETIME, clockTIME, TARRAY(2)
        clockTIME = ETIME(TARRAY)
        CALL idate(hoje)   ! hoje(1)=day, (2)=month, (3)=year
        CALL itime(agora)     ! agora(1)=hour, (2)=minute, (3)=second
        CALL get_environment_variable('PWD', pwd)

    --- Input data - Fill NAMELIST iterations and ref
    --- Read iterations.dat and reference.dat
    --- Compute Too, Tsub, Tinf 
    */
    init(); 

    too = ref.ts * ref.tn_too;
    tsup = ref.ts / ref.ts;
    tinf = 1.0;  //Temperatura ambiente

    /*--- Log de informações ---
        5 FORMAT ( 'Date ', i2.2, '/', i2.2, '/', i4.4, &
                '; time ',i2.2, ':', i2.2, ':', i2.2 )
        WRITE(*,5)  hoje(2), hoje(1), hoje(3), agora
        WRITE(*,*) 'Current working directory: ',trim(pwd)
        WRITE(*,*) '----------------------'
        WRITE(*,*) 'Tsup =', Tsup
        WRITE(*,*) 'Tinf =', Tinf
        WRITE(*,*) '----------------------'
        WRITE(*,*) 'Mesh =',imax,'x', jmax
        WRITE(*,*) 'imax * jmax =', int(imax*jmax)
        WRITE(*,*) 'eps =',eps
        WRITE(*,*) '----------------------'
    ver se tudo isso é necessário?
    WRITE(*,*) '----------------------'
    WRITE(*,*) 'v_c ='  ,v_c        , '[m/s]'
    WRITE(*,*) 'L_c ='  ,L_c        , '[m]'
    WRITE(*,*) 't_c ='  ,L_c/v_c    , '[s]'
    WRITE(*,*) 'v_i ='  ,v_i
    WRITE(*,*) 'V_idim ='  ,v_i*v_c , '[m/s]'
    WRITE(*,*) '----------------------'
    WRITE(*,*) 'g = '  ,g , '[m/s^2]'
    WRITE(*,*) 'S = '  ,S
    WRITE(*,*) 'q = '  ,q
    WRITE(*,*) 'Pr ='  ,Pr
        WRITE(*,*) 'Re =', Re
        WRITE(*,*) 'Pe =', Pe
        WRITE(*,*) 'Fr =', Fr
        WRITE(*,*) 'InvFr^2 =', InvFr2
        WRITE(*,*) '----------------------'
    */

    //--- Initializations ---
    itc = 1; //Initial iteration
    error = 100.0;
    tr = 1;
    tempo = 0.0;
    residual_p = 0.0;

    //--- Create mesh ---
    mesh();

    atualiza_tc<<<gridDim, blockDim>>>(dev_t, dev_c, dev_flag, temp_cylinder, concentracao_inicial, c_f, imax, jmax);

    //--- Set up initial flow field ---
    if(iterations.start_mode == 0){    
        IC(dev_um, dev_vm, dev_p, dev_t, dev_c, dev_pn);
    }else if(iterations.start_mode == 1){
        restart(dev_um, dev_vm, dev_p, dev_t, dev_c);
    }else if(iterations.start_mode == 2){        
        restart_dom(dev_um, dev_vm, dev_p, dev_t, dev_z, dev_h);
    }

    //--- Pseudo time step ---
    dtau = 5.e-2;
    dt = 0.5e-2;

    atualizar_matrizes_linearizadas<<<gridDimUm, blockDim>>>(dev_um, dev_um_tau, imax+1, jmax, 1, jmax+1);
    atualizar_matrizes_linearizadas<<<gridDimVm, blockDim>>>(dev_vm, dev_vm_tau, imax, jmax+1, 1, jmax+2);

    /*duration = omp_get_wtime() - duration;
    printf("init %lf/n", duration);
    duration = omp_get_wtime();
    */
    //--- Physical time step ---
    while(tempo < iterations.final_time){
        tempo = tempo + dt;

        // --- CUDA Events: início do timestep ---
        cudaEventRecord(ev_ts_start);
        acc_solveU = 0.0f; acc_solveV = 0.0f; acc_solveP = 0.0f;
        acc_solveZ = 0.0f; acc_solveC = 0.0f;
        acc_h2d = 0.0f; acc_d2h = 0.0f;
        
        //--- Pseudo-time calculation starts ---
        while(itc < iterations.itc_max){
            //--- Solve Momentum Equation with QUICK Scheme ---
            cudaEventRecord(ev_solveU_start);
            residual_u = solve_U(dev_um, dev_vm, dev_um_n, dev_um_tau, dev_vm_tau, dev_um_n_tau, dev_pn, dev_ui, dev_ru, dev_res_u);
            cudaEventRecord(ev_solveU_stop);
            cudaEventSynchronize(ev_solveU_stop);
            cudaEventElapsedTime(&time_solveU, ev_solveU_start, ev_solveU_stop);
            acc_solveU += time_solveU;

            cudaEventRecord(ev_solveV_start);
            residual_v = solve_V(dev_um, dev_vm, dev_vm_n, dev_um_tau, dev_vm_tau, dev_vm_n_tau, dev_pn, dev_t, dev_vi, dev_rv, dev_res_v);
            cudaEventRecord(ev_solveV_stop);
            cudaEventSynchronize(ev_solveV_stop);
            cudaEventElapsedTime(&time_solveV, ev_solveV_start, ev_solveV_stop);
            acc_solveV += time_solveV;
            
            //--- Solve Continuity Equation ---
            cudaEventRecord(ev_solveP_start);
            residual_p = solve_P(dev_p, dev_um_n_tau, dev_vm_n_tau, dev_pn);
            cudaEventRecord(ev_solveP_stop);
            cudaEventSynchronize(ev_solveP_stop);
            cudaEventElapsedTime(&time_solveP, ev_solveP_start, ev_solveP_stop);
            acc_solveP += time_solveP;
            
            //--- Solve Energy Equation ---
            cudaEventRecord(ev_solveZ_start);
            solve_Z(dev_um_n_tau, dev_vm_n_tau, dev_t, dev_t_n_tau, dev_t_tau, dev_rz);
            cudaEventRecord(ev_solveZ_stop);
            cudaEventSynchronize(ev_solveZ_stop);
            cudaEventElapsedTime(&time_solveZ, ev_solveZ_start, ev_solveZ_stop);
            acc_solveZ += time_solveZ;

            cudaEventRecord(ev_solveC_start);
            solve_C(dev_um_n_tau, dev_vm_n_tau, dev_c, dev_c_n_tau, dev_c_tau, dev_rc);
            cudaEventRecord(ev_solveC_stop);
            cudaEventSynchronize(ev_solveC_stop);
            cudaEventElapsedTime(&time_solveC, ev_solveC_start, ev_solveC_stop);
            acc_solveC += time_solveC;
            /*--- check convergence ---
            CALL convergence(itc, error, residual_p, residual_u, residual_v)
            itc = itc+1
            */

            error = fmax(residual_p, residual_u);
            error = fmax(residual_v, error);

            //--- Convergence criteria ---
            if(itc != 1 && error < iterations.eps)
                break;

            //--- Update variables ---
            atualizar_matrizes_linearizadas<<<gridDimUm, blockDim>>>(dev_um_n_tau, dev_um_tau, imax+1, jmax, 1, jmax+1);
            atualizar_matrizes_linearizadas<<<gridDimVm, blockDim>>>(dev_vm_n_tau, dev_vm_tau, imax, jmax+1, 1, jmax+2);
            atualizar_matrizes_linearizadas<<<gridDim, blockDim>>>(dev_pn, dev_p, imax, jmax, 1, jmax+1);
            atualizar_matrizes_linearizadas<<<gridDim, blockDim>>>(dev_t_n_tau, dev_t_tau, imax, jmax, 1, jmax+1);
            atualizar_matrizes_linearizadas<<<gridDim, blockDim>>>(dev_c_n_tau, dev_c_tau, imax, jmax, 1, jmax+1);

            itc++;
        }

        //--- End of pseudo-time calculation ---
        atualizar_matrizes_linearizadas<<<gridDimUm, blockDim>>>(dev_um_n_tau, dev_um, imax+1, jmax, 1, jmax+1);
        atualizar_matrizes_linearizadas<<<gridDimVm, blockDim>>>(dev_vm_n_tau, dev_vm, imax, jmax+1, 1, jmax+2);
        atualizar_matrizes_linearizadas<<<gridDim, blockDim>>>(dev_t_n_tau, dev_t, imax, jmax, 1, jmax+1);
        atualizar_matrizes_linearizadas<<<gridDim, blockDim>>>(dev_c_n_tau, dev_c, imax, jmax, 1, jmax+1);

        /*--- Logs of time and intermediate results
        !IF (MOD(tr, n_tr) .EQ. 0) THEN
        !    WRITE(*,*) '-----------------------------------------------------------'
        !    WRITE(*,*) 'Max Residual:', error
        !    WRITE(*,*) 'Physical time:', time
        !    WRITE(*,*) '-----------------------------------------------------------'
        !    WRITE(*,*) '         dtau:',dtau   
        !    WRITE(*,*) '         dt:',dt   
        !    WRITE(*,*) '    Residual U:',residual_u
        !    WRITE(*,*) '    Residual V:',residual_v
        !    WRITE(*,*) '    Residual P:',residual_p    
        !    !WRITE(*,*) ' Artificial viscosity:',artMAX    
        !    !WRITE(*,*) ' Art Compressibility Par:',c2    
        !    WRITE(*,*) '-----------------------------------------------------------'

            !--- Output preliminary results ---
        !    CALL comp_mean(u, v, um, vm)
        !    CALL transient(u, v, p, T ,C, tr)

            !--- Output data file ---
        !    CALL output(um, vm, u, v, p, T, C, itc)
        !END IF
        */
        itc = 0;
        error = 100.0;

        // --- CUDA Events: fim do timestep ---
        cudaEventRecord(ev_ts_stop);
        cudaEventSynchronize(ev_ts_stop);
        cudaEventElapsedTime(&time_timestep, ev_ts_start, ev_ts_stop);

        // NOTA: Este código usa Unified Memory (cudaMallocManaged).
        // Não há cudaMemcpy explícito Host→Device ou Device→Host.
        // As transferências são gerenciadas implicitamente pelo driver CUDA
        // via page faults. H2D e D2H medem 0.0 pois não há cópias explícitas.
        // Para medir page migrations, use: nsys profile --trace=cuda,um
        time_h2d = 0.0f;
        time_d2h = 0.0f;

        printf("[ITER %d][CUDA] kernel_solve_U = %.6f ms\n", tr, acc_solveU);
        printf("[ITER %d][CUDA] kernel_solve_V = %.6f ms\n", tr, acc_solveV);
        printf("[ITER %d][CUDA] kernel_solve_P = %.6f ms\n", tr, acc_solveP);
        printf("[ITER %d][CUDA] kernel_solve_Z = %.6f ms\n", tr, acc_solveZ);
        printf("[ITER %d][CUDA] kernel_solve_C = %.6f ms\n", tr, acc_solveC);
        printf("[ITER %d][CUDA] H2D = %.6f ms\n", tr, time_h2d);
        printf("[ITER %d][CUDA] D2H = %.6f ms\n", tr, time_d2h);
        printf("[ITER %d][CUDA] timestep_total = %.6f ms\n", tr, time_timestep);

        tr = tr + 1;
    }
    //--- End of physical calculation ---

    /*duration = omp_get_wtime() - duration;
    printf("loop %lf\n", duration);
    duration = omp_get_wtime();
    */

    //--- Final results ---
    /*!open (550,file='data/time.dat')            
    !write (550,*) time
    !close(550)
    */
    //--- Compute the velocity of mean points ---
    comp_mean(dev_u, dev_v, dev_um, dev_vm);
    
    #ifdef DEBUG
    transient(dev_u, dev_v, dev_p, dev_t, dev_c, itc);
    #endif
    
    //--- output data file ---
    output(dev_um, dev_vm, dev_u, dev_v, dev_p, dev_t, dev_c, itc);

    /*duration = omp_get_wtime() - duration;
    printf("post %lf\n", duration);
    */

    //=================================================================
    // CUDA Events — Destruição
    //=================================================================
    cudaEventDestroy(ev_ts_start);    cudaEventDestroy(ev_ts_stop);
    cudaEventDestroy(ev_solveU_start); cudaEventDestroy(ev_solveU_stop);
    cudaEventDestroy(ev_solveV_start); cudaEventDestroy(ev_solveV_stop);
    cudaEventDestroy(ev_solveP_start); cudaEventDestroy(ev_solveP_stop);
    cudaEventDestroy(ev_solveZ_start); cudaEventDestroy(ev_solveZ_stop);
    cudaEventDestroy(ev_solveC_start); cudaEventDestroy(ev_solveC_stop);
    cudaEventDestroy(ev_h2d_start);   cudaEventDestroy(ev_h2d_stop);
    cudaEventDestroy(ev_d2h_start);   cudaEventDestroy(ev_d2h_stop);

    //desalocando
    cudaFree(dev_um);
    cudaFree(dev_um_n);
    cudaFree(dev_um_tau);
    cudaFree(dev_um_n_tau);
    cudaFree(dev_vm);
    cudaFree(dev_vm_n);
    cudaFree(dev_vm_tau);
    cudaFree(dev_vm_n_tau);
    cudaFree(dev_u);
    cudaFree(dev_v);
    cudaFree(dev_rc);
    cudaFree(dev_p);
    cudaFree(dev_pn);
    cudaFree(dev_h);
    cudaFree(dev_t);
    cudaFree(dev_z);
    cudaFree(dev_rz);
    cudaFree(dev_c);
    cudaFree(dev_t_n_tau);
    cudaFree(dev_t_tau);
    cudaFree(dev_c_n_tau);
    cudaFree(dev_rv);
    cudaFree(dev_ru);
    cudaFree(dev_ui);
    cudaFree(dev_res_u);
    cudaFree(dev_vi);
    cudaFree(dev_res_v);
    cudaFree(dev_c_tau);
    desalocar_globais();
}