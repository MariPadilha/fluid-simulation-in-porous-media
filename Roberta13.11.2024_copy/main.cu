//Artificial Compressibility Methods
//Solve Momentum Equation with QUICK Scheme
//BASED ON:
//Versteeg, H. K., and W. Malalasekera. 
//"An introduction to computational Fluid Dynamics, The finite volume control, ed." (1995).
#include "comum.h"
#define N_IMAX 10
#define N_ITC 2

__global__ void atualizar_matrizes_linearizadas(double *origem, double *destino, int tamanhoLinha, int tamanhoColuna, int inicio, int coluna){
    int i = blockIdx.x * blockDim.x + threadIdx.x + inicio;
    int j = blockIdx.y * blockDim.y + threadIdx.y + inicio;
    
    if(i <= tamanhoLinha && j <= tamanhoColuna){
        destino[i*(coluna)+j] = origem[i*(coluna)+j];
    }
}

int main(int argc, char *argv[]){
    int itc, tr, i, j, threads = 256;
    int n_imax, n_itc;
    dim3 blocks;

    if(argc < 3){
        n_imax = N_IMAX;
        n_itc = N_ITC;
    }else{
        n_imax = atoi(argv[1]);
        n_itc = atoi(argv[2]);
    }

    calcular(n_imax, n_itc);         // define imax, jmax, dx_c, etc.
    alocar_globais(); 

    double *dev_um, *dev_vm;
    double *dev_um_n, *dev_vm_n;
    double *dev_um_tau, *dev_vm_tau;
    double *dev_um_n_tau, *dev_vm_n_tau;
    cudaMallocManaged((void**)&dev_um, sizeof(double)*(imax+2)*(jmax+1));
    cudaMallocManaged((void**)&dev_vm, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_um_n, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_um_tau, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_um_n_tau, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_vm_n, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_vm_tau, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_vm_n_tau, sizeof(double)*(imax+1)*(jmax+2));
    

    double **u = (double**)malloc(sizeof(double*)*(imax+1));
    double **v = (double**)malloc(sizeof(double*)*(imax+1));
    double **p = (double**)malloc(sizeof(double*)*(imax+1));
    double **pn = (double**)malloc(sizeof(double*)*(imax+1));
    double **h = (double**)malloc(sizeof(double*)*(imax+1));
    double **t = (double**)malloc(sizeof(double*)*(imax+1));
    double **z = (double**)malloc(sizeof(double*)*(imax+1));
    double **c = (double**)malloc(sizeof(double*)*(imax+1));
    double **t_n_tau = (double**)malloc(sizeof(double*)*(imax+1));
    double **t_tau = (double**)malloc(sizeof(double*)*(imax+1));
    double **c_n_tau = (double**)malloc(sizeof(double*)*(imax+1));
    double **c_tau = (double**)malloc(sizeof(double*)*(imax+1));
    for(int i = 0; i < imax+1; i++){
        u[i] = (double*)malloc(sizeof(double)*(jmax+1));
        v[i] = (double*)malloc(sizeof(double)*(jmax+1));
        p[i] = (double*)malloc(sizeof(double)*(jmax+1));
        pn[i] = (double*)malloc(sizeof(double)*(jmax+1));
        h[i] = (double*)malloc(sizeof(double)*(jmax+1));
        t[i] = (double*)malloc(sizeof(double)*(jmax+1));
        z[i] = (double*)malloc(sizeof(double)*(jmax+1));
        c[i] = (double*)malloc(sizeof(double)*(jmax+1));
        t_n_tau[i] = (double*)malloc(sizeof(double)*(jmax+1));
        t_tau[i] = (double*)malloc(sizeof(double)*(jmax+1));
        c_n_tau[i] = (double*)malloc(sizeof(double)*(jmax+1));
        c_tau[i] = (double*)malloc(sizeof(double)*(jmax+1)); 
    }


    double residual_p, residual_u, residual_v, error;

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
    time = 0.0;
    residual_p = 0.0;

    //--- Create mesh ---
    mesh();

    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){
            if(flag[i][j] != c_f){
                t[i][j] = temp_cylinder;
                c[i][j] = concentracao_inicial;
            }    
        }
    }

    //--- Set up initial flow field ---
    if(iterations.start_mode == 0){    
        IC(dev_um, dev_vm, p, t, c, pn);
    }else if(iterations.start_mode == 1){
        restart(dev_um, dev_vm, p, t, c);
    }else if(iterations.start_mode == 2){
        restart_dom(dev_um, dev_vm, p, t, z, h);
    }

    //--- Pseudo time step ---
    dtau = 5.e-2;
    dt = 0.5e-2;

    blocks = grid_1d(((imax+1)*jmax), threads);
    atualizar_matrizes_linearizadas<<<blocks, threads>>>(dev_um, dev_um_tau, imax+1, jmax, 1, jmax+1);

    blocks = grid_1d((imax*(jmax+1)), threads);
    atualizar_matrizes_linearizadas<<<blocks, threads>>>(dev_vm, dev_vm_tau, imax, jmax+1, 1, jmax+2);


    /*duration = omp_get_wtime() - duration;
    printf("init %lf/n", duration);
    duration = omp_get_wtime();
    */

    //--- Physical time step ---
    while(time < iterations.final_time){
        time = time + dt;

        //--- Pseudo-time calculation starts ---
        while(itc < iterations.itc_max){
            //--- Solve Momentum Equation with QUICK Scheme ---
            solve_U(dev_um, dev_vm, dev_um_n, dev_um_tau, dev_vm_tau, dev_um_n_tau, pn, &residual_u);
            solve_V(dev_um, dev_vm, dev_vm_n, dev_um_tau, dev_vm_tau, dev_vm_n_tau, pn, t, &residual_v);
            
            //--- Solve Continuity Equation ---
            solve_P(p, dev_um_n_tau, dev_vm_n_tau, pn, &residual_p);
            
            //--- Solve Energy Equation ---
            solve_Z(dev_um_n_tau, dev_vm_n_tau, t, t_n_tau, t_tau);
            solve_C(dev_um_n_tau, dev_vm_n_tau, c, c_n_tau, c_tau);

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
            blocks = grid_1d(((imax+1)*jmax), threads);
            atualizar_matrizes_linearizadas<<<blocks, threads>>>(dev_um_n_tau, dev_um_tau, imax+1, jmax, 1, jmax+1);

            blocks = grid_1d((imax*(jmax+1)), threads);
            atualizar_matrizes_linearizadas<<<blocks, threads>>>(dev_vm_n_tau, dev_vm_tau, imax, jmax+1, 1, jmax+2);
   

            for(i = 1; i <= imax; i++){
                for(j = 1; j <= jmax; j++){
                    p[i][j] = pn[i][j];
                    t_tau[i][j] = t_n_tau[i][j];
                    c_tau[i][j] = c_n_tau[i][j];
                }
            }

            itc++;
        }

        //--- End of pseudo-time calculation ---
        blocks = grid_1d(((imax+1)*jmax), threads);
        atualizar_matrizes_linearizadas<<<blocks, threads>>>(dev_um_n_tau, dev_um, imax+1, jmax, 1, jmax+1);

        blocks = grid_1d((imax*(jmax+1)), threads);
        atualizar_matrizes_linearizadas<<<blocks, threads>>>(dev_vm_n_tau, dev_vm, imax, jmax+1, 1, jmax+2);

        for(i = 1; i <= imax; i++){
            for(j = 1; j <= jmax; j++){
                t[i][j] = t_n_tau[i][j];
                c[i][j] = c_n_tau[i][j];
            }
        }
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
    comp_mean(u, v, dev_um, dev_vm);
    transient(u, v, p, t, c, itc);

    //--- output data file ---
    output(dev_um, dev_vm, u, v, p, t, c, itc);

    /*duration = omp_get_wtime() - duration;
    printf("post %lf\n", duration);
    */

    //desalocando
    cudaFree(dev_um);
    cudaFree(dev_um_n);
    cudaFree(dev_um_tau);
    cudaFree(dev_um_n_tau);
    cudaFree(dev_vm);
    cudaFree(dev_vm_n);
    cudaFree(dev_vm_tau);
    cudaFree(dev_vm_n_tau);

    for(int i = 0; i < imax+1; i++){
        free(u[i]);
        free(v[i]);
        free(p[i]);
        free(pn[i]);
        free(h[i]);
        free(t[i]);
        free(z[i]);
        free(c[i]);
        free(t_n_tau[i]);
        free(t_tau[i]);
        free(c_n_tau[i]);
        free(c_tau[i]);
    }
    free(u);
    free(v);
    free(p);
    free(pn);
    free(h);
    free(t);
    free(z);
    free(c);
    free(t_n_tau);
    free(t_tau);
    free(c_n_tau);
    free(c_tau);
    desalocar_globais();
}