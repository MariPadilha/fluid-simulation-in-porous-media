//Artificial Compressibility Methods
//Solve Momentum Equation with QUICK Scheme
//BASED ON:
//Versteeg, H. K., and W. Malalasekera. 
//"An introduction to computational Fluid Dynamics, The finite volume control, ed." (1995).
#include "comum.h"

int main(){
    calcular();         // define imax, jmax, dx_c, etc.
    alocar_globais(); 
    int itc, tr, i, j;

    double **um = (double**)malloc(sizeof(double*)*(imax+2));
    double **um_n = (double**)malloc(sizeof(double*)*(imax+2));
    double **res_u = (double**)malloc(sizeof(double*)*(imax+2));
    double **um_tau = (double**)malloc(sizeof(double*)*(imax+2));
    double **um_n_tau = (double**)malloc(sizeof(double*)*(imax+2));
    for(int i = 0; i < (imax+2); i++){
        um[i] = (double*)malloc(sizeof(double)*(jmax+1));
        um_n[i] = (double*)malloc(sizeof(double)*(jmax+1));
        res_u[i] = (double*)malloc(sizeof(double)*(jmax+1)); 
        um_tau[i] = (double*)malloc(sizeof(double)*(jmax+1));
        um_n_tau[i] = (double*)malloc(sizeof(double)*(jmax+1));
    }

    double **vm = (double**)malloc(sizeof(double*)*(imax+1));
    double **vm_n = (double**)malloc(sizeof(double*)*(imax+1));
    double **res_v = (double**)malloc(sizeof(double*)*(imax+1));
    double **vm_tau = (double**)malloc(sizeof(double*)*(imax+1));
    double **vm_n_tau = (double**)malloc(sizeof(double*)*(imax+1));
    for(int i = 0; i <= imax; i++){
        vm[i] = (double*)malloc(sizeof(double)*(jmax+2));
        vm_n[i] = (double*)malloc(sizeof(double)*(jmax+2));
        res_v[i] = (double*)malloc(sizeof(double)*(jmax+2));  
        vm_tau[i] = (double*)malloc(sizeof(double)*(jmax+2));
        vm_n_tau[i] = (double*)malloc(sizeof(double)*(jmax+2)); 
    }

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


    double residual_p, residual_u, residual_v, error, duration;

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

    for(j = 1; j <= jmax; j++){
        for(i = 1; i <= imax; i++){
            if(flag[i][j] != c_f){
                t[i][j] = temp_cylinder;
                c[i][j] = concentracao_inicial;
            }    
        }
    }

    //--- Set up initial flow field ---
    if(iterations.start_mode == 0){    
        IC(um, vm, p, t, c, pn);
    }else if(iterations.start_mode == 1){
        restart(um, vm, p, t, c);
    }else if(iterations.start_mode == 2){
        restart_dom(um, vm, p, t, z, h);
    }

    //--- Pseudo time step ---
    dtau = 5.e-2;
    dt = 0.5e-2;

    for(i = 1; i <= imax+1; i++){
        for(j = 1; j <= jmax; j++){
            um_tau[i][j] = um[i][j];
        }
    }
    for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax+1; j++){
            vm_tau[i][j] = vm[i][j];
        }
    }

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
            solve_U(um, vm, um_n, um_tau, vm_tau, um_n_tau, pn, &residual_u);
            solve_V(um, vm, vm_n, um_tau, vm_tau, vm_n_tau, pn, t, &residual_v);
            
            //--- Solve Continuity Equation ---
            solve_P(p, um_n_tau, vm_n_tau, pn, &residual_p);
            
            //--- Solve Energy Equation ---
            solve_Z(um_n_tau, vm_n_tau, t, t_n_tau, t_tau);
            solve_C(um_n_tau, vm_n_tau, c, c_n_tau, c_tau);

            /*--- check convergence ---
            CALL convergence(itc, error, residual_p, residual_u, residual_v)
            itc = itc+1
            */

            error = max(residual_p, residual_u);
            error = max(residual_v, error);

            //--- Convergence criteria ---
            if(itc != 1 && error < iterations.eps)
                break;

            //--- Update variables ---
            for(i = 1; i <= imax+1; i++){
                for(j = 1; j <= jmax; j++){
                    um_tau[i][j] = um_n_tau[i][j];
                }
            }

            for(i = 1; i <= imax; i++){
                for(j = 1; j <= jmax+1; j++){
                    vm_tau[i][j] = vm_n_tau[i][j];
                }
            }

            for(i = 1; i <= imax; i++){
                for(j = 1; j <= jmax; j++){
                    //printf("pn [%d][%d] = %lf\n", i, j, pn[i][j]);
                    p[i][j] = pn[i][j];
                    t_tau[i][j] = t_n_tau[i][j];
                    c_tau[i][j] = c_n_tau[i][j];
                }
            }

            itc += 1;
        }

        //--- End of pseudo-time calculation ---
        for(i = 1; i <= (imax+1); i++){
            for(j = 1; j <= jmax; j++){
                um[i][j] = um_n_tau[i][j];
                printf("[%i][%i] = %lf\n", i, j, um_n_tau[i][j]);
            }
        }
        
        for(i = 1; i <= imax; i++){
            for(j = 1; j <= jmax+1; j++){
                vm[i][j] = vm_n_tau[i][j];
            }
        }

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
    comp_mean(u, v, um, vm);
    transient(u, v, p, t, c, itc);

    //--- output data file ---
    output(um, vm, u, v, p, t, c, itc);

    /*duration = omp_get_wtime() - duration;
    printf("post %lf\n", duration);
    */
    for(int i = 0; i < (imax+2); i++){
        free(um[i]);
        free(um_n[i]);
        free(res_u[i]);
        free(um_tau[i]);
        free(um_n_tau[i]);
    }
    free(um);
    free(um_n);
    free(res_u);
    free(um_tau);
    free(um_n_tau);

    for(int i = 0; i <= imax; i++){
        free(vm[i]);
        free(vm_n[i]);
        free(res_v[i]);
        free(vm_tau[i]);
        free(vm_n_tau[i]);
    }
    free(vm);
    free(vm_n);
    free(res_v);
    free(vm_tau);
    free(vm_n_tau);

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