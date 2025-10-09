!Artificial Compressibility Methods
!Solve Momentum Equation with QUICK Scheme
!BASED ON:
!Versteeg, H. K., and W. Malalasekera. 
!"An introduction to computational Fluid Dynamics, The finite volume control, ed." (1995).


PROGRAM main
USE comum
USE omp_lib
IMPLICIT NONE

    INTEGER :: itc, tr, i, j
    INTEGER*4 today(3), now(3)   
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um, um_n, res_u
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm, vm_n, res_v
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau, um_n_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, vm_n_tau
    REAL(8), DIMENSION(1:imax,1:jmax)     :: u, v, P, Pn, H, T, Z, C
    REAL(8), DIMENSION(1:imax,1:jmax)     :: T_n_tau, T_tau, C_n_tau, C_tau 
    REAL(8) :: residual_p, residual_u, residual_v, error
    REAL(8) :: duration 
    duration = omp_get_wtime()
!    character(len=128) :: pwd
!    REAL ETIME, clockTIME, TARRAY(2)
!    clockTIME = ETIME(TARRAY)
!    CALL idate(today)   ! today(1)=day, (2)=month, (3)=year
!    CALL itime(now)     ! now(1)=hour, (2)=minute, (3)=second
!    CALL get_environment_variable('PWD', pwd)

    !--- Input data - Fill NAMELIST iterations and ref
    !--- Read iterations.dat and reference.dat
    !--- Compute Too, Tsub, Tinf 
    CALL init 

    Too = Ts * TnToo
    Tsup = Ts  / Ts
    Tinf = 1.0d0  ! Temperatura ambiente
 

    !--- Parameters ---
    !--- Compute cp_tot, rho_tot, k_tot, nu_tot, alpha_tot
    CALL properties 

    !--- Log de informações ---
!    5 FORMAT ( 'Date ', i2.2, '/', i2.2, '/', i4.4, &
!             '; time ',i2.2, ':', i2.2, ':', i2.2 )
!    WRITE(*,5)  today(2), today(1), today(3), now
!    WRITE(*,*) 'Current working directory: ',trim(pwd)
!    WRITE(*,*) '----------------------'
!    WRITE(*,*) 'Tsup =', Tsup
!    WRITE(*,*) 'Tinf =', Tinf
!    WRITE(*,*) '----------------------'
!    WRITE(*,*) 'Mesh =',imax,'x', jmax
!    WRITE(*,*) 'imax * jmax =', int(imax*jmax)
!    WRITE(*,*) 'eps =',eps
!    WRITE(*,*) '----------------------'
    !ver se tudo isso é necessário?
    !WRITE(*,*) '----------------------'
    !WRITE(*,*) 'v_c ='  ,v_c        , '[m/s]'
    !WRITE(*,*) 'L_c ='  ,L_c        , '[m]'
    !WRITE(*,*) 't_c ='  ,L_c/v_c    , '[s]'
    !WRITE(*,*) 'v_i ='  ,v_i
    !WRITE(*,*) 'V_idim ='  ,v_i*v_c , '[m/s]'
    !WRITE(*,*) '----------------------'
    !WRITE(*,*) 'g = '  ,g , '[m/s^2]'
    !WRITE(*,*) 'S = '  ,S
    !WRITE(*,*) 'q = '  ,q
    !WRITE(*,*) 'Pr ='  ,Pr
!    WRITE(*,*) 'Re =', Re
!    WRITE(*,*) 'Pe =', Pe
!    WRITE(*,*) 'Fr =', Fr
!    WRITE(*,*) 'InvFr^2 =', InvFr2
!    WRITE(*,*) '----------------------'

    !--- Initializations ---
    itc = 1 !Initial iteration
    error = 100.00d0
    tr = 1
    time = 0.d0
    residual_p = 0.d0

    !--- Create mesh ---
    CALL mesh

    DO j = 1, jmax
        DO i = 1, imax
            IF (flag(i,j) .NE. C_F) THEN
                T(i,j) = Temp_cylinder
                C(i,j) = concentracao_inicial
            END IF
        END DO
    END DO

    !--- Set up initial flow field ---
    IF (start_mode .EQ. 0) THEN    
        CALL IC(um, vm, p, T, C)
    ELSE IF (start_mode .EQ. 1) THEN
        CALL restart(um, vm, p, T, C)
    ELSE IF (start_mode .EQ. 2) THEN
        CALL restart_dom(um, vm, p, T, I, H)
    END IF

    !--- Pseudo time step ---
    dtau = 5.d-2
    dt = 0.5d-2

    um_tau = um
    vm_tau = vm

    duration = omp_get_wtime() - duration
    write(*,*) 'init', duration
    duration = omp_get_wtime()

    !--- Physical time step ---
    DO WHILE (time .LT. final_time)
        time = time + dt

        !--- Pseudo-time calculation starts ---
        DO WHILE (itc.LT.itc_max)
            !--- Solve Momentum Equation with QUICK Scheme ---
            CALL solve_U(um, vm, um_n, um_tau, vm_tau, um_n_tau, pn, residual_u)
            !write(*,*) 'Pos solve_U'
            CALL solve_V(um, vm, vm_n, um_tau, vm_tau, vm_n_tau, pn, T, residual_v)
            !write(*,*) 'Pos solve_V'
            !--- Solve Continuity Equation ---
            CALL solve_P(p, um_n_tau, vm_n_tau, pn, residual_p)
            !write(*,*) 'Pos solve_P'
            !--- Solve Energy Equation ---
            CALL solve_Z(um_n_tau, vm_n_tau, T, T_n_tau, T_tau)
            !write(*,*) 'Pos solve_Z'
            CALL solve_C(um_n_tau, vm_n_tau, C, C_n_tau, C_tau)
            !write(*,*) 'Pos solve_C'

            !--- check convergence ---
            !CALL convergence(itc, error, residual_p, residual_u, residual_v)
            !itc = itc+1

            error = MAX(residual_u, residual_v, residual_p)
            !write (*,*) itc, residual_u, residual_v, residual_p
            !--- Convergence criteria ---
            IF (itc .NE. 1 .AND. error .LT. eps) EXIT

            !--- Update variables ---
            um_tau = um_n_tau
            vm_tau = vm_n_tau
            p = pn
            T_tau = T_n_tau
            C_tau = C_n_tau
            itc = itc + 1
        END DO
        !--- End of pseudo-time calculation ---
        !write(*,*) 'End internal loop'
        um = um_n_tau
        vm = vm_n_tau
        T = T_n_tau
        C = C_n_tau

        !--- Logs of time and intermediate results
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

        itc = 0
        error = 100.d0

        tr = tr + 1
    END DO
    !--- End of physical calculation ---

    duration = omp_get_wtime() - duration
    write(*,*) 'loop', duration
    duration = omp_get_wtime()
    !--- Final results ---
    !open (550,file='data/time.dat')            
    !write (550,*) time
    !close(550)

    !--- Compute the velocity of mean points ---
    CALL comp_mean(u, v, um, vm)
            
    CALL transient(u, v, p, T, C, itc)

    !--- output data file ---
    CALL output(um, vm, u, v, p, T, C, itc)
    
    duration = omp_get_wtime() - duration
    write(*,*) 'post', duration
STOP
END PROGRAM main
