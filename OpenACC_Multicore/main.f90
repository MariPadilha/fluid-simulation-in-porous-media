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
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau, um_n_tau   !sao os temps de um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, vm_n_tau   !sao os temps de vm
    REAL(8), DIMENSION(1:imax,1:jmax)     :: u, v, P, T, Z, C
    REAL(8), DIMENSION(1:imax,1:jmax)     :: T_n_tau, C_n_tau 
    REAL(8) :: residual_p, residual_u, residual_v, error
    REAL(8) :: duration 
    duration = omp_get_wtime()

    !--- Input data / propriedades / malha ---
    CALL init 
    Too = Ts * TnToo
    Tsup = Ts / Ts
    Tinf = 1.0d0
    CALL properties 
    itc = 1
    error = 100.00d0
    tr = 1
    time = 0.d0
    residual_p = 0.d0
    CALL mesh

    IF (start_mode .EQ. 0) THEN    
        CALL IC(um, vm, p, T, C)
    ELSE IF (start_mode .EQ. 1) THEN
        CALL restart(um, vm, p, T, C)
    ELSE IF (start_mode .EQ. 2) THEN
        CALL restart_dom(um, vm, p, T, I)
    END IF

    dtau = 5.d-2
    dt = 0.5d-2
    um_tau = um
    vm_tau = vm

    !==============================================================
    ! OpenACC: alocação/cópias NA GPU ANTES do laço iterativo
    !==============================================================
    !inicializados na mesh.f90. Todos sao globais declarados em comum. dx e dy só usados em mesh.f90
!    !$acc enter data copyin(x,y,xm,ym,vol_u,vol_v,vol_p,areau_n,areau_s,areau_e,areau_w,areav_n,areav_s,areav_e,areav_w,liga_poros,epsilon1)
    !!flag

    !Temporários usados nas subrotinas de equations.f90
!    !$acc enter data create(fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap)
!    !$acc enter data create(Dn, Ds, De, Dw, Dp)
!    !$acc enter data create(u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P)
!    !$acc enter data create(v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P)
!    !$acc enter data create(afw, afe, afn, afs)
!    !$acc enter data create(artDivU, artDivV, q_art)
!    !$acc enter data create(dudxdx, dvdydy, dxdvdy, dydudx)
  
    !Temporários de solve_U
!    !$acc enter data create(RU, Ui, res_U)
!    !$acc enter data create(RV, Vi, res_V)
    
!    !$acc enter data copyin(um, vm, um_tau, vm_tau, p, u, v, T, Z, C)
!    !$acc enter data create(um_n_tau, vm_n_tau, T_n_tau, C_n_tau)
   
    !Temporários usados só em solve_P
!    !$acc enter data create(RP, Pi, res_P)
!    !$acc enter data create(dudx, dvdy)

!    !$acc enter data create(RZ, Zi, res_Z)
!    !$acc enter data create(dZudx, dZvdy)

!    !$acc enter data create(RC, Ci, res_C)
!    !$acc enter data create(dCudx, dCvdy)

    duration = omp_get_wtime() - duration
    write(*,*) 'init', duration
    duration = omp_get_wtime()

    !--- Physical time step ---
    DO WHILE (time .LT. final_time)
        time = time + dt

        !--- Pseudo-time calculation starts ---
        DO WHILE (itc.LT.itc_max)

            CALL solve_U(um, vm, um_tau, vm_tau, um_n_tau, p)

            CALL solve_V(um, vm, um_tau, vm_tau, vm_n_tau, p, T)

            CALL solve_P(um_n_tau, vm_n_tau, p)

            CALL solve_Z(um_n_tau, vm_n_tau, T, T_n_tau)!, T_tau)

            CALL solve_C(um_n_tau, vm_n_tau, C, C_n_tau)!, C_tau)

!            !$acc update self(res_u, res_v, res_p)
            error = MAX( MAXVAL(ABS(res_u)), MAXVAL(ABS(res_v)), MAXVAL(ABS(res_p)) )
            !error = MAX(residual_u, residual_v, residual_p)
            !write(*,*) res_u !, res_v, res_p
            !write(*,*) itc, MAXVAL(ABS(res_u)), MAXVAL(ABS(res_v)), MAXVAL(ABS(res_p))

            IF (itc .NE. 1 .AND. error .LT. eps) EXIT
            !IF(itc .EQ. 1) STOP            

            !um_tau = um_n_tau
            !$acc parallel loop collapse(2) default(present)
            do j=1,jmax
              do i=1,imax+1
                um_tau(i,j) = um_n_tau(i,j)
              enddo
            end do

            !vm_tau = vm_n_tau
            !$acc parallel loop collapse(2) default(present)
            do j=1,jmax+1
              do i=1,imax
                vm_tau(i,j) = vm_n_tau(i,j)
              enddo
            end do
            
            !p = pn !update no proprio solve_P
            !T_tau = T_n_tau
            !C_tau = C_n_tau
            itc = itc + 1
        END DO
        !--- End of pseudo-time calculation ---
        !um = um_n_tau
        !$acc parallel loop collapse(2) default(present)
            do j=1,jmax
              do i=1,imax+1
                um(i,j) = um_n_tau(i,j)
              enddo
            end do

        !vm = vm_n_tau
        !$acc parallel loop collapse(2) default(present)
            do j=1,jmax+1
              do i=1,imax
                vm(i,j) = vm_n_tau(i,j)
              enddo
            end do
        !T = T_n_tau
        !$acc parallel loop collapse(2) default(present)
            do j=1,jmax
              do i=1,imax
                T(i,j) = T_n_tau(i,j)
              enddo
            end do

        !C = C_n_tau
        !$acc parallel loop collapse(2) default(present)
            do j=1,jmax
              do i=1,imax
                C(i,j) = C_n_tau(i,j)
              enddo
            end do

        itc = 0
        error = 100.d0
        tr = tr + 1
    END DO
    !--- End of physical calculation ---

    duration = omp_get_wtime() - duration
    write(*,*) 'loop', duration
    duration = omp_get_wtime()

    !==============================================================
    ! Copiar resultados de volta para o host e desalocar na GPU
    !==============================================================
!    !$acc update self(um, vm, u, v, p, T, C)
    ! $acc exit data delete(flag)
!    !$acc exit data delete(res_u, res_v, um_n_tau, vm_n_tau, T_n_tau, C_n_tau)
!    !$acc exit data delete(um, vm, um_tau, vm_tau, p, u, v, T, Z, C)

    CALL comp_mean(u, v, um, vm)        
    CALL transient(u, v, p, T, C, itc)
    CALL output(um, vm, u, v, p, T, C, itc)

    duration = omp_get_wtime() - duration
    write(*,*) 'post', duration
STOP
END PROGRAM main
