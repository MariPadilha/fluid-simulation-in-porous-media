!==============================================================
! OpenACC: este arquivo pressupõe que TODOS os dados necessários
! já estão PRESENTES na GPU (alocados e copiados no main.f90).
! Evite data clauses aqui; use default(present) nos loops.
!==============================================================
!--- ResU ---
subroutine RESU(um,vm,p)
    use comum
    implicit none
    integer :: i, j   
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,  1:jmax  ) :: P

    !$acc parallel loop collapse(2) default(present)
    do j = 3, jmax-2
      do i = 3, imax-1
        fn(i,j) = 0.5d0 * ( vm(i  ,j+1) + vm(i-1,j+1) ) * areau_n(i) / (epsilon1(i,j) )
        fs(i,j) = 0.5d0 * ( vm(i  ,j  ) + vm(i-1,j  ) ) * areau_s(i) / (epsilon1(i,j) )
        fe(i,j) = 0.5d0 * ( um(i+1,j  ) + um(i  ,j  ) ) * areau_e(j) / (epsilon1(i,j) )
        fw(i,j) = 0.5d0 * ( um(i  ,j  ) + um(i-1,j  ) ) * areau_w(j) / (epsilon1(i,j) )
        df(i,j) = fe(i,j) - fw(i,j) + fn(i,j) - fs(i,j)
    
        Dn(i,j) = (epsilon1(i,j)/Re) * areau_n(i) / (y(j+1)-y(j  ))
        Ds(i,j) = (epsilon1(i,j)/Re) * areau_s(i) / (y(j  )-y(j-1))
        De(i,j) = (epsilon1(i,j)/Re) * areau_e(j) / (xm(i+1)-xm(i  ))
        Dw(i,j) = (epsilon1(i,j)/Re) * areau_w(j) / (xm(i  )-xm(i-1))
    
        !quick
        if(fw(i,j).GT.0.0d0) then
            afw(i,j) = 1.0d0
        elseif(fw(i,j).LT.0.0d0) then
            afw(i,j) = 0.0d0
        endif

        if(fe(i,j).GT.0.0d0) then
                afe(i,j) = 1.0d0
        elseif(fe(i,j).LT.0.0d0) then
                afe(i,j) = 0.0d0
        endif

        if(fn(i,j).GT.0.0d0) then
                afn(i,j) = 1.0d0
        elseif(fn(i,j).LT.0.0d0) then
                afn(i,j) = 0.0d0
        endif

        if(fs(i,j).GT.0.0d0) then
                afs(i,j) = 1.0d0
        elseif(fs(i,j).LT.0.0d0) then
                afs(i,j) = 0.0d0
        endif

        aw(i,j) = Dw(i,j) + 0.75d0  * afw(i,j) * fw(i,j) &
                + 0.125d0 * afe(i,j) * fe(i,j) &
                + 0.375d0 * ( 1.0d0 - afw(i,j) ) * fw(i,j)

        ae(i,j) =  De(i,j) - 0.375d0* afe(i,j) * fe(i,j) &
                - 0.75d0  * ( 1.0d0 - afe(i,j) ) * fe(i,j) &
                - 0.125d0 * ( 1.0d0 - afw(i,j) ) * fw(i,j)

        as(i,j) = Ds(i,j) + 0.75d0  * afs(i,j) * fs(i,j) &
                + 0.125d0 * afn(i,j) * fn(i,j) &
                + 0.375d0 * ( 1.0d0 - afs(i,j) ) * fs(i,j)

        an(i,j) =  Dn(i,j) - 0.375d0* afn(i,j) * fn(i,j) &
                - 0.75d0  * ( 1.0d0 - afn(i,j) ) * fn(i,j) &
                - 0.125d0 * ( 1.0d0 - afs(i,j) ) * fs(i,j)

        aww(i,j) = -0.125d0 *           afw(i,j)   * fw(i,j)
        aee(i,j) =  0.125d0 * ( 1.0d0 - afe(i,j) ) * fe(i,j)
        ass(i,j) = -0.125d0 *           afs(i,j)   * fs(i,j)
        ann(i,j) =  0.125d0 * ( 1.0d0 - afn(i,j) ) * fn(i,j)
        ap(i,j) = aw(i,j) + ae(i,j) + as(i,j) + an(i,j) + aww(i,j) + aee(i,j) + ass(i,j) + ann(i,j) + df(i,j)
        !end Quick

        u_W(i,j)  = um(i-1,j  )
        u_WW(i,j) = um(i-2,j  )
        u_E(i,j)  = um(i+1,j  )
        u_EE(i,j) = um(i+2,j  )
        u_S(i,j)  = um(i  ,j-1)
        u_SS(i,j) = um(i  ,j-2)
        u_N(i,j)  = um(i  ,j+1)
        u_NN(i,j) = um(i  ,j+2)        
        u_P(i,j)  = um(i  ,j  )
        v_P(i,j)  = vm(i  ,j  )

        dudxdx(i,j) = areau_e(j) * ( u_E(i,j) - u_P(i,j) ) / (xm(i+1)-xm(i  )) &
                -areau_w(j) * ( u_P(i,j) - u_W(i,j) ) / (xm(i  )-xm(i-1)) 
        dxdvdy(i,j) = areau_e(j) * (vm(i  ,j+1) - vm(i  ,j  )) / (ym(j+1)-ym(j)) &
                -areau_w(j) * (vm(i-1,j+1) - vm(i-1,j  )) / (ym(j+1)-ym(j))
            
        artDivU(i,j) = - b_art * (dudxdx(i,j) + dxdvdy(i,j))   
            
        !bulk artificial viscosity term from Ramshaw(1990)
        q_art(i,j) = epsilon1(i,j) * (p(i,j)-p(i-1,j)) / (x(i)-x(i-1)) + artDivU(i,j)
    
        RU(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * (-ap(i,j) * u_P(i,j) &
                + aww(i,j) * u_WW(i,j) + aw(i,j) * u_W(i,j)  &
                + aee(i,j) * u_EE(i,j) + ae(i,j) * u_E(i,j)  &
                + ass(i,j) * u_SS(i,j) + as(i,j) * u_S(i,j)  &
                + ann(i,j) * u_NN(i,j) + an(i,j) * u_N(i,j)) &
                - q_art(i,j) - epsilon1(i,j) * ( u_P(i,j)/(Re*Darcy_number) + &
                Cf/((epsilon1(i,j) * Darcy_number)**0.5d0) * u_p(i,j) * &
                ((u_p(i,j)**2.d0 + v_p(i,j)**2.d0)**0.5d0))*liga_poros(i,j) -g*epsilon1(i,j)
      enddo
    enddo

    CALL upwind_Ui(um,vm,p,2)
    CALL upwind_Ui(um,vm,p,jmax-1)

    CALL upwind_Uj(um,vm,p,2)
    CALL upwind_Uj(um,vm,p,imax)

return
end subroutine RESU


!--- upwind_U ---
subroutine upwind_Ui(um,vm,p,j)
    use comum
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um !, RU
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,  1:jmax  ) :: P

    !$acc parallel loop default(present)
    do i=2,imax
        !compute x-direction velocity component un
        fn(i,j) = 0.5d0 * ( vm(i  ,j+1) + vm(i-1,j+1) ) * areau_n(i) / (epsilon1(i,j) )
        fs(i,j) = 0.5d0 * ( vm(i  ,j  ) + vm(i-1,j  ) ) * areau_s(i) / (epsilon1(i,j) )
        fe(i,j) = 0.5d0 * ( um(i+1,j  ) + um(i  ,j  ) ) * areau_e(j) / (epsilon1(i,j) )
        fw(i,j) = 0.5d0 * ( um(i  ,j  ) + um(i-1,j  ) ) * areau_w(j) / (epsilon1(i,j) )

        df(i,j) = fe(i,j) - fw(i,j) + fn(i,j) - fs(i,j)

        Dn(i,j) = (epsilon1(i,j)/Re) * areau_n(i) / (y(j+1)-y(j  ))
        Ds(i,j) = (epsilon1(i,j)/Re) * areau_s(i) / (y(j  )-y(j-1))
        De(i,j) = (epsilon1(i,j)/Re) * areau_e(j) / (xm(i+1)-xm(i  ))
        Dw(i,j) = (epsilon1(i,j)/Re) * areau_w(j) / (xm(i  )-xm(i-1))

        !upwind
        aw(i,j) = Dw(i,j) + MAX(fw(i,j) , 0.0d0)
        as(i,j) = Ds(i,j) + MAX(fs(i,j) , 0.0d0)
        ae(i,j) = De(i,j) + MAX(0.0d0 , -fe(i,j))
        an(i,j) = Dn(i,j) + MAX(0.0d0 , -fn(i,j))

        ap(i,j) = aw(i,j) + ae(i,j) + as(i,j) + an(i,j) + df(i,j)

        u_W(i,j) = um(i-1,j  )
        u_E(i,j) = um(i+1,j  )
        u_S(i,j) = um(i  ,j-1)
        u_N(i,j) = um(i  ,j+1)
        u_P(i,j) = um(i  ,j  )
        v_P(i,j) = vm(i  ,j  )

        dudxdx(i,j) = areau_e(j) * ( u_E(i,j) - u_P(i,j) ) / (xm(i+1)-xm(i  )) &
                     -areau_w(j) * ( u_P(i,j) - u_W(i,j) ) / (xm(i  )-xm(i-1)) 
    
        dxdvdy(i,j) = areau_e(j) * (vm(i  ,j+1) - vm(i  ,j)) / (ym(j+1)-ym(j)) &
                     -areau_w(j) * (vm(i-1,j+1) - vm(i-1,j)) / (ym(j+1)-ym(j))

        !bulk artificial viscosity term from Ramshaw(1990)
        q_art(i,j) = epsilon1(i,j) * ( p(i,j)-p(i-1,j) ) / (x(i)-x(i-1)) - b_art * (dudxdx(i,j) + dxdvdy(i,j))

        RU(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - ap(i,j) * u_P(i,j) &
                + aw(i,j) * u_W(i,j) + ae(i,j) * u_E(i,j)  &
                + as(i,j) * u_S(i,j) + an(i,j) * u_N(i,j)) &
                - q_art(i,j) - epsilon1(i,j) * (  u_P(i,j)/(Re*Darcy_number) +&
                Cf/((epsilon1(i,j)*Darcy_number)**0.5d0) * u_P(i,j)*&
                ((u_P(i,j)**2.d0 + v_P(i,j)**2.d0)**0.5d0))*liga_poros(i,j) - g*epsilon1(i,j)
  enddo

return
end subroutine upwind_Ui


!--- upwind_U ---
subroutine upwind_Uj(um,vm,p,i)
    use comum
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um !, RU
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,1:jmax) :: P

    !$acc parallel loop default(present)
    do j=2,jmax-1
        !compute x-direction velocity component un
        fn(i,j) = 0.5d0 * ( vm(i  ,j+1) + vm(i-1,j+1) ) * areau_n(i) / (epsilon1(i,j) )
        fs(i,j) = 0.5d0 * ( vm(i  ,j  ) + vm(i-1,j  ) ) * areau_s(i) / (epsilon1(i,j) )
        fe(i,j) = 0.5d0 * ( um(i+1,j  ) + um(i  ,j  ) ) * areau_e(j) / (epsilon1(i,j) )
        fw(i,j) = 0.5d0 * ( um(i  ,j  ) + um(i-1,j  ) ) * areau_w(j) / (epsilon1(i,j) )

        df(i,j) = fe(i,j) - fw(i,j) + fn(i,j) - fs(i,j)

        Dn(i,j) = (epsilon1(i,j)/Re) * areau_n(i) / (y(j+1)-y(j  ))
        Ds(i,j) = (epsilon1(i,j)/Re) * areau_s(i) / (y(j  )-y(j-1))
        De(i,j) = (epsilon1(i,j)/Re) * areau_e(j) / (xm(i+1)-xm(i  ))
        Dw(i,j) = (epsilon1(i,j)/Re) * areau_w(j) / (xm(i  )-xm(i-1))

        !upwind
        aw(i,j) = Dw(i,j) + MAX(fw(i,j) , 0.0d0)
        as(i,j) = Ds(i,j) + MAX(fs(i,j) , 0.0d0)
        ae(i,j) = De(i,j) + MAX(0.0d0 , -fe(i,j))
        an(i,j) = Dn(i,j) + MAX(0.0d0 , -fn(i,j))

        ap(i,j) = aw(i,j) + ae(i,j) + as(i,j) + an(i,j) + df(i,j)

        u_W(i,j) = um(i-1,j  )
        u_E(i,j) = um(i+1,j  )
        u_S(i,j) = um(i  ,j-1)
        u_N(i,j) = um(i  ,j+1)
        u_P(i,j) = um(i  ,j  )
        v_P(i,j) = vm(i  ,j  )

        dudxdx(i,j) = areau_e(j) * ( u_E(i,j) - u_P(i,j) ) / (xm(i+1)-xm(i  )) &
                     -areau_w(j) * ( u_P(i,j) - u_W(i,j) ) / (xm(i  )-xm(i-1)) 
    
        dxdvdy(i,j) = areau_e(j) * (vm(i  ,j+1) - vm(i  ,j)) / (ym(j+1)-ym(j)) &
                     -areau_w(j) * (vm(i-1,j+1) - vm(i-1,j)) / (ym(j+1)-ym(j))

        !bulk artificial viscosity term from Ramshaw(1990)
        q_art(i,j) = epsilon1(i,j) * ( p(i,j)-p(i-1,j) ) / (x(i)-x(i-1)) - b_art * (dudxdx(i,j) + dxdvdy(i,j))

        RU(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - ap(i,j) * u_P(i,j) &
                + aw(i,j) * u_W(i,j) + ae(i,j) * u_E(i,j)  &
                + as(i,j) * u_S(i,j) + an(i,j) * u_N(i,j)) &
                - q_art(i,j) - epsilon1(i,j) * (  u_P(i,j)/(Re*Darcy_number) +&
                Cf/((epsilon1(i,j)*Darcy_number)**0.5d0) * u_P(i,j)*&
                ((u_P(i,j)**2.d0 + v_P(i,j)**2.d0)**0.5d0))*liga_poros(i,j) - g*epsilon1(i,j)
    enddo

return
end subroutine upwind_Uj


!--- solve_U ---
subroutine solve_U(um,vm,um_tau,vm_tau,um_n_tau,p)
    use comum
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau, um_n_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax  ) :: p
   
    CALL RESU(um_tau,vm_tau,p)
!$acc update self(RU)
!write (*,*) RU 
    !$acc parallel loop collapse(2) default(present)
    do j=2,jmax-1
        do i=3,imax-1
            res_u(i,j) = ( (um(i,j)-um_tau(i,j)) + RU(i,j)*dt) * dtau
            ui(i,j) = um_tau(i,j) + res_u(i,j)
        end do
    end do

    !call bcUV(ui,vm_tau)
    !------ contorno inferior -----------
    !$acc parallel loop default(present)
    DO i=1,imax
        vm_tau(i,1) = 0.5d0 !v_i
        vm_tau(i,2) = 0.5d0 !v_i !vm(i,1)     
    ENDDO
    !$acc parallel loop default(present)
    DO i=2,imax
        ui(i,1) = 0.d0    !!!DARCY
    ENDDO
    !------ contorno superior ----------
    !$acc parallel loop default(present)
    DO i=2,imax
        vm_tau(i,jmax  ) = 2.d0*vm_tau(i,jmax-2) - vm_tau(i,jmax-1)
        vm_tau(i,jmax+1) = 2.d0*vm_tau(i,jmax-1) - vm_tau(i,jmax)
    ENDDO
    !$acc parallel loop default(present)
    DO i=2,imax
        ui(i,jmax) = ui(i,jmax-1)
    ENDDO
    !------ contorno esquerdo --------
    !$acc parallel loop default(present)
    DO j=1,jmax
        ui(2,j) = 0.d0 !symmetry
        ui(1,j) = 0.d0 !symmetry
    ENDDO
    !$acc parallel loop default(present)
    DO j=1,jmax+1
        vm_tau(1,j) = vm_tau(2,j) !darcy
    ENDDO
    !------ contorno direito --------
    !$acc parallel loop default(present)
    DO j=1,jmax
        ui(imax  ,j) = 0.d0
        ui(imax+1,j) = 0.d0
    ENDDO
    !$acc parallel loop default(present)
    DO j=2,jmax
        vm_tau(imax,j) =  vm_tau(imax-1,j) !Darcy
    ENDDO
    !------------------------------------

    CALL RESU(ui,vm_tau,p)

    !$acc parallel loop collapse(2) default(present)
    do j=2,jmax-1
        do i=3,imax-1
            res_u(i,j) = ( (um(i,j)-um_tau(i,j)) + RU(i,j)*dt) * dtau
            ui(i,j) = 0.75d0 * um_tau(i,j) + 0.25d0 * (ui(i,j) + res_u(i,j))            
        enddo
    enddo

    !call bcUV(ui,vm_tau)
    !------ contorno inferior -----------
    !$acc parallel loop default(present)
    DO i=1,imax
        vm_tau(i,1) = 0.5d0 !v_i
        vm_tau(i,2) = 0.5d0 !v_i !vm(i,1)     
    ENDDO
    !$acc parallel loop default(present)
    DO i=2,imax
        ui(i,1) = 0.d0    !!!DARCY
    ENDDO
    !------ contorno superior ----------
    !$acc parallel loop default(present)
    DO i=2,imax
        vm_tau(i,jmax  ) = 2.d0*vm_tau(i,jmax-2) - vm_tau(i,jmax-1)
        vm_tau(i,jmax+1) = 2.d0*vm_tau(i,jmax-1) - vm_tau(i,jmax)
    ENDDO
    !$acc parallel loop default(present)
    DO i=2,imax
        ui(i,jmax) = ui(i,jmax-1)
    ENDDO
    !------ contorno esquerdo --------
    !$acc parallel loop default(present)
    DO j=1,jmax
        ui(2,j) = 0.d0 !symmetry
        ui(1,j) = 0.d0 !symmetry
    ENDDO
    !$acc parallel loop default(present)
    DO j=1,jmax+1
        vm_tau(1,j) = vm_tau(2,j) !darcy
    ENDDO
    !------ contorno direito --------
    !$acc parallel loop default(present)
    DO j=1,jmax
        ui(imax  ,j) = 0.d0
        ui(imax+1,j) = 0.d0
    ENDDO
    !$acc parallel loop default(present)
    DO j=2,jmax
        vm_tau(imax,j) =  vm_tau(imax-1,j) !Darcy
    ENDDO
    !------------------------------------

    CALL RESU(ui,vm_tau,p)

    !$acc parallel loop collapse(2) default(present)
    do j=2,jmax-1
        do i=3,imax-1
            res_u(i,j) = ( (um(i,j)-um_tau(i,j)) + RU(i,j)*dt) * dtau
            um_n_tau(i,j) = 1.0d0 / 3.0d0 * um_tau(i,j) + 2.0d0 / 3.0d0 * (ui(i,j) + res_u(i,j)) 
        enddo
    enddo

    !call bcUV(um_n_tau,vm_tau)
     !------ contorno inferior -----------
    !$acc parallel loop default(present)
    DO i=1,imax
        vm_tau(i,1) = 0.5d0 !v_i
        vm_tau(i,2) = 0.5d0 !v_i !vm(i,1)     
    ENDDO
    !$acc parallel loop default(present)
    DO i=2,imax
        um_n_tau(i,1) = 0.d0    !!!DARCY
    ENDDO
    !------ contorno superior ----------
    !$acc parallel loop default(present)
    DO i=2,imax
        vm_tau(i,jmax  ) = 2.d0*vm_tau(i,jmax-2) - vm_tau(i,jmax-1)
        vm_tau(i,jmax+1) = 2.d0*vm_tau(i,jmax-1) - vm_tau(i,jmax)
    ENDDO
    !$acc parallel loop default(present)
    DO i=2,imax
        um_n_tau(i,jmax) = um_n_tau(i,jmax-1)
    ENDDO
    !------ contorno esquerdo --------
    !$acc parallel loop default(present)
    DO j=1,jmax
        um_n_tau(2,j) = 0.d0 !symmetry
        um_n_tau(1,j) = 0.d0 !symmetry
    ENDDO
    !$acc parallel loop default(present)
    DO j=1,jmax+1
        vm_tau(1,j) = vm_tau(2,j) !darcy
    ENDDO
    !------ contorno direito --------
    !$acc parallel loop default(present)
    DO j=1,jmax
        um_n_tau(imax  ,j) = 0.d0
        um_n_tau(imax+1,j) = 0.d0
    ENDDO
    !$acc parallel loop default(present)
    DO j=2,jmax
        vm_tau(imax,j) =  vm_tau(imax-1,j) !Darcy
    ENDDO
    !------------------------------------
return
end subroutine solve_U


!--- ResV ---
subroutine RESV(um,vm,p,T)
    use comum
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,1:jmax) :: P,T

    !$acc parallel loop collapse(2) default(present)
    do j=3,jmax-1
      do i=3,imax-2
        fn(i,j) = 0.5d0 * ( vm(i  ,j  ) + vm(i  ,j+1) ) * areav_n(i) / epsilon1(i,j)
        fs(i,j) = 0.5d0 * ( vm(i  ,j  ) + vm(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)
        fe(i,j) = 0.5d0 * ( um(i+1,j  ) + um(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)
        fw(i,j) = 0.5d0 * ( um(i  ,j  ) + um(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)

        df(i,j) = fe(i,j) - fw(i,j) + fn(i,j) - fs(i,j)

        Dn(i,j) = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        Ds(i,j) = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        De(i,j) = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        Dw(i,j) = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))

        !quick
        if(fw(i,j).GT.0.0d0) then
                afw(i,j) = 1.0d0
        elseif(fw(i,j).LT.0.0d0) then
                afw(i,j) = 0.0d0
        endif

        if(fe(i,j).GT.0.0d0) then
                afe(i,j) = 1.0d0
        elseif(fe(i,j).LT.0.0d0) then
                afe(i,j) = 0.0d0
        endif

        if(fn(i,j).GT.0.0d0) then
                afn(i,j) = 1.0d0
        elseif(fn(i,j).LT.0.0d0) then
                afn(i,j) = 0.0d0
        endif

        if(fs(i,j).GT.0.0d0) then
                afs(i,j) = 1.0d0
        elseif(fs(i,j).LT.0.0d0) then
                afs(i,j) = 0.0d0
        endif

        aw(i,j) = Dw(i,j) + 0.75d0  * afw(i,j) * fw(i,j) &
                + 0.125d0 * afe(i,j) * fe(i,j) &
                + 0.375d0 * ( 1.0d0 - afw(i,j) ) * fw(i,j)

        ae(i,j) = De(i,j) - 0.375d0 * afe(i,j)   * fe(i,j) &
                - 0.75d0  * ( 1.0d0 - afe(i,j) ) * fe(i,j) &
                - 0.125d0 * ( 1.0d0 - afw(i,j) ) * fw(i,j)

        as(i,j) = Ds(i,j) + 0.75d0  * afs(i,j) * fs(i,j) &
                + 0.125d0 * afn(i,j) * fn(i,j) &
                + 0.375d0 * ( 1.0d0 - afs(i,j) ) * fs(i,j)

        an(i,j) = Dn(i,j) - 0.375d0 * afn(i,j)   * fn(i,j) &
                - 0.75d0  * ( 1.0d0 - afn(i,j) ) * fn(i,j) &
                - 0.125d0 * ( 1.0d0 - afs(i,j) ) * fs(i,j)

        aww(i,j) = -0.125d0 *           afw(i,j)   * fw(i,j)
        aee(i,j) =  0.125d0 * ( 1.0d0 - afe(i,j) ) * fe(i,j)
        ass(i,j) = -0.125d0 *           afs(i,j)   * fs(i,j)
        ann(i,j) =  0.125d0 * ( 1.0d0 - afn(i,j) ) * fn(i,j)

        ap(i,j) = aw(i,j) + ae(i,j) + as(i,j) + an(i,j) + aww(i,j) + aee(i,j) + ass(i,j) + ann(i,j) + df(i,j)
        !end Quick

        v_W(i,j)  = vm(i-1,j  )
        v_WW(i,j) = vm(i-2,j  )
        v_E(i,j)  = vm(i+1,j  )
        v_EE(i,j) = vm(i+2,j  )
        v_S(i,j)  = vm(i  ,j-1)
        v_SS(i,j) = vm(i  ,j-2)
        v_N(i,j)  = vm(i  ,j+1)
        v_NN(i,j) = vm(i  ,j+2)         
        v_P(i,j)  = vm(i  ,j  )
        u_P(i,j)  = um(i   ,j  )

        dvdydy(i,j) =  areav_n(i) * ( v_N(i,j) - v_P(i,j) ) / (ym(j+1)-ym(j  )) &
                -areav_s(i) * ( v_P(i,j) - v_S(i,j) ) / (ym(j  )-ym(j-1)) 
        dydudx(i,j) =  areav_n(i) * (um(i+1,j  )-um(i,j  )) / (xm(i+1)-xm(i)) &
                -areav_s(i) * (um(i+1,j-1)-um(i,j-1)) / (xm(i+1)-xm(i))

        artDivV(i,j) = - b_art *(dydudx(i,j) + dvdydy(i,j))            

        !bulk artificial viscosity term from Ramshaw(1990)
        q_art(i,j) = epsilon1(i,j)* ( p(i,j)-p(i,j-1) )/ (y(j)-y(j-1)) + artDivV(i,j)

        RV(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * (- ap(i,j) * v_P(i,j) &
                + aww(i,j) * v_WW(i,j) + aw(i,j) * v_W(i,j) &
                + aee(i,j) * v_EE(i,j) + ae(i,j) * v_E(i,j) &
                + ass(i,j) * v_SS(i,j) + as(i,j) * v_S(i,j) &
                + ann(i,j) * v_NN(i,j) + an(i,j) * v_N(i,j)) &
                - q_art(i,j) + InvFr2 * ( 1.d0 - 1.d0 / ( (T(i,j)+T(i,j-1)) * 0.5d0 ) ) &
                - epsilon1(i,j)*(   v_P(i,j)/(Re*Darcy_number) + &
                Cf/((epsilon1(i,j)*Darcy_number)**0.5d0)* v_P(i,j) * &
                ((u_P(i,j)**2.d0 +v_P(i,j)**2.d0)**0.5d0))*liga_poros(i,j) 
      enddo
    enddo

    CALL upwind_Vi(um,vm,p,T,2)
    CALL upwind_Vi(um,vm,p,T,jmax)
        
    CALL upwind_Vj(um,vm,p,T,2)
    CALL upwind_Vj(um,vm,p,T,imax-1)

return
end subroutine RESV


!--- upwind_V ---
subroutine upwind_Vi(um, vm, p, T, j)
    use comum
    implicit none
    integer :: i, j, itc
    REAL(8), DIMENSION(1:imax+1, 1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  , 1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax  , 1:jmax  ) :: P,T

    !$acc parallel loop default(present)
    do i=2, imax-1
        fn(i,j) = 0.5d0 * ( vm(i  ,j  ) + vm(i  ,j+1) ) * areav_n(i) / epsilon1(i,j) 
        fs(i,j) = 0.5d0 * ( vm(i  ,j  ) + vm(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)  
        fe(i,j) = 0.5d0 * ( um(i+1,j  ) + um(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)  
        fw(i,j) = 0.5d0 * ( um(i  ,j  ) + um(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)  

        df(i,j) = fe(i,j) - fw(i,j) + fn(i,j) - fs(i,j)

        Dn(i,j) = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        Ds(i,j) = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        De(i,j) = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        Dw(i,j) = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))

        !upwind
        aw(i,j) = Dw(i,j) + MAX(fw(i,j) , 0.0d0)
        as(i,j) = Ds(i,j) + MAX(fs(i,j) , 0.0d0)
        ae(i,j) = De(i,j) + MAX(0.0d0 , -fe(i,j))
        an(i,j) = Dn(i,j) + MAX(0.0d0 , -fn(i,j))

        ap(i,j) = aw(i,j) + ae(i,j) + as(i,j) + an(i,j) + df(i,j)

        v_W(i,j) = vm(i-1,j  )
        v_E(i,j) = vm(i+1,j  )
        v_S(i,j) = vm(i  ,j-1)
        v_N(i,j) = vm(i  ,j+1)
        v_P(i,j) = vm(i  ,j  )
        u_P(i,j) = um(i  ,j  )

        dvdydy(i,j) = areav_n(i) * ( v_N(i,j) - v_P(i,j) ) / (ym(j+1)-ym(j  )) &
                     -areav_s(i) * ( v_P(i,j) - v_S(i,j) ) / (ym(j  )-ym(j-1)) 
            
        dydudx(i,j) = areav_n(i) * (um(i+1,j  )-um(i,j  )) / (xm(i+1)-xm(i)) &
                     -areav_s(i) * (um(i+1,j-1)-um(i,j-1)) / (xm(i+1)-xm(i))
            
        !bulk artificial viscosity term from Ramshaw(1990)
        q_art(i,j) = epsilon1(i,j) * ( p(i,j)-p(i,j-1) ) / (y(j)-y(j-1)) - b_art * (dydudx(i,j) + dvdydy(i,j))
        RV(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * (- ap(i,j) * v_P(i,j) &
                + aw(i,j) * v_W(i,j) + ae(i,j) * v_E(i,j)  &
                + as(i,j) * v_S(i,j) + an(i,j) * v_N(i,j))  &
                - q_art(i,j) + InvFr2 * ( 1.d0 - 1.d0 / ( (T(i,j)+T(i,j-1)) * 0.5d0 ) ) &
                - epsilon1(i,j) * ( v_P(i,j)/(Re*Darcy_number) + &
                Cf/((epsilon1(i,j)*Darcy_number)**0.5d0)* v_P(i,j) * &
                ((u_p(i,j)**2.d0 + v_P(i,j)**2.d0)**0.5d0))*liga_poros(i,j) 
    enddo
return
end subroutine upwind_Vi


!--- upwind_V ---
subroutine upwind_Vj(um, vm, p, T, i)
    use comum
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1, 1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  , 1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax  , 1:jmax  ) :: P,T
    
    !$acc parallel loop default(present)
    do j=2, jmax
        fn(i,j) = 0.5d0 * ( vm(i  ,j  ) + vm(i  ,j+1) ) * areav_n(i) / epsilon1(i,j) 
        fs(i,j) = 0.5d0 * ( vm(i  ,j  ) + vm(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)  
        fe(i,j) = 0.5d0 * ( um(i+1,j  ) + um(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)  
        fw(i,j) = 0.5d0 * ( um(i  ,j  ) + um(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)  

        df(i,j) = fe(i,j) - fw(i,j) + fn(i,j) - fs(i,j)

        Dn(i,j) = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        Ds(i,j) = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        De(i,j) = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1) -x(i  ))
        Dw(i,j) = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  ) -x(i-1))

        !upwind
        aw(i,j) = Dw(i,j) + MAX(fw(i,j) , 0.0d0)
        as(i,j) = Ds(i,j) + MAX(fs(i,j) , 0.0d0)

        ae(i,j) = De(i,j) + MAX(0.0d0 , -fe(i,j))
        an(i,j) = Dn(i,j) + MAX(0.0d0 , -fn(i,j))

        ap(i,j) = aw(i,j) + ae(i,j) + as(i,j) + an(i,j) + df(i,j)

        v_W(i,j) = vm(i-1,j  )
        v_E(i,j) = vm(i+1,j  )
        v_S(i,j) = vm(i  ,j-1)
        v_N(i,j) = vm(i  ,j+1)
        v_P(i,j) = vm(i  ,j  )
        u_P(i,j) = um(i  ,j  )

        dvdydy(i,j) = areav_n(i) * ( v_N(i,j) - v_P(i,j) ) / (ym(j+1)-ym(j  )) &
                -areav_s(i) * ( v_P(i,j) - v_S(i,j) ) / (ym(j  )-ym(j-1)) 
            
        dydudx(i,j) = areav_n(i) * (um(i+1,j  )-um(i,j  )) / (xm(i+1)-xm(i)) &
                -areav_s(i) * (um(i+1,j-1)-um(i,j-1)) / (xm(i+1)-xm(i))
            
        !bulk artificial viscosity term from Ramshaw(1990)
        q_art(i,j) = epsilon1(i,j) * ( p(i,j)-p(i,j-1) ) / (y(j)-y(j-1)) - b_art * (dydudx(i,j) + dvdydy(i,j))
        RV(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * (- ap(i,j) * v_P(i,j) &
                + aw(i,j) * v_W(i,j) + ae(i,j) * v_E(i,j)  &
                + as(i,j) * v_S(i,j) + an(i,j) * v_N(i,j))  &
                - q_art(i,j) + InvFr2 * ( 1.d0 - 1.d0 / ( (T(i,j)+T(i,j-1)) * 0.5d0 ) ) &
                - epsilon1(i,j) * ( v_P(i,j)/(Re*Darcy_number) + &
                Cf/((epsilon1(i,j)*Darcy_number)**0.5d0)* v_P(i,j) * &
                ((u_p(i,j)**2.d0 + v_P(i,j)**2.d0)**0.5d0))*liga_poros(i,j) 
    enddo
return
end subroutine upwind_Vj


!--- solve_V ---
subroutine solve_V(um,vm,um_tau,vm_tau,vm_n_tau,p,T)
    use comum
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, vm_n_tau
    REAL(8), DIMENSION(1:imax,1:jmax) :: P,T

    CALL RESV(um_tau,vm_tau,p,T)

    !$acc parallel loop collapse(2) default(present)
    do j=3,jmax-1
        do i=2,imax-1
            res_v(i,j) = ( (vm(i,j)-vm_tau(i,j)) + RV(i,j)*dt) * dtau
            vi(i,j) = ( vm_tau(i,j) + res_v(i,j) )       
        enddo
    enddo

    !call bcUV(um_tau,vi)
    !------ contorno inferior -----------
    !$acc parallel loop default(present)
    DO i=1,imax
        vi(i,1) = 0.5d0 !v_i
        vi(i,2) = 0.5d0 !v_i !vm(i,1)     
    ENDDO
    !$acc parallel loop default(present)
    DO i=2,imax
        um_tau(i,1) = 0.d0    !!!DARCY
    ENDDO
    !------ contorno superior ----------
    !$acc parallel loop default(present)
    DO i=2,imax
        vi(i,jmax  ) = 2.d0*vi(i,jmax-2) - vi(i,jmax-1)
        vi(i,jmax+1) = 2.d0*vi(i,jmax-1) - vi(i,jmax)
    ENDDO
    !$acc parallel loop default(present)
    DO i=2,imax
        um_tau(i,jmax) = um_tau(i,jmax-1)
    ENDDO
    !------ contorno esquerdo --------
    !$acc parallel loop default(present)
    DO j=1,jmax
        um_tau(2,j) = 0.d0 !symmetry
        um_tau(1,j) = 0.d0 !symmetry
    ENDDO
    !$acc parallel loop default(present)
    DO j=1,jmax+1
        vi(1,j) = vi(2,j) !darcy
    ENDDO
    !------ contorno direito --------
    !$acc parallel loop default(present)
    DO j=1,jmax
        um_tau(imax  ,j) = 0.d0
        um_tau(imax+1,j) = 0.d0
    ENDDO
    !$acc parallel loop default(present)
    DO j=2,jmax
        vi(imax,j) =  vi(imax-1,j) !Darcy
    ENDDO
    !------------------------------------

    CALL RESV(um_tau,vi,p,T)

    !$acc parallel loop collapse(2) default(present)
    do j=3,jmax-1
        do i=2,imax-1
            res_v(i,j) =( (vm(i,j)-vm_tau(i,j)) +  RV(i,j)*dt) * dtau
            vi(i,j) =( 0.75d0 * vm_tau(i,j) + 0.25d0 * ( vi(i,j) + res_v(i,j)) )
        enddo
    enddo

    !call bcUV(um_tau,vi)
    !------ contorno inferior -----------
    !$acc parallel loop default(present)
    DO i=1,imax
        vi(i,1) = 0.5d0 !v_i
        vi(i,2) = 0.5d0 !v_i !vm(i,1)     
    ENDDO
    !$acc parallel loop default(present)
    DO i=2,imax
        um_tau(i,1) = 0.d0    !!!DARCY
    ENDDO
    !------ contorno superior ----------
    !$acc parallel loop default(present)
    DO i=2,imax
        vi(i,jmax  ) = 2.d0*vi(i,jmax-2) - vi(i,jmax-1)
        vi(i,jmax+1) = 2.d0*vi(i,jmax-1) - vi(i,jmax)
    ENDDO
    !$acc parallel loop default(present)
    DO i=2,imax
        um_tau(i,jmax) = um_tau(i,jmax-1)
    ENDDO
    !------ contorno esquerdo --------
    !$acc parallel loop default(present)
    DO j=1,jmax
        um_tau(2,j) = 0.d0 !symmetry
        um_tau(1,j) = 0.d0 !symmetry
    ENDDO
    !$acc parallel loop default(present)
    DO j=1,jmax+1
        vi(1,j) = vi(2,j) !darcy
    ENDDO
    !------ contorno direito --------
    !$acc parallel loop default(present)
    DO j=1,jmax
        um_tau(imax  ,j) = 0.d0
        um_tau(imax+1,j) = 0.d0
    ENDDO
    !$acc parallel loop default(present)
    DO j=2,jmax
        vi(imax,j) =  vi(imax-1,j) !Darcy
    ENDDO
    !------------------------------------

    CALL RESV(um_tau,vi,p,T)

    !$acc parallel loop collapse(2) default(present)
    do j=3,jmax-1
        do i=2,imax-1
            res_v(i,j) =( (vm(i,j)-vm_tau(i,j)) +  RV(i,j)*dt) * dtau
            vm_n_tau(i,j) = 1.0d0 / 3.0d0 * vm_tau(i,j) + 2.0d0 / 3.0d0 * ( vi(i,j) + res_v(i,j)) 
        enddo
    enddo
   
    !call bcUV(um_tau,vm_n_tau)
    !------ contorno inferior -----------
    !$acc parallel loop default(present)
    DO i=1,imax
        vm_n_tau(i,1) = 0.5d0 !v_i
        vm_n_tau(i,2) = 0.5d0 !v_i !vm(i,1)     
    ENDDO
    !$acc parallel loop default(present)
    DO i=2,imax
        um_tau(i,1) = 0.d0    !!!DARCY
    ENDDO
    !------ contorno superior ----------
    !$acc parallel loop default(present)
    DO i=2,imax
        vm_n_tau(i,jmax  ) = 2.d0*vm_n_tau(i,jmax-2) - vm_n_tau(i,jmax-1)
        vm_n_tau(i,jmax+1) = 2.d0*vm_n_tau(i,jmax-1) - vm_n_tau(i,jmax)
    ENDDO
    !$acc parallel loop default(present)
    DO i=2,imax
        um_tau(i,jmax) = um_tau(i,jmax-1)
    ENDDO
    !------ contorno esquerdo --------
    !$acc parallel loop default(present)
    DO j=1,jmax
        um_tau(2,j) = 0.d0 !symmetry
        um_tau(1,j) = 0.d0 !symmetry
    ENDDO
    !$acc parallel loop default(present)
    DO j=1,jmax+1
        vm_n_tau(1,j) = vm_n_tau(2,j) !darcy
    ENDDO
    !------ contorno direito --------
    !$acc parallel loop default(present)
    DO j=1,jmax
        um_tau(imax  ,j) = 0.d0
        um_tau(imax+1,j) = 0.d0
    ENDDO
    !$acc parallel loop default(present)
    DO j=2,jmax
        vm_n_tau(imax,j) =  vm_n_tau(imax-1,j) !Darcy
    ENDDO
    !------------------------------------
return
end subroutine solve_V


!--- Solve Continuity Equation ---
subroutine solve_P(um_n_tau,vm_n_tau,p)
    use comum
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_n_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_n_tau
    REAL(8), DIMENSION(1:imax,1:jmax) :: P

    ! RALSTON'S METHOD (Second Order Runge-Kutta) mass conservation
    !call res_P()
    !$acc parallel loop collapse(2) default(present)
    do j=2,jmax-1 
        do i=2,imax-1
            dudx(i,j) = um_n_tau(i+1,j) * areau_e(j) - um_n_tau(i,j) * areau_w(j)
            dvdy(i,j) = vm_n_tau(i,j+1) * areav_n(i) - vm_n_tau(i,j) * areav_s(i)
            RP(i,j) = - ( dudx(i,j) + dvdy(i,j) ) 
            pi(i,j) = p(i,j) + dtau * RP(i,j) * beta
        enddo
    enddo

    !CALL bcP(pi)
    !--------contorno inferior e superior --------
    !$acc parallel loop default(present)
    do i=1, imax
        pi(i,1  ) = pi(i,2)
        pi(i,jmax) = pi(i,jmax-1) + 1.d0*(pi(i,jmax-1)-pi(i,jmax-2))
    enddo
    !-------contorno esquerdo e direito ------- 
    !$acc parallel loop default(present)
    do j=1, jmax
        pi(1  ,j) =  pi(2,j)
        pi(imax,j) = pi(imax-1,j)
    enddo

    !call res_P()
    !$acc parallel loop collapse(2) default(present)
    do j=2,jmax-1
        do i=2,imax-1
            !dudx(i,j) = um_n_tau(i+1,j) * areau_e(j) - um_n_tau(i,j) * areau_w(j)
            !dvdy(i,j) = vm_n_tau(i,j+1) * areav_n(i) - vm_n_tau(i,j) * areav_s(i)
            !RP(i,j) = - ( dudx(i,j) + dvdy(i,j) )
            pi(i,j) = 0.75d0 * p(i,j) + 0.25d0 * (pi(i,j) + dtau * beta * RP(i,j)) 
        enddo
    enddo

    !CALL bcP(pi)
    !--------contorno inferior e superior --------
    !$acc parallel loop default(present)
    do i=1, imax
        pi(i,1  ) = pi(i,2)
        pi(i,jmax) = pi(i,jmax-1) + 1.d0*(pi(i,jmax-1)-pi(i,jmax-2))
    enddo
    !-------contorno esquerdo e direito ------- 
    !$acc parallel loop default(present)
    do j=1, jmax
        pi(1  ,j) =  pi(2,j)
        pi(imax,j) = pi(imax-1,j)
    enddo

    !call res_P()
    !$acc parallel loop collapse(2) default(present)
    do j=2,jmax-1
        do i=2,imax-1
            !dudx(i,j) = um_n_tau(i+1,j) * areau_e(j) - um_n_tau(i,j) * areau_w(j)
            !dvdy(i,j) = vm_n_tau(i,j+1) * areav_n(i) - vm_n_tau(i,j) * areav_s(i)
            RP(i,j) = - ( dudx(i,j) + dvdy(i,j) )
            res_p(i,j) = dtau * RP(i,j) * beta
            p(i,j) = 1.0d0 / 3.0d0 * p(i,j) + 2.0d0 / 3.0d0 * (pi(i,j) + res_p(i,j))
        enddo
    enddo

    !CALL bcP(p)
    !--------contorno inferior e superior --------
    !$acc parallel loop default(present)
    do i=1, imax
        p(i,1  ) = p(i,2)
        p(i,jmax) = p(i,jmax-1) + 1.d0*(p(i,jmax-1)-p(i,jmax-2))
    enddo
    !-------contorno esquerdo e direito ------- 
    !$acc parallel loop default(present)
    do j=1, jmax
        p(1  ,j) =  p(2,j)
        p(imax,j) = p(imax-1,j)
    enddo
return
end subroutine solve_P


!--- Solve Mixture Fraction ---
subroutine solve_Z(um_n_tau,vm_n_tau,Z,Z_n_tau)
    use comum
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_n_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_n_tau
    REAL(8), DIMENSION(1:imax,  1:jmax)   :: Z, Z_n_tau

    ! RALSTON'S METHOD (Second Order Runge-Kutta)
    CALL RESZ(um_n_tau,vm_n_tau,Z_n_tau)

    !$acc parallel loop collapse(2) default(present)
    do j=2,jmax-1
      do i=2,imax-1    
        res_Z(i,j) = ( (Z(i,j)-Z_n_tau(i,j)) + RZ(i,j)*dt) * dtau 
        Zi(i,j) = Z_n_tau(i,j) + res_Z(i,j)             
      enddo
    enddo

    !CALL bcZ(Zi)
    !--------contorno inferior e superior --------
    !$acc parallel loop default(present)
    do i=1, imax
        Zi(i,1  )  = Tinf !Zn(i,2)
        Zi(i,jmax) = Zi(i,jmax-1) + 1.d0*(Zi(i,jmax-1)-Zi(i,jmax-2))
    enddo
    !-------contorno esquerdo e direito -------
    !$acc parallel loop default(present)
    do j=1,jmax
        Zi(1  ,j)  = Zi(2,j)
        Zi(imax,j) = Zi(imax-1,j)
    enddo

    CALL RESZ(um_n_tau,vm_n_tau,Zi)

    !$acc parallel loop collapse(2) default(present)
    do j=2,jmax-1
      do i=2,imax-1   
        res_Z(i,j) = ( (Z(i,j)-Z_n_tau(i,j)) + RZ(i,j)*dt) * dtau
        Zi(i,j) = 0.75d0 * Z_n_tau(i,j) + 0.25d0 * (Zi(i,j) + res_Z(i,j))
      enddo
    enddo

    !CALL bcZ(Zi)
    !--------contorno inferior e superior --------
    !$acc parallel loop default(present)
    do i=1, imax
        Zi(i,1  )  = Tinf !Zn(i,2)
        Zi(i,jmax) = Zi(i,jmax-1) + 1.d0*(Zi(i,jmax-1)-Zi(i,jmax-2))
    enddo
    !-------contorno esquerdo e direito -------
    !$acc parallel loop default(present)
    do j=1,jmax
        Zi(1  ,j)  = Zi(2,j)
        Zi(imax,j) = Zi(imax-1,j)
    enddo

    CALL RESZ(um_n_tau,vm_n_tau,Zi)

    !$acc parallel loop collapse(2) default(present)
    do j=2,jmax-1
      do i=2,imax-1
        res_Z(i,j) = ((Z(i,j)-Z_n_tau(i,j)) + RZ(i,j)*dt) * dtau 
        Z_n_tau(i,j) = 1.0d0 / 3.0d0 * Z_n_tau(i,j) + 2.0d0 / 3.0d0 * (Zi(i,j) + res_Z(i,j))
      enddo
    enddo

    !CALL bcZ(Z_n_tau)
    !--------contorno inferior e superior --------
    !$acc parallel loop default(present)
    do i=1, imax
        Z_n_tau(i,1  )  = Tinf !Zn(i,2)
        Z_n_tau(i,jmax) = Z_n_tau(i,jmax-1) + 1.d0*(Z_n_tau(i,jmax-1)-Z_n_tau(i,jmax-2))
    enddo
    !-------contorno esquerdo e direito -------
    !$acc parallel loop default(present)
    do j=1,jmax
        Z_n_tau(1  ,j)  = Z_n_tau(2,j)
        Z_n_tau(imax,j) = Z_n_tau(imax-1,j)
    enddo
return
end subroutine solve_Z


!--- ResZ ---
subroutine RESZ(um_n_tau,vm_n_tau,Z)
    use comum
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1, 1:jmax  ) :: um_n_tau
    REAL(8), DIMENSION(1:imax  , 1:jmax+1) :: vm_n_tau
    REAL(8), DIMENSION(1:imax  , 1:jmax  ) :: Z !, RZ

    !$acc parallel loop collapse(2) default(present)
    do j=2,jmax-1
      do i=2,imax-1

        dZudx(i,j) = 0.5d0 * (Z(i+1,j)+Z(i,j)) * um_n_tau(i+1,j) * areau_e(j) &
                - 0.5d0 * (Z(i-1,j)+Z(i,j)) * um_n_tau(i  ,j) * areau_w(j)
        dZvdy(i,j) = 0.5d0 * (Z(i,j+1)+Z(i,j)) * vm_n_tau(i,j+1) * areav_n(i) &
                - 0.5d0 * (Z(i,j-1)+Z(i,j)) * vm_n_tau(i,j  ) * areav_s(i)

        De(i,j) = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i+1)-x(i  ))  
        Dw(i,j) = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i  )-x(i-1))  
        Dn(i,j) = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j+1)-y(j  ))  
        Ds(i,j) = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j  )-y(j-1))  
        Dp(i,j) = De(i,j) + Dw(i,j) + Dn(i,j) + Ds(i,j) 

        RZ(i,j) = 1.d0 / (xm(i+1)-xm(i)) / (ym(j+1)-ym(j)) * &
                (-Dp(i,j)*Z(i,j) + De(i,j)*Z(i+1,j) + Dw(i,j)*Z(i-1,j) + &
                Dn(i,j)*Z(i,j+1) + Ds(i,j)*Z(i,j-1) - &
                (1.d0-liga_poros(i,j))*(dZudx(i,j) + dZvdy(i,j)) ) / &
                (liga_poros(i,j)*(epsilon1(i,j)-1.d0)+1.d0)
       enddo
     enddo

return
end subroutine RESZ


!--- solve_C - concentration ---
subroutine solve_C(um_n_tau,vm_n_tau,C,C_n_tau)
    use comum
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_n_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_n_tau
    REAL(8), DIMENSION(1:imax,  1:jmax)   :: C, C_n_tau

    ! RALSTON'S METHOD (Second Order Runge-Kutta)
    CALL RESC(um_n_tau,vm_n_tau,C_n_tau)

    !$acc parallel loop collapse(2) default(present)
    do j=2,jmax-1
      do i=2,imax-1
        res_C(i,j) = ((C(i,j)-C_n_tau(i,j)) + RC(i,j)*dt) * dtau 
        Ci(i,j) = C_n_tau(i,j) + res_C(i,j)
      enddo
    enddo
    
    !CALL bcC(Ci)
    !--------contorno inferior e superior --------
    !$acc parallel loop default(present)
    do i=1,imax
        Ci(i,1  ) = 0.d0 !Zn(i,2)
        Ci(i,jmax) = Ci(i,jmax-1) + 1.d0*(Ci(i,jmax-1)-Ci(i,jmax-2))
    enddo
    !-------contorno esquerdo e direito -------
    !$acc parallel loop default(present)
    do j=1,jmax
        Ci(1  ,j) = Ci(2,j)
        Ci(imax,j) = Ci(imax-1,j)
    enddo

    CALL RESC(um_n_tau,vm_n_tau,Ci)
    
    !$acc parallel loop collapse(2) default(present)
    do j=2,jmax-1
      do i=2,imax-1
        res_C(i,j) = ((C(i,j)-C_n_tau(i,j)) + RC(i,j)*dt) * dtau 
        Ci(i,j) = 0.75d0 * C_n_tau(i,j) + 0.25d0 * (Ci(i,j) + res_C(i,j))   
      enddo        
    enddo

    !CALL bcC(Ci)
    !--------contorno inferior e superior --------
    !$acc parallel loop default(present)
    do i=1,imax
        Ci(i,1  ) = 0.d0 !Zn(i,2)
        Ci(i,jmax) = Ci(i,jmax-1) + 1.d0*(Ci(i,jmax-1)-Ci(i,jmax-2))
    enddo
    !-------contorno esquerdo e direito -------
    !$acc parallel loop default(present)
    do j=1,jmax
        Ci(1  ,j) = Ci(2,j)
        Ci(imax,j) = Ci(imax-1,j)
    enddo

    CALL RESC(um_n_tau,vm_n_tau,Ci)

    !$acc parallel loop collapse(2) default(present)
    do j=2,jmax-1
      do i=2,imax-1
        res_C(i,j) = ((C(i,j)-C_n_tau(i,j)) + RC(i,j)*dt) * dtau 
        C_n_tau(i,j) = 1.0d0 / 3.0d0 * C_n_tau(i,j) + 2.0d0 / 3.0d0 * (Ci(i,j) + res_C(i,j))
      enddo
    enddo

    !CALL bcC(C_n_tau)
    !--------contorno inferior e superior --------
    !$acc parallel loop default(present)
    do i=1,imax
        C_n_tau(i,1  ) = 0.d0 !Zn(i,2)
        C_n_tau(i,jmax) = C_n_tau(i,jmax-1) + 1.d0*(C_n_tau(i,jmax-1)-C_n_tau(i,jmax-2))
    enddo
    !-------contorno esquerdo e direito -------
    !$acc parallel loop default(present)
    do j=1,jmax
        C_n_tau(1  ,j) = C_n_tau(2,j)
        C_n_tau(imax,j) = C_n_tau(imax-1,j)
    enddo
return
end subroutine solve_C


!--- ResC ---
subroutine RESC(um_n_tau, vm_n_tau, C)
    use comum
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1, 1:jmax  ) :: um_n_tau
    REAL(8), DIMENSION(1:imax  , 1:jmax+1) :: vm_n_tau
    REAL(8), DIMENSION(1:imax,   1:jmax)   :: C

    !$acc parallel loop collapse(2) default(present)
    do j=2,jmax-1
      do i=2,imax-1
        dCudx(i,j) = 0.5d0 * (C(i+1,j)+C(i,j)) * um_n_tau(i+1,j) * areau_e(j) &
                   - 0.5d0 * (C(i-1,j)+C(i,j)) * um_n_tau(i  ,j) * areau_w(j)
        dCvdy(i,j) = 0.5d0 * (C(i,j+1)+C(i,j)) * vm_n_tau(i,j+1) * areav_n(i) &
                   - 0.5d0 * (C(i,j-1)+C(i,j)) * vm_n_tau(i,j  ) * areav_s(i)

        De(i,j) = (ym(j+1)-ym(j)) * (1.d0/Re/Sc) / (x(i+1)-x(i  ))
        Dw(i,j) = (ym(j+1)-ym(j)) * (1.d0/Re/Sc) / (x(i  )-x(i-1))
        Dn(i,j) = (xm(i+1)-xm(i)) * (1.d0/Re/Sc) / (y(j+1)-y(j  ))
        Ds(i,j) = (xm(i+1)-xm(i)) * (1.d0/Re/Sc) / (y(j  )-y(j-1))
        Dp(i,j) = De(i,j) + Dw(i,j) + Dn(i,j) + Ds(i,j)

        RC(i,j) = 1.d0 / (xm(i+1)-xm(i)) / (ym(j+1)-ym(j)) *&
                (-Dp(i,j)*C(i,j) + De(i,j)*C(i+1,j) + Dw(i,j)*C(i-1,j) + &
                  Dn(i,j)*C(i,j+1) + Ds(i,j)*C(i,j-1) - &
                (dCudx(i,j) + dCvdy(i,j)) / &
                (liga_poros(i,j)*(epsilon1(i,j)-1.d0)+1.d0) )
       enddo
     enddo

return
end subroutine RESC
