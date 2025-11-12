subroutine bcUV(um,vm)
    use comum
    implicit none
    integer :: i, j
    real(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    real(8), DIMENSION(1:imax  ,1:jmax+1) :: vm

    !------ contorno inferior -----------
    !$acc parallel loop default(present)
    DO i=1,imax
        vm(i,1) = 0.5d0 !v_i
        vm(i,2) = 0.5d0 !v_i !vm(i,1)     
    ENDDO
    !$acc parallel loop default(present)
    DO i=2,imax
        um(i,1) = 0.d0    !!!DARCY
    ENDDO

    !------ contorno superior ----------
    !$acc parallel loop default(present)
    DO i=2,imax
        vm(i,jmax  ) = 2.d0*vm(i,jmax-2) - vm(i,jmax-1)
        vm(i,jmax+1) = 2.d0*vm(i,jmax-1) - vm(i,jmax)
    ENDDO
    !$acc parallel loop default(present)
    DO i=2,imax
        um(i,jmax) = um(i,jmax-1) 
    ENDDO
    
    !------ contorno esquerdo --------
    !$acc parallel loop default(present)
    DO j=1,jmax
        um(2,j) = 0.d0 !symmetry
        um(1,j) = 0.d0 !symmetry
    ENDDO
    !$acc parallel loop default(present)
    DO j=1,jmax+1
        vm(1,j) = vm(2,j) !darcy
    ENDDO

    !------ contorno direito --------
    !$acc parallel loop default(present)
    DO j=1,jmax
        um(imax  ,j) = 0.d0
        um(imax+1,j) = 0.d0        
    ENDDO
    !$acc parallel loop default(present)
    DO j=2,jmax
        vm(imax,j) =  vm(imax-1,j) !Darcy
    ENDDO

return
end subroutine bcUV


!SUBROUTINE bcP(pn)
!    USE comum
!    IMPLICIT NONE
!    INTEGER :: i, j
!    REAL(8), DIMENSION(1:imax,1:jmax) :: pn

!    !boundary condition
!    !$acc parallel loop default(present)    
!    !--------contorno inferior e superior --------
!    do i=1, imax
!        pn(i,1  ) = pn(i,2)
!        pn(i,jmax) = pn(i,jmax-1) + 1.d0*(pn(i,jmax-1)-pn(i,jmax-2))
!    enddo
!    !$acc parallel loop default(present)
!    !-------contorno esquerdo e direito -------    
!    do j=1, jmax
!        pn(1  ,j) =  pn(2,j)
!        pn(imax,j) = pn(imax-1,j)
!    enddo
!RETURN
!END SUBROUTINE bcP


!SUBROUTINE bcZ(Zt)
!    USE comum
!    IMPLICIT NONE
!    INTEGER :: i, j
!    REAL(8), DIMENSION(1:imax,1:jmax) :: Zt,gradZ_x,gradZ_y
!    !boundary condition   
!    !--------contorno inferior e superior --------
!    !$acc parallel loop default(present)
!    do i=1, imax
!        Zt(i,1  )  = Tinf !Zn(i,2)
!        Zt(i,jmax) = Zt(i,jmax-1) + 1.d0*(Zt(i,jmax-1)-Zt(i,jmax-2))
!    enddo
!    !-------contorno esquerdo e direito -------    
!    !$acc parallel loop default(present)
!    do j=1,jmax
!        Zt(1  ,j)  = Zt(2,j)
!        Zt(imax,j) = Zt(imax-1,j)
!    enddo

    !i=1
    ! $acc parallel loop default(present)
    !do j=2,jmax-1
    !    gradZ_x(1,j) = x(i) * ( Zt(1+1,j) - Zt(1,j  ) ) * dx(1)
    !    gradZ_y(1,j) = y(j) * ( Zt(1,j+1) - Zt(1,j  ) ) * dy(j)
    !enddo

    !j=1
    ! $acc parallel loop default(present)
    !do i=1,imax-1
    !    gradZ_x(i,1) = x(i) * ( Zt(i+1,1) - Zt(i,1  ) ) * dx(i+1)
    !    gradZ_y(i,1) = y(i) * ( Zt(i,1+1) - Zt(i,1  ) ) * dy(1+1)
    !enddo
    ! $acc parallel loop default(present)
    !do j=2,jmax
    !    do i=2,imax
    !        gradZ_x(i,1) = x(i) * ( Zt(i,1) - Zt(i-1,1  ) ) * dx(i)
    !        gradZ_y(i,1) = y(1) * ( Zt(i,1) - Zt(i  ,1-1) ) * dy(1)
    !    enddo
    !enddo
  
!RETURN
!END SUBROUTINE bcZ


!SUBROUTINE bcC(Cn)
!    USE comum
!    IMPLICIT NONE
!    INTEGER :: i, j
!    REAL(8), DIMENSION(1:imax,1:jmax) :: Cn,gradC_x,gradC_y
!    !boundary condition    
!    !--------contorno inferior e superior --------
!    !$acc parallel loop default(present)
!    do i=1,imax
!        Cn(i,1  ) = 0.d0 !Zn(i,2)
!        Cn(i,jmax) = Cn(i,jmax-1) + 1.d0*(Cn(i,jmax-1)-Cn(i,jmax-2))
!    enddo
!    !-------contorno esquerdo e direito -------    
!    !$acc parallel loop default(present)
!    do j=1,jmax
!        Cn(1  ,j) = Cn(2,j)
!        Cn(imax,j) = Cn(imax-1,j)
!    enddo

    !i=1
    ! $acc parallel loop default(present)
    !do j=2,jmax-1
    !    gradC_x(1,j) = x(1) * ( Cn(1+1,j) - Cn(1,j  ) ) * dx(1)
    !    gradC_y(1,j) = y(j) * ( Cn(1,j+1) - Cn(1,j  ) ) * dy(j)
    !enddo

    !j=1
    ! $acc parallel loop default(present)
    !do i=1,imax-1
    !    gradC_x(i,1) = x(i) * ( Cn(i+1,1) - Cn(i,1  ) ) * dx(i+1)
    !    gradC_y(i,1) = y(i) * ( Cn(i,1+1) - Cn(i,1  ) ) * dy(1+1)        
    !enddo
    ! $acc parallel loop default(present)    
    !do j=2, jmax
    !    do i=2, imax
    !        gradC_x(i,j) = x(i) * ( Cn(i,j) - Cn(i-1,j  ) ) * dx(i)
    !        gradC_y(i,j) = y(j) * ( Cn(i,j) - Cn(i  ,j-1) ) * dy(j)
    !    enddo
    !enddo
        
!RETURN
!END SUBROUTINE bcC
