subroutine bcUV(um,vm)
    use comum
    implicit none
    integer :: i, j
    real(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    real(8), DIMENSION(1:imax  ,1:jmax+1) :: vm

    !------ contorno inferior -----------
    DO i=1,imax
        vm(i,1) = v_i
        vm(i,2) = vm(i,1)     
    ENDDO

    DO i=2,imax
        um(i,1) = 0.d0    !!!DARCY
    ENDDO

    !------ contorno superior ----------
    DO i=2,imax
        vm(i,jmax  ) =2.d0*vm(i,jmax-2) - vm(i,jmax-1)
        vm(i,jmax+1) = 2.d0*vm(i,jmax-1) - vm(i,jmax)
    ENDDO

    DO i=2,imax
        um(i,jmax) = um(i,jmax-1) 
    ENDDO
    
    !------ contorno esquerdo --------
    DO j=1,jmax
        um(2,j) = 0.d0 !symmetry
        um(1,j) = 0.d0 !symmetry
    ENDDO

    DO j=1,jmax+1
        vm(1,j) = vm(2,j) !darcy
    ENDDO

    !------ contorno direito --------
    DO j=1,jmax
        um(imax  ,j) = 0.d0
        um(imax+1,j) = 0.d0        
    ENDDO

    DO j=2,jmax
        vm(imax,j) =  vm(imax-1,j) !Darcy
    ENDDO

return
end subroutine bcUV


SUBROUTINE bcP(pn)
    USE comum
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax,1:jmax) :: pn

    !boundary condition
    
    !--------contorno inferior e superior --------
    do i=1, imax
        pn(i,1  ) = pn(i,2)
        pn(i,jmax) = pn(i,jmax-1) + 1.d0*(pn(i,jmax-1)-pn(i,jmax-2))
    enddo
    
    !-------contorno esquerdo e direito -------    
    do j=1, jmax
        pn(1  ,j) =  pn(2,j)
        pn(imax,j) = pn(imax-1,j)
    enddo

RETURN
END SUBROUTINE bcP


SUBROUTINE bcZ(Zt)
    USE comum
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax,1:jmax) :: Zt,gradZ_x,gradZ_y

    !boundary condition
    
    !--------contorno inferior e superior --------
    do i=1, imax
        Zt(i,1  )  = Tinf !Zn(i,2)
        Zt(i,jmax) = Zt(i,jmax-1) + 1.d0*(Zt(i,jmax-1)-Zt(i,jmax-2))
    enddo
    
    !-------contorno esquerdo e direito -------    
    do j=1,jmax
        Zt(1  ,j)  = Zt(2,j)
        Zt(imax,j) = Zt(imax-1,j)
    enddo

    i=1
    do j=2,jmax-1
        gradZ_x(i,j) = x(i) * ( Zt(i+1,j) - Zt(i,j  ) ) * dx(i)
        gradZ_y(i,j) = y(j) * ( Zt(i,j+1) - Zt(i,j  ) ) * dy(j)
    enddo

    j=1
    do i=1,imax-1
        gradZ_x(i,j) = x(i) * ( Zt(i+1,j) - Zt(i,j  ) ) * dx(i+1)
        gradZ_y(i,j) = y(i) * ( Zt(i,j+1) - Zt(i,j  ) ) * dy(j+1)
    enddo
        
    do j=2,jmax
        do i=2,imax
            gradZ_x(i,j) = x(i) * ( Zt(i,j) - Zt(i-1,j  ) ) * dx(i)
            gradZ_y(i,j) = y(j) * ( Zt(i,j) - Zt(i  ,j-1) ) * dy(j)
        enddo
    enddo
  
RETURN
END SUBROUTINE bcZ


SUBROUTINE bcC(Cn)
    USE comum
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax,1:jmax) :: Cn,gradC_x,gradC_y

    !boundary condition
    
    !--------contorno inferior e superior --------
    do i=1,imax
        Cn(i,1  ) = 0.d0 !Zn(i,2)
        Cn(i,jmax) = Cn(i,jmax-1) + 1.d0*(Cn(i,jmax-1)-Cn(i,jmax-2))
    enddo
    
    !-------contorno esquerdo e direito -------    
    do j=1,jmax
        Cn(1  ,j) = Cn(2,j)
        Cn(imax,j) = Cn(imax-1,j)
    enddo

    i=1
    do j=2,jmax-1
        gradC_x(i,j) = x(i) * ( Cn(i+1,j) - Cn(i,j  ) ) * dx(i)
        gradC_y(i,j) = y(j) * ( Cn(i,j+1) - Cn(i,j  ) ) * dy(j)
    enddo

    j=1
    do i=1,imax-1
        gradC_x(i,j) = x(i) * ( Cn(i+1,j) - Cn(i,j  ) ) * dx(i+1)
        gradC_y(i,j) = y(i) * ( Cn(i,j+1) - Cn(i,j  ) ) * dy(j+1)        
    enddo
        
    do j=2, jmax
        do i=2, imax
            gradC_x(i,j) = x(i) * ( Cn(i,j) - Cn(i-1,j  ) ) * dx(i)
            gradC_y(i,j) = y(j) * ( Cn(i,j) - Cn(i  ,j-1) ) * dy(j)
        enddo
    enddo
        
RETURN
END SUBROUTINE bcC
