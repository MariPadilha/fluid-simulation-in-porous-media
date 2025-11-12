subroutine transient(u,v,p,T,C,tr)
    use comum
    implicit none
    integer :: i, j, tr
    REAL(8), DIMENSION(1:imax,1:jmax) :: u, v, P, C, T
    
    character(len=100) :: filename2

    write(filename2,*) (tr/n_tr)+100
    filename2 = adjustl(filename2)   

    
    open(unit=550, file='transient/data/u'//trim(filename2)//'.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(u(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)

    open(unit=550, file='transient/data/v'//trim(filename2)//'.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(v(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)
    
    open(unit=550, file='transient/data/P'//trim(filename2)//'.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(P(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)
    
    open(unit=550, file='transient/data/T'//trim(filename2)//'.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(T(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)
    
    open(unit=550, file='transient/data/C'//trim(filename2)//'.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(C(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)

    open (550,file='transient/data/time'//trim(filename2)//'.dat')

                write (550,*) time

    close(550)
    
    open (550,file='transient/data/grid.dat')


         do j=1,jmax
                do i = 1,imax
            
                write (550,*) X(i) , Y(j)

                enddo
         enddo    

    close(550)

!    i = 1
!    do j=int(y_down/dy),jmax
!                
!        if (Zc(i,j) .gt. 1.d0) then
!
!        yf = ((1.d0-Zc(i,j-1))*( Y(j)- Y(j-1) )) / ( Zc(i,j) - Zc(i,j-1) ) + Y(j-1)
!        
!        go to 200        
!
!        endif    
!        
!    200 enddo    
!
!    open(unit=550,file='data/flametip.dat',status='unknown',position='append')
!        
!        if (yf .lt. y_up) then
!            write(550,*) itc,yf
!        else
!            write(550,*) itc, y_up
!        endif
!              
!        
!    close(550)
!

!########## OUTPUT TO USE IN PARAVIEW APP ################################
    OPEN(unit=1,file='data/paraview_output'//trim(filename2)//'.vtk',status='unknown')

    WRITE(1,'(a)')'# vtk DataFile Version 3.0' 
    WRITE(1,'(a)')'Droplet Combustion'   
    WRITE(1,'(a)')'ASCII'   
    WRITE(1,'(a)')'DATASET STRUCTURED_GRID'  
    WRITE(1,'(A10,A1,I3,A1,I3,2A1)')'DIMENSIONS',' ',jmax,' ',imax,' ','1'
!WRITE(1,fmt='(2i4.2)') nj,ni

!WRITE(1,'(A10,I2,I2,A1)')'DIMENSIONS', nj,ni,'1'
!WRITE(1,'(A,I8.0)')'DIMENSIONS', nj

    WRITE(1,'(A6,A1,I6,A1,A5)')'POINTS',' ', jmax*imax,' ','float'
    DO i=1,imax
        DO j=1,jmax
            WRITE(1,'(3F14.4)') x(i), y(j), 0.d0
        ENDDO
    ENDDO 
    WRITE(1,'(A10,A1,I6)')'POINT_DATA',' ', imax*jmax
    WRITE(1,'(a)')'SCALARS VEL_MAGNITUDE float'
    WRITE(1,'(a)')'LOOKUP_TABLE default'
    DO i=1,imax
        DO j=1,jmax
            WRITE(1,'(F14.4)') dsqrt((u(i,j))**2+(v(i,j))**2)
        ENDDO
    ENDDO

    WRITE(1,'(a)')'SCALARS T float'
    WRITE(1,'(a)')'LOOKUP_TABLE default'
    DO i=1,imax
        DO j=1,jmax
            WRITE(1,'(F14.4)') 	T(i,j)
        ENDDO
    ENDDO

    WRITE(1,'(a)')'SCALARS P float'
    WRITE(1,'(a)')'LOOKUP_TABLE default'
    DO i=1,imax
        DO j=1,jmax
            WRITE(1,'(F14.4)') 	P(i,j)
        ENDDO
    ENDDO

 WRITE(1,'(a)')'VECTORS U float'
    WRITE(1,'(a)')'LOOKUP_TABLE default'
    DO i=1,imax
        DO j=1,jmax
            WRITE(1,'(F14.4)') 	U(i,j)
        ENDDO
    ENDDO
    
WRITE(1,'(a)')'VECTORS V float'
    WRITE(1,'(a)')'LOOKUP_TABLE default'
    DO i=1,imax
        DO j=1,jmax
            WRITE(1,'(F14.4)') 	V(i,j)
        ENDDO
    ENDDO




!    WRITE(1,'(a)')'SCALARS Z float'
!!    WRITE(1,'(a)')'LOOKUP_TABLE default'
!    DO i=1,imax
!        DO j=1,jmax
!            WRITE(1,'(F14.4)') 	Z(i,j)
!        ENDDO
!    ENDDO

!    WRITE(1,'(a)')'SCALARS H float'
!    WRITE(1,'(a)')'LOOKUP_TABLE default'
!    DO i=1,imax
!        DO j=1,jmax
!            WRITE(1,'(F14.4)') 	H(i,j)
!        ENDDO
!    ENDDO

  !  WRITE(1,'(a)')'SCALARS Yi float'
  !  WRITE(1,'(a)')'LOOKUP_TABLE default'
   ! DO i=1,imax
    !    DO j=1,jmax
     !       WRITE(1,'(F14.4)') 	C(i,j)
      !  ENDDO
    !ENDDO

   ! WRITE(1,'(a)')'SCALARS U_V float'
    !WRITE(1,'(a)')'LOOKUP_TABLE default'
    !WRITE(1,'(a)')'VECTORS Vectors float'

  !  DO i=1,imax
   !     DO j=1,jmax		
    !        WRITE(1,'(3F14.4)') 	u(i,j), v(i,j), 0.d0
    !    ENDDO
    ! ENDDO

   ! CLOSE(1)


return
end subroutine transient


