SUBROUTINE output(um,vm,u,v,p,T,C,k)
    USE COMUM
    IMPLICIT NONE
    INTEGER :: i,j,k

    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,1:jmax) :: u,v,P,Z,T,H,Yi,C
    character*16 filename
    character(len=100) :: filename2 

    open (550,file='data/output_variables.dat')
    do i = 1,imax
        do j = 1,jmax
            write (550,*) X(i), Y(j), u(i,j), v(i,j), p(i,j), T(i,j), C(i,j)
        enddo
    enddo
    close(550)

    !--- RESTART/RESTART.dat ---
    open (550,file='data/restart/restartU.dat',status='unknown')
    do i=1,imax+1
        do j = 1,jmax
            write (550,*) um(i,j)
        enddo
    enddo
    close(550)

   open (550,file='data/restart/restartV.dat',status='unknown')
   do i=1,imax
       do j = 1,jmax+1
            write (550,*) vm(i,j)    
       enddo
   enddo
   close(550)

   open (550,file='data/restart/restartPTC.dat',status='unknown')
    do i=1,imax
        do j = 1,jmax
            write (550,*) p(i,j),T(i,j),C(i,j)       
        enddo
    enddo
    close(550)

RETURN
END SUBROUTINE output
