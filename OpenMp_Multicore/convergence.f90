!--- check convergence ---
SUBROUTINE convergence(itc, error, residual_p, residual_u, residual_v)
    use comum
    IMPLICIT NONE
    
    integer :: i, j, itc
    real(8) :: error, residual_p, residual_u, residual_v

    itc = itc+1
    error = MAX(residual_u, residual_v, residual_p)

    !IF (MOD(itc,1000).EQ.0) THEN
    !    WRITE(*,*) itc, ' ', residual_u, residual_v, residual_p
    !ENDIF

    !OPEN(unit=549,file='data/error.dat',status='unknown',position='appEND')
    !WRITE(549,*) itc, ' ', residual_u, residual_v, residual_p
    !close(549)
   
RETURN
END SUBROUTINE convergence
