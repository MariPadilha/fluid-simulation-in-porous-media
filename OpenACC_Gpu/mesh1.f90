SUBROUTINE mesh1
    USE comum
    IMPLICIT NONE
    INTEGER :: i, j, ii, jj
    !parametros de refinamento, max = 2
    REAL(8) :: etax, etay, sx, sy

    DO i=1,imax                              
        X(i) = ( (dfloat(i) - 1.d0) * dx_c )  
    ENDDO

    DO j=1,jmax
        Y(j) = ( (dfloat(j) - 1.d0) * dx_c ) - y_down
        
        write(*,*) y(j)
    ENDDO

!############### REFINAMENTO NA MALHA PARA X ###############
    DO i = 1,imax
            etax = ( X(i) - X(imax) ) / ( X(1) - X(imax) )
            sx = Px_grid * etax + (1.d0 - Px_grid) &
                 * (1.d0 - ( (tanh(Q_grid * (1.d0 - etax)) ) / tanh(Q_grid) ) )
            X(i) = X(imax) - sx * (X(imax) - X(1) )
    ENDDO 

!############### REFINAMENTO NA MALHA PARA Y ###############
    DO j = int(y_down/dx_c)+1,jmax
        etay = ( Y(j) - Y(jmax) ) / ( Y(int(y_down/dx_c)+1) - Y(jmax) )
        sy = Py_grid * etay + (1.d0 - Py_grid) &
             * (1.d0 - ( (tanh(Q_grid * (1.d0 - etay)) ) / tanh(Q_grid) ) )
        Y(j) = Y(jmax) - sy * (Y(jmax) - Y(int(y_down/dx_c)+1) )
    ENDDO
    

    DO j = int(y_down/dx_c)+1,1,-1
        etay = ( Y(j) + Y(jmax) ) / ( Y(int(y_down/dx_c)+1) + Y(jmax) )
        sy = Py_grid * etay + (1.d0 - Py_grid) &
             * (1.d0 - ( (tanh(Q_grid * (1.d0 - etay)) ) / tanh(Q_grid) ) )
        Y(j) = -Y(jmax) - sy * (-Y(jmax) - Y(int(y_down/dx_c)+1) )
    ENDDO

!#############################################################
 
    OPEN (1,file="data/grid.dat")

    DO j=1,jmax,2 !plotar na direção de i
        DO i=1,imax
            WRITE(1,*) x(i),y(j)
        END DO

        IF (j<jmax) THEN
            jj=j+1
        DO i=imax,1,-1
            WRITE(1,*) x(i),y(jj)
        END DO
        END IF
    END DO

    DO i=1,imax,2 !plotar na direção de j
        DO j=jmax,1,-1
            WRITE(1,*) x(i),y(j)
        END DO

        IF (i<imax) THEN
            ii=i+1
            DO j=1,jmax
            WRITE(1,*) x(ii),y(j)
            END DO
        END IF
    END DO
