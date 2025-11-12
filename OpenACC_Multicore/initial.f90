SUBROUTINE init
    USE comum
    IMPLICIT NONE
    
    !--- iterations values ---
    OPEN(unit=550,file='input/iterations.dat',status='old')
    READ(550,nml=iterations)
    WRITE(*,nml=iterations)
    CLOSE(unit=550)
    
    !--- reference values ---
    OPEN(unit=550,file='input/reference.dat',status='old')
    READ(550,nml=ref)
    WRITE(*,nml=ref)
    CLOSE(unit=550)

END SUBROUTINE init


SUBROUTINE IC(um,vm,p,T,C) !condicoes iniciais
    USE comum
    IMPLICIT NONE

    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1, 1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  , 1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax  , 1:jmax)   :: P, T, C
    REAL(8) :: eta_e
    
    vm = 1.0d-3*v_i
    um = 0.d0
    p = 1.0d0
    T = Tinf
    C = 0.d0

    DO j = 1, jmax
        DO i = 1, imax
            IF (flag(i,j) .NE. C_F) THEN
                T(i,j) = Temp_cylinder
                C(i,j) = concentracao_inicial
            ENDIF
        ENDDO
    ENDDO 

    call system ('rm data/flametip.dat')
    call system ('rm data/error.dat')

RETURN
END SUBROUTINE IC


SUBROUTINE restart(um,vm,p,T,C)
    USE comum
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1, 1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  , 1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax  , 1:jmax)   :: P,T,C

    WRITE(*,*) 'RESTARTING PROGRAM'

    OPEN (550,file='data/restart/restartU.dat',status='old',access='sequential')
    do j = 1, jmax
        do i = 1, imax+1
            read (550, *) um(i, j)    
        enddo
    enddo
    CLOSE(550)

    OPEN (550,file='data/restart/restartV.dat',status='old',access='sequential')
    do j=1, jmax+1
        do i = 1, imax
            read (550, *) vm(i, j)
        enddo
    enddo
    CLOSE(550)

    OPEN (550,file='data/restart/restartPTC.dat',status='old',access='sequential')
    do j=1, jmax
        do i = 1, imax
            read (550,*) p(i, j), T(i, j), C(i, j)
        enddo
    enddo
    CLOSE(550)

RETURN
END SUBROUTINE restart


SUBROUTINE restart_dom(um,vm,p,T,Z)
    USE comum
    IMPLICIT NONE
    INTEGER :: i, j
    INTEGER, parameter :: rimax = 41  !em x
    INTEGER, parameter :: rjmax = 321 !em y
    REAL(8), DIMENSION(1:rimax+1, 1:rjmax  ) :: umr
    REAL(8), DIMENSION(1:rimax  , 1:rjmax+1) :: vmr
    REAL(8), DIMENSION(1:rimax  , 1:rjmax)   :: Pres, Zr, Tr, Hr, H_res
    REAL(8), DIMENSION(1:imax+1, 1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  , 1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax  , 1:jmax)   :: P,Z,T

    WRITE(*,*) 'RESTARTING PROGRAM'

    OPEN (550,file='data/restart/restartU.dat',status='old',access='sequential')
    do j = 1, rjmax
        do i = 1, rimax + 1
            read (550,*) umr(i,j)
        enddo
    enddo
    CLOSE(550)

    OPEN (550,file='data/restart/restartV.dat',status='old',access='sequential')
    do j = 1, rjmax + 1
        do i = 1, rimax
            read (550,*) vmr(i,j)
        enddo
    enddo
    CLOSE(550)

    OPEN (550,file='data/restart/restartPTZH.dat',status='old',access='sequential')
    do j = 1, rjmax
        do i = 1, rimax    
            read (550,*) Pres(i,j), Tr(i,j), Zr(i,j), H_res(i,j)
            Hr(i,j) = H_res(i,j) + ( ((S + 1.d0) * Lf * Tinf / q + 1.d0) - H_res(i,j) )
        enddo
    enddo
    CLOSE(550)

    DO j=1, jmax
        DO i = 1, imax
            IF (j .LE. rjmax) THEN
                p(i,j) = pres(i,j)
                T(i,j) = Tr(i,j)
                Z(i,j) = Zr(i,j)
                !H(i,j) = Hr(i,j)
            ELSE
                p(i,j) = p(i,j-1)
                T(i,j) = T(i,j-1)
                Z(i,j) = Z(i,j-1)
                !H(i,j) = H(i,j-1)
            ENDIF
        ENDDO
    ENDDO

    DO j = 1, jmax
        DO i = 1, imax+1
            IF (j .LE. rjmax) THEN
                um(i,j) = umr(i,j)
            ELSE
                um(i,j) = um(i,j-1)
            ENDIF
        ENDDO
    ENDDO

    DO j = 1, jmax+1
        DO i = 1, imax
            IF (j .LE. rjmax) THEN
               vm(i,j) = vmr(i,j)
            ELSE
               vm(i,j) = vm(i,j-1)
            ENDIF
        ENDDO
    ENDDO

RETURN
END SUBROUTINE restart_dom 
