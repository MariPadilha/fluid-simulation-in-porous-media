MODULE comum

    ! Declarações de NAMELIST para iterações e referências	
    NAMELIST /iterations/ itc_max, nc, n_tr, n_out, n_vort, beta, b_art, &
                          dtau_f,final_time, eps, eps_mass, start_mode
    NAMELIST /ref/ Tnu, YF_b, YO_oo, Ts, TnToo

    ! Parâmetros de iteração e controle
    INTEGER :: itc_max        !numero de iterações
        
    ! Frequência dos outputs:
    INTEGER :: nc             !erros
    INTEGER :: n_tr           !plota parte transiente
    INTEGER :: n_out          !salva resultados preliminares
    INTEGER :: n_vort         !salva dados do vortice

    REAL(8) :: beta           !parâmetro de compressibilidade        
    REAL(8) :: b_art          !coeficiente da dissipação artificial
    REAL(8) :: dtau_f         !fator de correcao para calc de dt        
    REAL(8) :: final_time     !tempo máximo de duração do tempo da simulação
    REAL(8) :: eps, eps_mass  !criterio de convergencia  
    INTEGER :: restart_mode   !tipo de start, se eh CI ou solucao anterior

    ! Passo de tempo
    REAL(8) :: dtau, dt
    REAL(8) :: time

    ! Parâmetros físicos e geométricos
    REAL(8), PARAMETER :: porosidade=0.5d0 !Lido em main, mesh, nonsymetric_mesh
    REAL(8), PARAMETER :: Darcy_number = 1.0d-2 !Lido em equations
    REAL(8), PARAMETER :: Cf = 1.75d0 / ((150.d0 * porosidade**3.d0)**0.5d0)
    REAL(8), PARAMETER :: Temp_cylinder = 1.d0, concentracao_inicial=1.0d0

    ! Parâmetros de refinamento P e Q
    !colocar P=1 para a malha uniforme
    !REAL(8), PARAMETER :: Px_grid = 1.6d0, Py_grid = 1.6d0, Q_grid = 1.8d0 
    REAL(8), PARAMETER :: Px_grid = 1.d0, Py_grid = 1.6d0, Q_grid = 1.0d0 

    ! Tamanho do domínio e da malha
    REAL(8), PARAMETER :: Lhori = 12.    !largura total
    REAL(8), PARAMETER :: y_up = 15.     !altura em y+
    REAL(8), PARAMETER :: y_down = 5.    !altura em y-
    REAL(8), PARAMETER :: Hvert = y_up + y_down !altura total
    INTEGER, PARAMETER :: imax = 51*3 !204 !153 !10 !51 !*10 !numero de pontos da malha em x
    REAL(8), PARAMETER :: dx_c = Lhori / (imax-1) !dita o tamanho de dy e dx
    INTEGER, PARAMETER :: jmax = (Hvert / dx_c) + 1!numero de pontos da malha em y
    !INTEGER, PARAMETER :: imax = int_points + 1 !(Lhori/dx_c)+1!numero de pontos da malha em x

    ! Main data structures for control of the mesh
    REAL(8) :: x(1:imax), y(1:jmax)       !malha principal
    REAL(8) :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada

    REAL(8) :: vol_u(2:imax,1:jmax)       !volume de controle de u
    REAL(8) :: vol_v(1:imax,2:jmax)       !volume de controle de v
    REAL(8) :: vol_p(1:imax,1:jmax)       !volume de controle de p

    ! Áreas para u e v
    REAL(8) :: areau_n(2:imax)            !area n de u
    REAL(8) :: areau_s(2:imax)            !area s de u
    REAL(8) :: areau_e(1:jmax)            !area e de u
    REAL(8) :: areau_w(1:jmax)            !area w de u

    REAL(8) :: areav_n(1:imax)            !area n de v
    REAL(8) :: areav_s(1:imax)            !area s de v
    REAL(8) :: areav_e(2:jmax)            !area e de v
    REAL(8) :: areav_w(2:jmax)            !area w de v

    ! Variáveis auxiliares        
    REAL(8) :: dx(2:imax), dy(2:jmax)
    REAL(8) :: epsilon1(imax,jmax)
    REAL(8) :: liga_poros(imax,jmax)   
    REAL(8), PARAMETER :: rad1 = 1.d0 !raio do cilindro
    
    ! flags for obstacle interior, boundary, fluid cells, and close to the boundary
    INTEGER, PARAMETER :: C_I = 2, C_B = 1, C_F = 0, C_BS = 3     
    INTEGER :: flag(imax,jmax)

    ! Constantes físicas e parâmetros de fluidos
    REAL(8), PARAMETER :: g = 9.80665d0         !gravitational constant [m/s^2]
    REAL(8), PARAMETER :: ao  = 1.d-3           !initial radius [m]
    REAL(8), PARAMETER :: L_c = ao
    REAL(8), PARAMETER :: v_i = 0.5d0
    REAL(8), PARAMETER :: v_c = v_i    
    REAL(8), PARAMETER :: Fr = 1.0d0
    REAL(8), PARAMETER :: InvFr2 = 1.d0 / (Fr * Fr)
    REAL(8) :: S
    REAL(8) :: Lf = 1.d0
    REAL(8) :: Lo = 1.d0

    !NAMELIST ref
    REAL(8) :: Tnu    ! for Methane , 3.51d0 for n-Heptane
    REAL(8) :: YF_b 
    REAL(8) :: YO_oo 
    REAL(8) :: Ts    ! Tb  = boiling temperature [k]
    REAL(8) :: TnToo

    !Compute in main, after init. Depends of Ts
    REAL(8) :: Too   ! Too = dimen ambient temp [k]
    REAL(8) :: Tsup
    REAL(8) :: Tinf

    REAL(8), PARAMETER :: q_dim = 50.15d6 !used only to q calculation  ! q    = combustion heat [J/kg] for Methane CH4
    REAL(8) :: q !Compute in initial, equations, boundary, Depends of properties

    !Parameters computed in subroutine properties 
    REAL(8) :: cp_tot    != 1937.3540735d0! [J/kgK]   ! Affects q
    REAL(8) :: rho_tot   != 1.1950251341d0  ! [kg/m^3]! No affect
    REAL(8) :: k_tot     != 7.3404587d-2    ! [W/mK]  ! No affect
    REAL(8) :: nu_tot    != 9.30492d-5      ! [m^2/s] ! Affects Pr
    REAL(8) :: alpha_tot                    ! [m^2/s] ! Affects Pr

    REAL(8), PARAMETER :: Re = 1.d0
    REAL(8) :: Pr ! Affects Pe
    REAL(8), PARAMETER :: Pe = 1.d0
    REAL(8), PARAMETER :: Sc = 1.d0


    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um, um_n, ui
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm, vm_n, vi

    REAL(8), DIMENSION(3:imax-1,2:jmax-1) :: res_u
    REAL(8), DIMENSION(2:imax-1,3:jmax-1) :: res_v

    !Temporary variables
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: RU !,um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: RV!,vm
    !REAL(8), DIMENSION(1:imax,1:jmax) :: P
    
    real(8), DIMENSION(2:imax,1:jmax) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    REAL(8), DIMENSION(2:imax,1:jmax) :: Dn, Ds, De, Dw
    REAL(8), DIMENSION(2:imax-1,1:jmax) :: Dp
    real(8), DIMENSION(2:imax,1:jmax) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8), DIMENSION(2:imax,1:jmax) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8), DIMENSION(2:imax-1,1:jmax) :: afw, afe, afn, afs !alpha
    real(8), DIMENSION(2:imax,1:jmax) :: q_art
    real(8) :: artDivU(imax,jmax), artDivV(imax,jmax)
    real(8), DIMENSION(2:imax,1:jmax) :: dudxdx, dvdydy, dxdvdy, dydudx
    
    REAL(8), DIMENSION(2:imax-1,1:jmax) :: dudx, dvdy
    REAL(8), DIMENSION(1:imax,1:jmax) :: RP, Pi
    REAL(8), DIMENSION(2:imax-1,2:jmax-1) :: res_p

    REAL(8), DIMENSION(2:imax-1,1:jmax) :: dZudx, dZvdy
    

    REAL(8), DIMENSION(1:imax,1:jmax)     :: RZ, Zi
    REAL(8), DIMENSION(2:imax-1,2:jmax-1) :: res_Z
    
    REAL(8), DIMENSION(1:imax,  1:jmax)   :: RC, Ci
    REAL(8), DIMENSION(2:imax-1,2:jmax-1) :: res_C

    REAL(8), DIMENSION(2:imax-1,1:jmax) :: dCudx, dCvdy

end MODULE comum
