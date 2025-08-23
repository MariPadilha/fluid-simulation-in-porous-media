subroutine properties
    use comum
    implicit none
    REAL(8) :: cp_air, rho_air, k_air, nu_air, mu_air, alpha_air
    REAL(8) :: cp_met, rho_met, k_met, mu_met, nu_met, alpha_met
    REAL(8), parameter :: X_met=0.25d0, X_air=0.75d0

    rho_air = 356.6123426031d0 * Too ** (- 1.0013371777d0)
    write(*,*)'rho_air',rho_air

    cp_air = 4.73251788718275d-11 * Too ** 4.d0 - 2.5966626136914d-7 * &
             Too ** 3.d0 + 0.0004309589d0 * Too ** 2.d0 - 0.0748376208d0 * Too + 993.1365094124d0

    k_air = 4.63985599705762d-12 * Too ** 3.d0 - 0.000000029d0 * &
            Too ** 2.d0 + 9.04946015844389d-5 * Too + 0.0009128206d0

    alpha_air =  - 1.50263396909776d-14 * Too ** 3.d0 + 1.0946104054489d-10 * Too**2.d0 &
                 + 8.25165965131297d-8 * Too - 1.20726065874961d-5

    mu_air = 4.40778022994337d-15 * Too ** 3.d0 - 2.31633149512214d-11 * Too ** 2.d0 &
             + 5.78517134041232d-8 * Too + 2.98966575978449d-6

    nu_air = mu_air / rho_air

!!!!NIST webbook for Methane (C H 4 ) at 1atm

    cp_met = (3.42925052210868d-6*Too**2.d0 + 0.0003463592d0*Too + 1.8482748861d0) * 10**3.d0![J/kgK]

    rho_met =  2.60191696337177d-11*Too**4.d0 - 5.26952306479639d-8*Too**3.d0 + &
               4.12957366234466d-5*Too**2.d0 - 0.0155954244d0*Too + 2.8262044923d0 ![kg/m3]

    k_met =  1.30126895841081d-7*Too**2.d0 + 0.000064197d0*Too + 0.0037188796d0 ![W/mK]

    mu_met = 2.84900946271051d-8 * Too + 2.63131387329591d-6 ![Ns/m2 ]
    
    nu_met = mu_met/rho_met 

    alpha_met = (k_met)/(cp_met*rho_met)

    rho_tot   = X_met * rho_met + X_air * rho_air

    cp_tot    = X_met * cp_met +  X_air * cp_air

    k_tot     = X_met * k_met +   X_air * k_air
 
    nu_tot    = X_met * nu_met +  X_air * nu_air

    alpha_tot = X_met * alpha_met +  X_air * alpha_air

    !Remove it?
    ! S = (nu * YF_b) / YO_oo

    !Cf = 1.75d0 / ((150.d0 * porosidade**3.d0)**0.5d0)
    !kp = 1.d0 / 180.d0 * (porosidade**3.d0 * 100.d-6**2.0) / ((1.d0 - porosidade)**2.d0)

    !write(*,*) "Darcy number", Darcy_number

    q = (q_dim * YF_b) / (cp_tot * Too)
    Pr = nu_tot / alpha_tot

end subroutine properties
