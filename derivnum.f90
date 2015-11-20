!> \file derivnum.f90
!! \BRIEF 
!> Module with derivnum subroutine - compute numerical derivatives of carbonate system vars 
!> with respect to DIC, Alk, total phosphorus, total silicon, T, S
MODULE mderivnum
CONTAINS


!>    Computes numerical derivatives of standard carbonate system variables  
!>    (H+, pCO2, fCO2, CO2*, HCO3- and CO3--, OmegaA, OmegaC) with respect to one given input variable
!!    FROM
!!    temperature, salinity, pressure,
!!    total alkalinity (ALK), dissolved inorganic carbon (DIC),
!!    silica and phosphate concentrations
SUBROUTINE derivnum (dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx,                   &
                     dco3_dx, dOmegaA_dx, dOmegaC_dx,                                &
                     temp, sal, alk, dic, sil, phos, Patm, depth, lat, N, derivar,   &
                     optCON, optT, optP, optB, optK1K2, optKf, optGAS                )
  !   Purpose:
  !     Computes numerical derivatives of standard carbonate system variables 
  !     (H+, pCO2, fCO2, CO2*, HCO3- and CO3--, OmegaA, OmegaC) with respect to one given input variable
  !     FROM:
  !     temperature, salinity, pressure,
  !     total alkalinity (ALK), dissolved inorganic carbon (DIC),
  !     silica and phosphate concentrations

  !     INPUT variables:
  !     ================
  !     Patm    = atmospheric pressure [atm]
  !     depth   = depth [m]     (with optP='m', i.e., for a z-coordinate model vertical grid is depth, not pressure)
  !             = pressure [db] (with optP='db')
  !     lat     = latitude [degrees] (needed to convert depth to pressure, i.e., when optP='m')
  !             = dummy array (unused when optP='db')
  !     temp    = potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not in situ temp)
  !             = in situ   temperature [degrees C] (with optT='Tinsitu', e.g., for data)
  !     sal     = salinity in [psu]
  !     alk     = total alkalinity in [eq/m^3] with optCON = 'mol/m3'
  !             =               [eq/kg]  with optCON = 'mol/kg'
  !     dic     = dissolved inorganic carbon [mol/m^3] with optCON = 'mol/m3'
  !             =                            [mol/kg]  with optCON = 'mol/kg'
  !     sil     = silica    [mol/m^3] with optCON = 'mol/m3'
  !             =           [mol/kg]  with optCON = 'mol/kg'
  !     phos    = phosphate [mol/m^3] with optCON = 'mol/m3'
  !             =           [mol/kg]  with optCON = 'mol/kg'
  !     derivar = 3-character identifier of input variable with respect to which derivative is requested
  !               possibilities are 'alk', 'dic', 'pho', 'sil', 'tem', or 'sal'
  !
  !     INPUT options:
  !     ==============
  !     -----------
  !     optCON: choose input & output concentration units - mol/kg (data) vs. mol/m^3 (models)
  !     -----------
  !       -> 'mol/kg' for DIC, ALK, sil, & phos given on mokal scale, i.e., in mol/kg  (std DATA units)
  !       -> 'mol/m3' for DIC, ALK, sil, & phos given in mol/m^3 (std MODEL units)
  !     -----------
  !     optT: choose in situ vs. potential temperature as input
  !     ---------
  !     NOTE: Carbonate chem calculations require IN-SITU temperature (not potential Temperature)
  !       -> 'Tpot' means input is pot. Temperature (in situ Temp "tempis" is computed)
  !       -> 'Tinsitu' means input is already in-situ Temperature, not pot. Temp ("tempis" not computed)
  !     ---------
  !     optP: choose depth (m) vs pressure (db) as input
  !     ---------
  !       -> 'm'  means "depth" input is in "m" (thus in situ Pressure "p" [db] is computed)
  !       -> 'db' means "depth" input is already in situ pressure [db], not m (thus p = depth)
  !     ---------
  !     optB: choose total boron formulation - Uppström (1974) vs. Lee et al. (2010)
  !     ---------
  !       -> 'u74' means use classic formulation of Uppström (1974) for total Boron
  !       -> 'l10' means use newer   formulation of Lee et al. (2010) for total Boron
  !     ---------
  !     optK1K2:
  !     ---------
  !       -> 'l'   means use Lueker et al. (2000) formulations for K1 & K2 (recommended by Dickson et al. 2007)
  !                **** BUT this should only be used when 2 < T < 35 and 19 < S < 43
  !       -> 'm10' means use Millero (2010) formulation for K1 & K2 (see Dickson et al., 2007)
  !                **** Valid for 0 < T < 50°C and 1 < S < 50 psu
  !     ----------
  !     optKf:
  !     ----------
  !       -> 'pf' means use Perez & Fraga (1987) formulation for Kf (recommended by Dickson et al., 2007)
  !               **** BUT Valid for  9 < T < 33°C and 10 < S < 40.
  !       -> 'dg' means use Dickson & Riley (1979) formulation for Kf (recommended by Dickson & Goyet, 1994)
  !     -----------
  !     optGAS: choose in situ vs. potential fCO2 and pCO2
  !     ---------
  !       PRESSURE corrections for K0 and the fugacity coefficient (Cf) 
  !       -> 'Pzero'   = 'zero order' fCO2 and pCO2 (typical approach, which is flawed)
  !                      considers in situ T & only atm pressure (hydrostatic=0)
  !       -> 'Ppot'    = 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !                      considers potential T & only atm pressure (hydrostatic press = 0)
  !       -> 'Pinsitu' = 'in situ' fCO2 and pCO2 (accounts for huge effects of pressure)
  !                      considers in situ T & total pressure (atm + hydrostatic)
  !     ---------

  !     OUTPUT variables:
  !     =================
  !     dh_dx      = derivative of ion [H+] concentration on total scale in [mol/kg]
  !     dpco2_dx   = derivative of CO2 partial pressure (uatm)
  !     dfco2_dx   = derivative of CO2 fugacity (uatm)
  !     dco2_dx    = derivative of aqueous CO2 concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     dhco3_dx   = derivative of bicarbonate (HCO3-) concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     dco3_dx    = derivative of carbonate (CO3--) concentration in [mol/kg] or [mol/m^3] depending on optCON
  !     dOmegaA_dx = derivative of Omega for aragonite, i.e., the aragonite saturation state
  !     dOmegaC_dx = derivative of Omega for calcite, i.e., the   calcite saturation state


  USE msingledouble
  USE mvars
  
  IMPLICIT NONE

! Input variables
  !>     number of records
  INTEGER, INTENT(in) :: N
  !> either <b>in situ temperature</b> (when optT='Tinsitu', typical data) 
  !! OR <b>potential temperature</b> (when optT='Tpot', typical models) <b>[degree C]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: temp
  !> salinity <b>[psu]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: sal
  !> total alkalinity in <b>[eq/m^3]</b> (when optCON = 'mol/m3') OR in <b>[eq/kg]</b>  (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: alk
  !> dissolved inorganic carbon in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: dic
  !> SiO2 concentration in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: sil
  !> phosphate concentration in <b>[mol/m^3]</b> (when optCON = 'mol/m3') OR in <b>[mol/kg]</b> (when optCON = 'mol/kg')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: phos
!f2py optional , depend(sal) :: n=len(sal)
  !> atmospheric pressure <b>[atm]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: Patm
  !> depth in \b meters (when optP='m') or \b decibars (when optP='db')
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: depth
  !> latitude <b>[degrees north]</b>
  REAL(kind=rx), INTENT(in), DIMENSION(N) :: lat
  ! 3-character identifier of input variable with respect to which derivative is requested
  ! ('alk', 'dic', 'pho', 'sil', 'tem', 'sal'
  CHARACTER(3), INTENT(in) ::  derivar

  !> choose either \b 'mol/kg' (std DATA units) or \b 'mol/m3' (std MODEL units) to select 
  !! concentration units for input (for alk, dic, sil, phos) & output (co2, hco3, co3)
  CHARACTER(6), INTENT(in) :: optCON
  !> choose \b 'Tinsitu' for in situ temperature or \b 'Tpot' for potential temperature (in situ Temp is computed, needed for models)
  CHARACTER(7), INTENT(in) :: optT
  !> for depth input, choose \b "db" for decibars (in situ pressure) or \b "m" for meters (pressure is computed, needed for models)
  CHARACTER(2), INTENT(in) :: optP
  !> for total boron, choose either \b 'u74' (Uppstrom, 1974) or \b 'l10' (Lee et al., 2010).
  !! The 'l10' formulation is based on 139 measurements (instead of 20), 
  !! uses a more accurate method, and
  !! generally increases total boron in seawater by 4% 
!f2py character*3 optional, intent(in) :: optB='l10'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optB
  !> for Kf, choose either \b 'pf' (Perez & Fraga, 1987) or \b 'dg' (Dickson & Riley, 1979)
!f2py character*2 optional, intent(in) :: optKf='pf'
  CHARACTER(2), OPTIONAL, INTENT(in) :: optKf
  !> for K1,K2 choose either \b 'l' (Lueker et al., 2000) or \b 'm10' (Millero, 2010) 
!f2py character*3 optional, intent(in) :: optK1K2='l'
  CHARACTER(3), OPTIONAL, INTENT(in) :: optK1K2
  !> for K0,fugacity coefficient choose either \b 'Ppot' (no pressure correction) or \b 'Pinsitu' (with pressure correction) 
  !! 'Ppot'    - for 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
  !! 'Pinsitu' - for 'in situ' values of fCO2 and pCO2, accounting for pressure on K0 and Cf
  !! with 'Pinsitu' the fCO2 and pCO2 will be many times higher in the deep ocean
!f2py character*7 optional, intent(in) :: optGAS='Pinsitu'
  CHARACTER(7), OPTIONAL, INTENT(in) :: optGAS

! Output variables:
  !> derivative of H on the <b>total scale</b>
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dh_dx
  !> derivative of CO2 partial pressure <b>[uatm]</b>
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dpco2_dx
  !> derivative of CO2 fugacity <b>[uatm]</b>
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dfco2_dx
  !> derivative of aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dco2_dx
  !> derivative of (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dhco3_dx
  !> derivative of (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dco3_dx
  !> derivative of Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dOmegaA_dx
  !> derivative of Omega for calcite, i.e., the calcite saturation state
  REAL(kind=rx), INTENT(OUT), DIMENSION(N) :: dOmegaC_dx

! Local variables
  !> pH on the <b>total scale</b>
  REAL(kind=rx), DIMENSION(1,2) :: ph
  !> ion (H+) concentration on the <b>total scale</b>
  REAL(kind=rx), DIMENSION(1,2) :: h
  !> CO2 partial pressure <b>[uatm]</b>
  REAL(kind=rx), DIMENSION(1,2) :: pco2
  !> CO2 fugacity <b>[uatm]</b>
  REAL(kind=rx), DIMENSION(1,2) :: fco2
  !> aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optCON
  REAL(kind=rx), DIMENSION(1,2) :: co2
  !> bicarbonate ion (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), DIMENSION(1,2) :: hco3
  !> carbonate ion (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optCON
  REAL(kind=rx), DIMENSION(1,2) :: co3
  !> Omega for aragonite, i.e., the aragonite saturation state
  REAL(kind=rx), DIMENSION(1,2) :: OmegaA
  !> Omega for calcite, i.e., the calcite saturation state
  REAL(kind=rx), DIMENSION(1,2) :: OmegaC
  !> Revelle factor, i.e., dpCO2/pCO2 / dDIC/DIC
  REAL(kind=rx), DIMENSION(1) :: BetaD
  !> in-situ density of seawater; rhoSW = f(s, t, p) in <b>[kg/m3]</b>
  REAL(kind=rx), DIMENSION(1) :: rhoSW
  !> pressure <b>[decibars]</b>; p = f(depth, latitude) if computed from depth [m] (when optP='m') OR p = depth [db] (when optP='db')
  REAL(kind=rx), DIMENSION(1) :: p
  !> in-situ temperature \b <b>[degrees C]</b>
  REAL(kind=rx), DIMENSION(1) :: tempis

  ! Arrays that are copies of input parameters
  REAL(kind=rx), DIMENSION(1) :: atemp
  REAL(kind=rx), DIMENSION(1) :: asal
  REAL(kind=rx), DIMENSION(1) :: aalk
  REAL(kind=rx), DIMENSION(1) :: adic
  REAL(kind=rx), DIMENSION(1) :: asil
  REAL(kind=rx), DIMENSION(1) :: aphos
  REAL(kind=rx), DIMENSION(1) :: aPatm
  REAL(kind=rx), DIMENSION(1) :: adepth
  REAL(kind=rx), DIMENSION(1) :: alat
 
  ! value of small delta to apply to input variable when computing numerical derivative
  ! it is actually the ratio of delta relative to input variable value
  REAL(kind=rx) :: rel_delta_x
  ! Value of input variable with respect to which we derive
  REAL(kind=rx) :: input_value
  REAL(kind=rx), DIMENSION(1) :: ainput1, ainput2
  REAL(kind=rx) :: abs_delta, dX
  
! Arrays to pass optional arguments into or use defaults (Dickson et al., 2007)
  CHARACTER(3) :: opB
  CHARACTER(2) :: opKf
  CHARACTER(3) :: opK1K2
  CHARACTER(7) :: opGAS

  INTEGER :: i

! Set defaults for optional arguments (in Fortran 90)
! Note:  Optional arguments with f2py (python) are set above with 
!        the !f2py statements that precede each type declaraion
  IF (PRESENT(optB)) THEN
!   print *,"optB present:"
!   print *,"optB = ", optB 
    opB = optB
  ELSE
!   Default is Lee et al (2010) for total boron
!   print *,"optB NOT present:"
    opB = 'l10'
!   print *,"opB = ", opB 
  ENDIF
  IF (PRESENT(optKf)) THEN
!   print *,"optKf = ", optKf
    opKf = optKf
  ELSE
!   print *,"optKf NOT present:"
!   Default is Perez & Fraga (1987) for Kf
    opKf = 'pf'
!   print *,"opKf = ", opKf
  ENDIF
  IF (PRESENT(optK1K2)) THEN
!   print *,"optK1K2 = ", optK1K2
    opK1K2 = optK1K2
  ELSE
!   print *,"optK1K2 NOT present:"
!   Default is Lueker et al. 2000) for K1 & K2
    opK1K2 = 'l'
!   print *,"opK1K2 = ", opK1K2
  ENDIF
  IF (PRESENT(optGAS)) THEN
    opGAS = optGAS
  ELSE
    opGAS = 'Pinsitu'
  ENDIF

  DO i = 1, N
  
    ! Copy input values
    atemp(1)  = temp(i)
    asal(1)   = sal(i)
    aalk(1)   = alk(i)
    adic(1)   = dic(i)
    asil(1)   = sil(i)
    aphos(1)  = phos(i)
    aPatm(1)  = Patm(i)
    adepth(1) = depth(i)
    alat(1)   = lat(i)

    ! Select input value to vary slightly
    ! values for individual rel_delta_x determined by minimizing difference
    ! between numerical derivatives and mocsy-dnad automatic derivatives
    ! (mocsy_dnad), which gives results as good as analytic solutions
    SELECT CASE (derivar)
        CASE ('alk')
            input_value = alk(i)
            rel_delta_x = 1.0d-6
        CASE ('dic')
            input_value = dic(i)
            rel_delta_x = 1.0d-6
        CASE ('pho')
            input_value = phos(i)
            rel_delta_x = 1.0d-3
        CASE ('sil')
            input_value = sil(i)
            rel_delta_x = 1.0d-3
        CASE ('tem')
            input_value = temp(i)
            rel_delta_x = 1.0d-4
        CASE ('sal')
            input_value = sal(i)
            rel_delta_x = 1.0d-4
        CASE DEFAULT
            PRINT *,"derivar must be 3-char variable: 'alk', 'dic', 'pho', 'sil', 'tem', or 'sal'"
            STOP
    END SELECT
    ! Determine two slightly different values of selected input value
    abs_delta = input_value * rel_delta_x
    ainput1(1) = input_value - abs_delta
    ainput2(1) = input_value + abs_delta
    ! Compute total absolue delta
    dX = ainput2(1) - ainput1(1)
    
    ! Call routine vars() twice with two slightly different values of selected input
    SELECT CASE (derivar)
        CASE ('alk')
            call vars(ph(:,1), pco2(:,1), fco2(:,1), co2(:,1), hco3(:,1), co3(:,1),               &
                            OmegaA(:,1), OmegaC(:,1), BetaD, rhoSW, p, tempis,                    &
                            atemp, asal, ainput1, adic, asil, aphos, aPatm, adepth, alat, 1,      &
                            optCON, optT, optP, opB, opK1K2, opKf, opGAS                        ) 
            call vars(ph(:,2), pco2(:,2), fco2(:,2), co2(:,2), hco3(:,2), co3(:,2),               &
                            OmegaA(:,2), OmegaC(:,2), BetaD, rhoSW, p, tempis,                    &
                            atemp, asal, ainput2, adic, asil, aphos, aPatm, adepth, alat, 1,      &
                            optCON, optT, optP, opB, opK1K2, opKf, opGAS                        ) 
        CASE ('dic')
            call vars(ph(:,1), pco2(:,1), fco2(:,1), co2(:,1), hco3(:,1), co3(:,1),               &
                            OmegaA(:,1), OmegaC(:,1), BetaD, rhoSW, p, tempis,                    &
                            atemp, asal, aalk, ainput1, asil, aphos, aPatm, adepth, alat, 1,      &
                            optCON, optT, optP, opB, opK1K2, opKf, opGAS                        ) 
            call vars(ph(:,2), pco2(:,2), fco2(:,2), co2(:,2), hco3(:,2), co3(:,2),               &
                            OmegaA(:,2), OmegaC(:,2), BetaD, rhoSW, p, tempis,                    &
                            atemp, asal, aalk, ainput2, asil, aphos, aPatm, adepth, alat, 1,      &
                            optCON, optT, optP, opB, opK1K2, opKf, opGAS                        ) 
        CASE ('pho')
            call vars(ph(:,1), pco2(:,1), fco2(:,1), co2(:,1), hco3(:,1), co3(:,1),               &
                            OmegaA(:,1), OmegaC(:,1), BetaD, rhoSW, p, tempis,                    &
                            atemp, asal, aalk, adic, asil, ainput1, aPatm, adepth, alat, 1,       &
                            optCON, optT, optP, opB, opK1K2, opKf, opGAS                        ) 
            call vars(ph(:,2), pco2(:,2), fco2(:,2), co2(:,2), hco3(:,2), co3(:,2),               &
                            OmegaA(:,2), OmegaC(:,2), BetaD, rhoSW, p, tempis,                    &
                            atemp, asal, aalk, adic, asil, ainput2, aPatm, adepth, alat, 1,       &
                            optCON, optT, optP, opB, opK1K2, opKf, opGAS                        ) 
        CASE ('sil')
            call vars(ph(:,1), pco2(:,1), fco2(:,1), co2(:,1), hco3(:,1), co3(:,1),               &
                            OmegaA(:,1), OmegaC(:,1), BetaD, rhoSW, p, tempis,                    &
                            atemp, asal, aalk, adic, ainput1, aphos, aPatm, adepth, alat, 1,      &
                            optCON, optT, optP, opB, opK1K2, opKf, opGAS                        ) 
            call vars(ph(:,2), pco2(:,2), fco2(:,2), co2(:,2), hco3(:,2), co3(:,2),               &
                            OmegaA(:,2), OmegaC(:,2), BetaD, rhoSW, p, tempis,                    &
                            atemp, asal, aalk, adic, ainput2, aphos, aPatm, adepth, alat, 1,      &
                            optCON, optT, optP, opB, opK1K2, opKf, opGAS                        ) 
        CASE ('tem')
            call vars(ph(:,1), pco2(:,1), fco2(:,1), co2(:,1), hco3(:,1), co3(:,1),               &
                            OmegaA(:,1), OmegaC(:,1), BetaD, rhoSW, p, tempis,                    &
                            ainput1, asal, aalk, adic, asil, aphos, aPatm, adepth, alat, 1,       &
                            optCON, optT, optP, opB, opK1K2, opKf, opGAS                        ) 
            call vars(ph(:,2), pco2(:,2), fco2(:,2), co2(:,2), hco3(:,2), co3(:,2),               &
                            OmegaA(:,2), OmegaC(:,2), BetaD, rhoSW, p, tempis,                    &
                            ainput2, asal, aalk, adic, asil, aphos, aPatm, adepth, alat, 1,       &
                            optCON, optT, optP, opB, opK1K2, opKf, opGAS                        ) 
        CASE ('sal')
            call vars(ph(:,1), pco2(:,1), fco2(:,1), co2(:,1), hco3(:,1), co3(:,1),               &
                            OmegaA(:,1), OmegaC(:,1), BetaD, rhoSW, p, tempis,                    &
                            atemp, ainput1, aalk, adic, asil, aphos, aPatm, adepth, alat, 1,      &
                            optCON, optT, optP, opB, opK1K2, opKf, opGAS                        ) 
            call vars(ph(:,2), pco2(:,2), fco2(:,2), co2(:,2), hco3(:,2), co3(:,2),               &
                            OmegaA(:,2), OmegaC(:,2), BetaD, rhoSW, p, tempis,                    &
                            atemp, ainput2, aalk, adic, asil, aphos, aPatm, adepth, alat, 1,      &
                            optCON, optT, optP, opB, opK1K2, opKf, opGAS                        ) 
        CASE DEFAULT
            PRINT *,"derivar must be 3-char variable: 'alk', 'dic', 'pho', 'sil', 'tem', or 'sal'"
            STOP
    END SELECT
    
    ! H+ concentration from ph
    h(1,1) = 10**(-ph(1,1))
    h(1,2) = 10**(-ph(1,2))
    
    ! Compute derivatives by method of centered difference

    dh_dx(i)      = (h(1,2)    - h(1,1))    / dx
    dpco2_dx(i)   = (pco2(1,2) - pco2(1,1)) / dx
    dfco2_dx(i)   = (fco2(1,2) - fco2(1,1)) / dx
    dco2_dx(i)    = (co2(1,2)  - co2(1,1))  / dx
    dhco3_dx(i)   = (hco3(1,2) - hco3(1,1)) / dx
    dco3_dx(i)    = (co3(1,2)  - co3(1,1))  / dx
    dOmegaA_dx(i) = (OmegaA(1,2) - OmegaA(1,1)) / dx
    dOmegaC_dx(i) = (OmegaC(1,2) - OmegaC(1,1)) / dx
  END DO

  RETURN
END SUBROUTINE derivnum
END MODULE mderivnum
